#' Dichotimize a training expression set and fit a logistic ridge regression model which is applied to the test expression matirx.
#'
#' Dichotimize a training expression set and fit a logistic ridge regression model which is applied to the test expression matirx.
#' This function will return a set of probabilities.
#'
#' @importFrom genefilter rowttests
#' @import ridge
#' @param trainingExprData Gene expression matrix for samples for which we the phenotype is already known.
#' @param trainingPtype The known phenotype, a vector in the same order as the columns of "trainingExprData" or with the same names as colnames of "trainingExprData".
#' @param testExprData Gene expression matrix for samples on which we wish to predict a phenotype. Gene names as rows, samples names as columns.
#' @param batchCorrect The type of batch correction to be used. Options are "eb", "none", .....
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @param printOutput Set to FALSE to supress output
#' @param numGenesSelected Specifies how genes are selected for "variableSelectionMethod". Options are "tTests", "pearson" and "spearman".
#' @param numSens The number of sensitive cell lines to be fit in the logistic regression model.
#' @param numRes The number of resistant cell lines fit in the logistic regression model.
#' @param minNumSamples The minimum number of test samples, print an error if the number of columns of "testExprData" is below this threshold. A large number of test samples may be necessary to correct for batch effects.
#' @keywords internal
#' @return classifySamples
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
classifySamples <- function(trainingExprData, trainingPtype, testExprData, batchCorrect="eb", minNumSamples=10, selection=-1, printOutput=TRUE, numGenesSelected=1000, numSens=15, numRes=55)
{
  sensInd <- order(trainingPtype)[1:numSens]
  resInd <- order(trainingPtype)[(length(trainingPtype)-numRes):length(trainingPtype)]

  homDataErlot <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)

  tTests <- genefilter::rowttests(data.matrix(cbind(homDataErlot$train[, sensInd], homDataErlot$train[, resInd])), as.factor(c(rep("sens", length(sensInd)), rep("res", length(resInd)))))
  topCorVarGenes <- rownames(tTests[order(tTests[, "p.value"]),])[1:numGenesSelected] # this makes more sense than spearman of pearson....

  trainExpr <- t(homDataErlot$train[ topCorVarGenes, c(sensInd, resInd)])
  trainPtyle <- as.numeric(as.factor(c(rep("sens", length(sensInd)), rep("res", length(resInd))))) - 1

  trainDat <- data.frame(trainPtyle, trainExpr)
  cat("Fitting model, may take some time....\n")
  ridgeLogisticModel_all <- logisticRidge(trainPtyle~., data=trainDat)

  preDataLRR <- data.frame(t(homDataErlot$test[topCorVarGenes, ]))
  predsLRR <- predict(ridgeLogisticModel_all, preDataLRR)

  return(predsLRR)
}

## Functions to predict a phenotype from microarray expression data of different platforms.
##

#' Calculates phenotype from microarray data.
#'
#' This function uses ridge regression to calculate a phenotype from an gene expression,
#' given a gene expression matrix where the phenotype is already known. The function
#' also integrates the two datasets using a user-defined procedure, power transforms
#' the known phenotype and provides several other options for flexible and powerful prediction
#' from a gene expression matrix.
#'
#' @param trainingExprData The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
#' @param trainingPtype The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
#' @param testExprData The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#' @param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#' @param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent be default
#' @param minNumSamples How many training and test samples are requried. Print an error if below this threshold
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @param printOutput Set to FALSE to supress output
#' @param removeLowVaringGenesFrom what kind of genes should be removed
#' @return A vector of the estimated phenotype, in the same order as the columns of "testExprData".
#' @import sva
#' @import ridge
#' @importFrom car powerTransform
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
calcPhenotype <- function(trainingExprData, trainingPtype, testExprData, batchCorrect="eb", powerTransformPhenotype=TRUE, removeLowVaryingGenes=.2, minNumSamples=10, selection=-1, printOutput=TRUE, removeLowVaringGenesFrom="homogenizeData")
{

  # check if the supplied data are of the correct classes
  testExprData <- as.matrix(testExprData)
  if(class(testExprData)[1] != "matrix") stop("ERROR: \"testExprData\" must be a matrix.");
  if(class(trainingExprData)[1] != "matrix") stop("ERROR: \"trainingExprData\" must be a matrix.");
  if(class(trainingPtype)[1] != "numeric") stop("ERROR: \"trainingPtype\" must be a numeric vector.");
  if(ncol(trainingExprData) != length(trainingPtype)) stop("The training phenotype must be of the same length as the number of columns of the training expressin matrix.");

  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }

  # Get the homogenized data
  homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection, printOutput=printOutput)

  # Do variable selection if specified. By default we remove 20% of least varying genes from the homogenized dataset.
  # We can also remove the intersection of the lowest 20% from both training and test sets (for the gene ids remaining in the homogenized data)
  # Otherwise, keep all genes.
  # Check batchCorrect paramter
  if(!(removeLowVaringGenesFrom %in% c("homogenizeData", "rawData")))
  {
    stop("\"removeLowVaringGenesFrom\" must be one of \"homogenizeData\", \"rawData\"")
  }

  keepRows <- seq(1:nrow(homData$train)) # by default we will keep everything
  if(removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1) # if the proportion of things to keep is between 0 and 1
  {
    if(removeLowVaringGenesFrom == "homogenizeData") # if we are filtering based on the homoginized data
    {
      keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes=removeLowVaryingGenes)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if(printOutput) cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
    else if(removeLowVaringGenesFrom == "rawData") # if we are filtering based on the raw data, i.e. the intersection of the things filtered from both datasets.
    {
      evaluabeGenes <- rownames(homData$test) # pull the gene names from the genes remaining in the homoginized dataset
      keepRowsTrain <- doVariableSelection(trainingExprData[evaluabeGenes, ], removeLowVaryingGenes=removeLowVaryingGenes)
      keepRowsTest <- doVariableSelection(testExprData[evaluabeGenes, ], removeLowVaryingGenes=removeLowVaryingGenes)
      keepRows <- intersect(keepRowsTrain, keepRowsTest) # only keep things that are kept (i.e. variable) in both datasets
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if(printOutput) cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
  }

  # PowerTranform phenotype if specified.
  offset = 0
  if(powerTransformPhenotype)
  {
    if(min(trainingPtype) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
    {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }

    transForm <- car::powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }

  # create the Ridge Regression model on our training data
  if(printOutput) cat("\nFitting Ridge Regression model... ");
  trainFrame <- data.frame(Resp=trainingPtype, t(homData$train[keepRows, ]))
  rrModel <- ridge::linearRidge(Resp ~ ., data=trainFrame)
  if(printOutput) cat("Done\n\nCalculating predicted phenotype...");

  # calculate the relative contribution of each gene to the prediction
  # i might report these, I don't know if there's any point.
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta

  # predict the new phenotype for the test data.
  # if there is a single test sample, there may be problems in predicting using the predict() function for the linearRidge package
  # This "if" statement provides a workaround
  if(class(homData$test)[1] == "numeric")
  {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol=1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=rbind(testFrame, testFrame))[1]
  }
  else #predict for more than one test sample
  {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata=testFrame)
  }

  # if the response variable was transformed, untransform it now...
  if(powerTransformPhenotype)
  {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if(printOutput) cat("Done\n\n");

  return(preds)
}

#' Cross validation on training dataset
#'
#' This function does cross validation on a training set to estimate prediction accuracy on a training set.
#' If the actual test set is provided, the two datasets can be subsetted and homogenized before the
#' cross validation analysis is preformed. This may improve the estimate of prediction accuracy.
#'
#' @param trainingExprData The training data. A matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "testExprData"
#' @param trainingPtype The known phenotype for "trainingExprData". A numeric vector which MUST be the same length as the number of columns of "trainingExprData".
#' @param testExprData The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#' @param cvFold Specify the "fold" requried for cross validation. "-1" will do leave one out cross validation (LOOCV)
#' @param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#' @param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 precent be default
#' @param minNumSamples How many training and test samples are requried. Print an error if below this threshold
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @return An object of class "pRRopheticCv", which is a list with two members, "cvPtype" and "realPtype", which correspond to the cross valiation predicted phenotype and the  user provided measured phenotype respectively.
#' @import ridge
#' @importFrom sva ComBat
#' @importFrom car powerTransform
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
predictionAccuracyByCv <- function(trainingExprData, trainingPtype, testExprData=-1, cvFold=-1, powerTransformPhenotype=TRUE, batchCorrect="eb", removeLowVaryingGenes=.2, minNumSamples=10, selection=1)
{

  # check if an adequate number of training and test samples have been supplied.
  if((ncol(trainingExprData) < minNumSamples))
  {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for intrinsic difference in your training and test sets and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }

  # check if a test matrix was supplied, if not, don't homogenize the data (but store it in a list called homData anyway for convenience)
  if(is.null(testExprData))
  {
    homData <- list()
    homData$selection <- selection
    homData$train <- trainingExprData
  } else if(!is.null(testExprData)) {
    homData <- homogenizeData(testExprData, trainingExprData, batchCorrect=batchCorrect, selection=selection) # homogenize the data.
  }

  nTrain <- ncol(trainingExprData)
  predPtype <- numeric() # a numeric vector to hold the predicted phenotypes for the CV subgroups

  # Perform either N fold cross validation or LOOCV, depending on the "cvFold" variable.
  if(cvFold == -1) # if we are doing leave-one-out cross-validation (LOOCV).
  {
    for(i in 1:nTrain)
    {
      testMatTemp <- matrix(homData$train[,i], ncol=1)
      rownames(testMatTemp) <- rownames(homData$train)
      #predPtype[i] <- calcPhenotype(testMatTemp, trainCvSet[,-i], trainingPtype[-i], batchCorrect="none", minNumSamples=0, selection=homData$selection, removeLowVaryingGenes=removeLowVaryingGenes, powerTransformPhenotype=powerTransformPhenotype)
      predPtype[i] <- calcPhenotype(homData$train[,-i], trainingPtype[-i], testMatTemp, batchCorrect="none", minNumSamples=0, selection=homData$selection, removeLowVaryingGenes=removeLowVaryingGenes, powerTransformPhenotype=powerTransformPhenotype, printOutput=FALSE)

      # print an update for each 20% of the this, this is slow, so we should give some updates....
      if(i %% as.integer(nTrain/5) == 0)
        cat(paste(i, "of" , nTrain, "iterations complete. \n"))
    }
  }
  else if(cvFold > 1) # if we are doing N-fold cross validation
  {
    randTestSamplesIndex <- sample(1:nTrain) # create a random order for samples

    # create a vector which indicates which samples are in which group... and split into list of groups
    sampleGroup <- rep(cvFold, nTrain)
    groupSize <- as.integer(nTrain / cvFold)
    for(j in 1:(cvFold-1)) { sampleGroup[(((j-1)*groupSize)+1):(j*groupSize)] <- rep(j, groupSize) }
    cvGroupIndexList <- split(randTestSamplesIndex, sampleGroup)

    # predict on each of the groups....
    for(j in 1:cvFold)
    {

      # create the ranomdly chosen "training" and "test" sets for cross validation
      testCvSet <- homData$train[, cvGroupIndexList[[j]]]
      trainCvSet <- homData$train[, unlist(cvGroupIndexList[-j])]
      trainPtypeCv <- trainingPtype[unlist(cvGroupIndexList[-j])]

      predPtype <- c(predPtype, calcPhenotype(trainCvSet, trainPtypeCv, testCvSet, batchCorrect="none", minNumSamples=0, selection=homData$selection, removeLowVaryingGenes=removeLowVaryingGenes, powerTransformPhenotype=powerTransformPhenotype, printOutput=FALSE))

      cat(paste("\n", j, " of ", cvFold, " iterations complete.", sep=""))
    }

    # re-order the predicted phenotypes correctly, as they were ranomized when this started.
    predPtype <- predPtype[order(randTestSamplesIndex)]

  } else {
    stop("Unrecognised value of \"cvFold\"")
  }

  finalData <- list(cvPtype=predPtype, realPtype=trainingPtype)
  class(finalData) <- "pRRopheticCv"

  return(finalData)
}

#' R^2 from "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", i.e. the output of cross validation, calculate
#' the R^2 value for the prediction (an estimate of prediction accuracy).
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#' @param powerTranform powerTranform phenotype or not.
#' @return a numeric vector containing the R squared value from the cross validation.
#' @keywords internal
#' @importFrom car powerTransform
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
estimateRsqr.pRRopheticCv <- function(cvOutData, powerTranform=TRUE)
{
  # calculate the R^2
  return(summary(lm(cvOutData$realPtype~cvOutData$cvPtype))$r.squared)
}


#' Confidence intervals from "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", i.e. the output of cross validation, calculate
#' an average confidence interval for the predictions.
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#' @param conf the confidence interval required, by default 95 precent confidence interval.
#' @return a numeric vector containing the average upper and lower confidence interval.
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
estimateCI.pRRopheticCv <- function(cvOutData, conf=.95)
{

  # report 95% confidence intervals and 50% confidence intervals (estimated from the untransformed data)
  allDeviationsRaw <- cvOutData$cvPtype - cvOutData$realPtype
  allDeviations <- abs(allDeviationsRaw)
  inCi <- which(allDeviations < quantile(allDeviations, conf))
  ci <- c(min(allDeviationsRaw[inCi]), max(allDeviationsRaw[inCi]))
  return(ci)
}


#' Mean prediction error from "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", estiamte the mean prediction error,
#' i.e. the mean difference between the predicted and measured phenotype.
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#' @return a numeric vector containing the mean prediction error from the cross validation.
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
estimateMeanPredictionError.pRRopheticCv <- function(cvOutData)
{
  allDeviationsRaw <- abs(cvOutData$cvPtype - cvOutData$realPtype)
  return(mean(allDeviationsRaw))
}

#' Median prediction error from "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", estiamte the median prediction error,
#' i.e. the median difference between the predicted and measured phenotype.
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#' @return a numeric vector containing the median prediction error from the cross validation.
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
estimateMedianPredictionError.pRRopheticCv <- function(cvOutData)
{
  allDeviationsRaw <- abs(cvOutData$cvPtype - cvOutData$realPtype)
  return(mean(median))
}


#' Summary of "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", print various metrics that
#' summarize the performance of the cross validataion analysis
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#' @keywords internal
#' @return summary.pRRopheticCv
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
summary.pRRopheticCv <- function(cvOutData)
{
  cat("\nSummary of cross-validation results:\n\n")

  corOut <- cor.test(cvOutData[[1]], cvOutData[[2]])
  cat(paste("Pearsons correlation:", round(corOut$estimate, digits=2), ", P = ", corOut$p.value, "\n"))
  cat(paste("R-squared value: ", round(estimateRsqr.pRRopheticCv(cvOutData), digits=2), "\n", sep=""))
  cis <- estimateCI.pRRopheticCv(cvOutData)
  cat(paste("Estimated 95% confidence intervals: ", round(cis[1], digits=2), ", ", round(cis[2], digits=2), "\n", sep=""))
  cat(paste("Mean prediction error: ", round(estimateMeanPredictionError.pRRopheticCv(cvOutData), digits=2), "\n\n", sep=""))
}

#' Plot "pRRopheticCv" object.
#'
#' Given an object of class "pRRopheticCv", plot the cross validation
#' predicted values against the measured values. Also plots a regression
#' line.
#'
#' @param cvOutData an object of class "pRRopheticCv", i.e. the outpout of the "predictionAccuracyByCv()" function
#' @keywords internal
#' @return plot.pRRopheticCv
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
plot.pRRopheticCv <- function(cvOutData)
{
  coefs <- coef(lm(cvOutData$realPtype~cvOutData$cvPtype))
  plot(cvOutData$cvPtype, cvOutData$realPtype, main="Measured phenotype Vs. C.V. predicted phenotype", xlab="Predicted Phenotype", ylab="Measured Phenotype")
  abline(coefs[1], coefs[2], col="red")
}

#' a function to do variable selection on expression matrices.
#'
#' This funtino will I.e. remove genes with low variation
#' It returns a vector of row ids to keep. Note, rownames() must be specified.
#'
#' @param exprMat a matrix of expression levels, rows contain genes and columns contain samples.
#' @param removeLowVaryingGenes the proportion of low varying genes to be removed.
#' @return a vector of row ids to keep
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
doVariableSelection <- function(exprMat, removeLowVaryingGenes)
{
  vars <- apply(exprMat, 1, var)
  return(order(vars, decreasing=TRUE)[seq(1:as.integer(nrow(exprMat)*(1-removeLowVaryingGenes)))])
}

#' Take two expression matrices and return homogenized versions of the matrices.
#'
#' This function accepts two expression matrices, with gene ids as rownames() and
#' sample ids as colnames(). It will deal with duplicate gene ids.
#' subset and order the matrices correctly.
#' and perform homogenize the data using whatever method is specified (by default Combat from the sva library).
#'
#' @param testExprMat Gene expression matrix for samples on which we wish to predict a phenotype. Gene names as rows, samples names as columns.
#' @param trainExprMat Gene expression matrix for samples for which we the phenotype is already known.
#' @param batchCorrect The type of batch correction to be used. Options are "eb" for Combat, "none", or "qn" for quantile normalization.
#' @param selection parameter can be used to specify how duplicates are handled, by default value -1 means ask the user. 1 means summarize duplictes by their mean and 2 means to disguard all duplicate genes.
#' @param printOutput Set to FALSE to supress output
#' @return a list containing two entries $train and $test, which are the homogenized input matrices.
#' @importFrom sva ComBat
#' @importFrom preprocessCore normalize.quantiles
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
homogenizeData <- function(testExprMat, trainExprMat, batchCorrect="eb", selection=-1, printOutput=TRUE)
{
  # Check batchCorrect paramter
  if(!(batchCorrect %in% c("eb", "qn", "none", "rank", "rank_then_eb", "standardize")))
    stop("\"batchCorrect\" must be one of \"eb\", \"qn\", \"rank\", \"rank_then_eb\", \"standardize\" or \"none\"")

  # check if both row and column names have been specified
  if(is.null(rownames(trainExprMat)) || is.null(rownames(testExprMat)))
  {
    stop("ERROR: Gene identifiers must be specified as \"rownames()\" on both training and test expression matrices. Both matices must have the same type of gene identifiers.")
  }

  # check that some of the row names overlap between both datasets (print an error if none overlap.
  if(sum(rownames(trainExprMat) %in% rownames(testExprMat)) == 0)
  {
    stop("ERROR: The rownames() of the supplied expression matrices do not match. Note that these are case-sensitive.")
  } else {
    if(printOutput) cat(paste("\n", sum(rownames(trainExprMat) %in% rownames(testExprMat)), " gene identifiers overlap between the supplied expression matrices... \n", paste=""));
  }

  # if there are duplicate gene names, give the option of removing them or summarizing them by their mean.
  if((sum(duplicated(rownames(trainExprMat))) > 0) || sum(sum(duplicated(rownames(testExprMat))) > 0))
  {
    if(selection == -1) #print the following if we're asking the user how to proceed.
    {
      cat("\nExpression matrix contain duplicated gene identifiers (i.e. duplicate rownames()), how would you like to proceed:")
      cat("\n1. Summarize duplicated gene ids by their mean value (acceptable in most cases)")
      cat("\n2. Disguard all duplicated genes (recommended if unsure)")
      cat("\n3. Abort (if you want to deal with duplicate genes ids manually)\n")
    }

    while(is.na(selection) | selection <= 0 | selection > 3 )
    {
      selection <- readline("Selection: ")
      selection <- ifelse(grepl("[^1-3.]", selection), -1 , as.numeric(selection))
    }

    cat("\n")

    if(selection == 1) # summarize duplicates by their mean
    {
      if((sum(duplicated(rownames(trainExprMat))) > 0))
      {
        trainExprMat <- summarizeGenesByMean(trainExprMat)
      }
      if((sum(duplicated(rownames(testExprMat))) > 0))
      {
        testExprMat <- summarizeGenesByMean(testExprMat)
      }
    }
    else if(selection == 2) # disguard all duplicated genes
    {
      if((sum(duplicated(rownames(trainExprMat))) > 0))
      {
        keepGenes <- names(which(table(rownames(trainExprMat)) == 1))
        trainExprMat <- trainExprMat[keepGenes, ]
      }

      if((sum(duplicated(rownames(testExprMat))) > 0))
      {
        keepGenes <- names(which(table(rownames(testExprMat)) == 1))
        testExprMat <- testExprMat[keepGenes, ]
      }
    } else {
      stop("Exectution Aborted!")
    }
  }

  # subset and order gene ids on the expression matrices
  commonGenesIds <- rownames(trainExprMat)[rownames(trainExprMat) %in% rownames(testExprMat)]
  trainExprMat <- trainExprMat[commonGenesIds, ]
  testExprMat <- testExprMat[commonGenesIds, ]

  # subset and order the two expresison matrices
  if(batchCorrect == "eb")
  {
    # subset to common genes andbatch correct using ComBat
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame("(Intercept)"=rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)), rep("test", ncol(testExprMat))))
    combatout <- sva::ComBat(dataMat, whichbatch, mod=mod)
    return(list(train=combatout[, whichbatch=="train"], test=combatout[, whichbatch=="test"], selection=selection))
  }
  else if(batchCorrect == "standardize") # standardize to mean 0 and variance 1 in each dataset, using a non EB based approach
  {
    for(i in 1:nrow(trainExprMat))
    {
      row <- trainExprMat[i, ]
      trainExprMat[i, ] <- ((row - mean(row)) / sd(row))
    }
    for(i in 1:nrow(testExprMat))
    {
      row <- testExprMat[i, ]
      testExprMat[i, ] <- ((row - mean(row)) / sd(row))
    }
    return(list(train=trainExprMat, test=testExprMat, selection=selection))
  }
  else if(batchCorrect == "rank") # the random-rank transform approach, that may be better when applying models to RNA-seq data....
  {
    for(i in 1:nrow(trainExprMat))
    {
      trainExprMat[i, ] <- rank(trainExprMat[i, ], ties.method="random")
    }
    for(i in 1:nrow(testExprMat))
    {
      testExprMat[i, ] <- rank(testExprMat[i, ], ties.method="random")
    }
    return(list(train=trainExprMat, test=testExprMat, selection=selection))
  }
  else if(batchCorrect == "rank_then_eb") # rank-transform the RNA-seq data, then apply ComBat....
  {
    # first, rank transform the RNA-seq data
    for(i in 1:nrow(testExprMat))
    {
      testExprMat[i, ] <- rank(testExprMat[i, ], ties.method="random")
    }

    # subset to common genes andbatch correct using ComBat
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame("(Intercept)"=rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)), rep("test", ncol(testExprMat))))
    combatout <- sva::ComBat(dataMat, whichbatch, mod=mod)
    return(list(train=combatout[, whichbatch=="train"], test=combatout[, whichbatch=="test"], selection=selection))
  }
  else if(batchCorrect == "qn")
  {
    dataMat <- cbind(trainExprMat, testExprMat)
    dataMatNorm <- preprocessCore::normalize.quantiles(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)), rep("test", ncol(testExprMat))))
    return(list(train=dataMatNorm[, whichbatch=="train"], test=dataMatNorm[, whichbatch=="test"], selection=selection))
  } else {
    return(list(train=trainExprMat, test=testExprMat, selection=selection))
  }
}


## This file contains functions for prediction and classification from the CGP cell line data....

#' Given a gene expression matrix, predict drug senstivity for a drug in CGP
#'
#' Given a gene expression matrix, predict drug senstivity for a drug in CGP.
#'
#' @param testMatrix a gene expression matrix with gene names as row ids and sample names as column ids.
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param tissueType specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated). This be one any of "all" (for everything, default option), "allSolidTumors" (everything except for blood), "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive"
#' @param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#' @param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent be default
#' @param minNumSamples How many training and test samples are requried. Print an error if below this threshold
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @param printOutput Set to FALSE to supress output
#' @param removeLowVaringGenesFrom what kind of genes should be removed
#' @param dataset version of GDSC dataset
#' @return a gene expression matrix that does not contain duplicate gene ids
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
pRRopheticPredict <- function(testMatrix, drug, tissueType="all", batchCorrect="eb", powerTransformPhenotype=TRUE, removeLowVaryingGenes=.2, minNumSamples=10, selection=-1, printOutput=TRUE, removeLowVaringGenesFrom="homogenizeData", dataset="cgp2014")
{
  cgpTrainData <- getCGPinfo(drug, tissueType, dataset) # get the IC50 and expression data for this drug/tissueType

  if(!all(!duplicated(colnames(cgpTrainData$trainDataOrd)))) { # remove duplicated gene expression and IC50 in the same cell lines
    cgpTrainData$trainDataOrd <- cgpTrainData$trainDataOrd[,!duplicated(colnames(cgpTrainData$trainDataOrd))]
    cgpTrainData$ic50sOrd <- cgpTrainData$ic50sOrd[!duplicated(names(cgpTrainData$ic50sOrd))]
  }

  predictedPtype <- calcPhenotype(cgpTrainData$trainDataOrd, cgpTrainData$ic50sOrd, testMatrix, batchCorrect=batchCorrect, powerTransformPhenotype=powerTransformPhenotype, removeLowVaryingGenes=removeLowVaryingGenes, minNumSamples=minNumSamples, selection=selection, printOutput=printOutput, removeLowVaringGenesFrom=removeLowVaringGenesFrom)

  return(predictedPtype)

}


#' This function uses X fold cross validation on the TrainingSet to estimate the accuracy of the
#' phenotype prediction fold: How many fold cross-validation to use.
#'
#' This function does cross validation on a training set to estimate prediction accuracy on a training set.
#' If the actual test set is provided, the two datasets can be subsetted and homogenized before the
#' cross validation analysis is preformed. This may improve the estimate of prediction accuracy.
#'
#' @param testExprData The test data where the phenotype will be estimted. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#' @param tissueType specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated). This be one any of "all" (for everything, default option), "allSolidTumors" (everything except for blood), "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive"
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param cvFold Specify the "fold" requried for cross validation. "-1" will do leave one out cross validation (LOOCV)
#' @param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#' @param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent by default.
#' @param minNumSamples How many training and test samples are requried. Print an error if below this threshold
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @return An object of class "pRRopheticCv", which is a list with two members, "cvPtype" and "realPtype", which correspond to the cross valiation predicted phenotype and the  user provided measured phenotype respectively.
#' @import ridge
#' @importFrom sva ComBat
#' @importFrom car powerTransform
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
pRRopheticCV <- function(drug, tissueType="all", testExprData=NULL, cvFold=-1, powerTransformPhenotype=TRUE, batchCorrect="eb", removeLowVaryingGenes=.2, minNumSamples=10, selection=1)
{
  cgpTrainData <- getCGPinfo(drug, tissueType) # get the IC50 and expression data for this drug/tissueType

  # I may need to alter this function so it can either take the test data or not take the test data...
  cvOut <- predictionAccuracyByCv(cgpTrainData$trainDataOrd, cgpTrainData$ic50sOrd, testExprData=testExprData, cvFold=cvFold, powerTransformPhenotype=powerTransformPhenotype, batchCorrect=batchCorrect, removeLowVaryingGenes=removeLowVaryingGenes, minNumSamples=minNumSamples, selection=selection)
  return(cvOut)
}

#' Given a drug and tissue type, return CGP expression and drug sensitivity data.
#'
#' Given a drug and tissue type, return CGP expression and drug sensitivity data.
#'
#' @param drug The name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param tissueType Specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated). This be one any of "all" (for everything, default option), "allSolidTumors" (everything except for blood), "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive"
#' @param dataset The dataset from which you wish to build the predictive models. Default is "cgp2012", also available "cgp2016", comming soon "ctrp".
#' @return a list with two entries, trainDataOrd the ordered expression data and ic50sOrd the drug sensitivity data.
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
getCGPinfo <-  function(drug, tissueType="all", dataset="cgp2014")
{
  # list of possible datasets that can be used.
  possibleDatasets <- c("cgp2014", "cgp2016")
  if(!dataset %in% possibleDatasets) stop(paste("ERROR: the dataset specified was not found. Note dataset names are case sensitive. Please select from: ", (paste(possibleDatasets, collapse=", "))));

  if(dataset == "cgp2014")
  {
    # was a valid tissue type specified; tissue types represeted by > 40 cell lines
    if(!tissueType %in% c("all", "allSolidTumors", "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive")) stop("ERROR: the tissue type specified must be one of \"all\", \"allSolidTumors\", \"blood\", \"breast\", \"CNS\", \"GI tract\", \"lung\", \"skin\" or \"upper aerodigestive\". These tissue types are represented by at least 40 cell lines (although not all 40 may have been screened with each drug).");

    # was a valid drug specified
    possibleDrugs <- c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")
    if(!drug %in% possibleDrugs) stop(paste("ERROR: the drug specified was not found. Note drug names are case sensitive. Please select from: ", (paste(possibleDrugs, collapse=", "))));

    # load("/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/data/drugAndPhenoCgp.RData") # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)
    data("drugAndPhenoCgp") # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)

    colIc50Name <- paste(drug, "_IC_50", sep="")
    ic50s <- as.numeric(drugSensitivityDataCgp[, colIc50Name])
    names(ic50s) <- drugSensitivityDataCgp[ ,"Cell.Line"]
    whichNas <- which(is.na(ic50s))
    ic50s <- ic50s[-whichNas]
    tissue <- drugSensitivityDataCgp[ ,"Cancer.Type"]
    names(tissue) <- drugSensitivityDataCgp[ ,"Cell.Line"]
    tissue <- tissue[-whichNas]

    # if a tissue type has been specified, use only tissues of that type.
    if(tissueType != "all")
    {
      if(tissueType == "allSolidTumors")
      {
        tissueType <- ic50s <- ic50s[!(tissue %in% "blood")]
      }
      else
      {
        ic50s <- ic50s[tissue %in% tissueType]
      }
    }

    # map the drug sensitivity and expression data
    pDataUnique <- drugToCellLineDataCgp[drugToCellLineDataCgp$Source.Name %in%
                                           names(which(table(drugToCellLineDataCgp$Source.Name) == 1)), ]
    rownames(pDataUnique) <- pDataUnique$Source.Name
    commonCellLines <- rownames(pDataUnique)[rownames(pDataUnique) %in% names(ic50s)]
    pDataUniqueOrd <- pDataUnique[commonCellLines, ]
    ic50sOrd <- ic50s[commonCellLines]
    trainDataOrd <- gdsc_brainarray_syms[, pDataUniqueOrd$"Array.Data.File"]

    return(list(ic50sOrd=ic50sOrd, trainDataOrd=trainDataOrd))
  }
  else if(dataset == "cgp2016")
  {
    cat("\nUsing updated CGP 2016 datsets for prediction\n\n")

    # was a valid tissue type specified; tissue types represeted by > 40 cell lines
    if(!tissueType %in% c("all", "aero_digestive_tract", "blood", "bone", "breast", "digestive_system", "lung", "nervous_system", "skin", "urogenital_system")) stop("ERROR: the tissue type specified must be one of \"all\", \"aero_digestive_tract\", \"blood\", \"bone\", \"breast\", \"digestive_system\", \"lung\", \"nervous_system\", \"skin\", \"urogenital_system\". These tissue types are represented by at least 40 cell lines (although not all 40 may have been screened with each drug).");

    # was a valid drug specified
    possibleDrugs2016 <- c("Erlotinib", "Rapamycin", "Sunitinib", "PHA-665752", "MG-132", "Paclitaxel", "Cyclopamine", "AZ628", "Sorafenib", "VX-680", "Imatinib", "TAE684", "Crizotinib", "Saracatinib", "S-Trityl-L-cysteine", "Z-LLNle-CHO", "Dasatinib", "GNF-2", "CGP-60474", "CGP-082996", "A-770041", "WH-4-023", "WZ-1-84", "BI-2536", "BMS-536924", "BMS-509744", "CMK", "Pyrimethamine", "JW-7-52-1", "A-443654", "GW843682X", "MS-275", "Parthenolide", "KIN001-135", "TGX221", "Bortezomib", "XMD8-85", "Roscovitine", "Salubrinal", "Lapatinib", "GSK269962A", "Doxorubicin", "Etoposide", "Gemcitabine", "Mitomycin C", "Vinorelbine", "NSC-87877", "Bicalutamide", "QS11", "CP466722", "Midostaurin", "CHIR-99021", "AP-24534", "AZD6482", "JNK-9L", "PF-562271", "HG-6-64-1", "JQ1", "JQ12", "DMOG", "FTI-277", "OSU-03012", "Shikonin", "AKT inhibitor VIII", "Embelin", "FH535", "PAC-1", "IPA-3", "GSK-650394", "BAY 61-3606", "5-Fluorouracil", "Thapsigargin", "Obatoclax Mesylate", "BMS-754807", "Lisitinib", "Bexarotene", "Bleomycin", "LFM-A13", "GW-2580", "AUY922", "Phenformin", "Bryostatin 1", "Pazopanib", "LAQ824", "Epothilone B", "GSK1904529A", "BMS345541", "Tipifarnib", "BMS-708163", "Ruxolitinib", "AS601245", "Ispinesib Mesylate", "TL-2-105", "AT-7519", "TAK-715", "BX-912", "ZSTK474", "AS605240", "Genentech Cpd 10", "GSK1070916", "KIN001-102", "LY317615", "GSK429286A", "FMK", "QL-XII-47", "CAL-101", "UNC0638", "XL-184", "WZ3105", "XMD14-99", "AC220", "CP724714", "JW-7-24-1", "NPK76-II-72-1", "STF-62247", "NG-25", "TL-1-85", "VX-11e", "FR-180204", "Tubastatin A", "Zibotentan", "YM155", "NSC-207895", "VNLG/124", "AR-42", "CUDC-101", "Belinostat", "I-BET-762", "CAY10603", "Linifanib ", "BIX02189", "CH5424802", "EKB-569", "GSK2126458", "KIN001-236", "KIN001-244", "KIN001-055", "KIN001-260", "KIN001-266", "Masitinib", "MP470", "MPS-1-IN-1", "BHG712", "OSI-930", "OSI-027", "CX-5461", "PHA-793887", "PI-103", "PIK-93", "SB52334", "TPCA-1", "TG101348", "Foretinib", "Y-39983", "YM201636", "Tivozanib", "GSK690693", "SNX-2112", "QL-XI-92", "XMD13-2", "QL-X-138", "XMD15-27", "T0901317", "EX-527", "THZ-2-49", "KIN001-270", "THZ-2-102-1", "AICAR", "Camptothecin", "Vinblastine", "Cisplatin", "Cytarabine", "Docetaxel", "Methotrexate", "ATRA", "Gefitinib", "Navitoclax", "Vorinostat", "Nilotinib", "RDEA119", "CI-1040", "Temsirolimus", "Olaparib", "Veliparib", "Bosutinib", "Lenalidomide", "Axitinib", "AZD7762", "GW 441756", "CEP-701", "SB 216763", "17-AAG", "VX-702", "AMG-706", "KU-55933", "Elesclomol", "Afatinib", "GDC0449", "PLX4720", "BX-795", "NU-7441", "SL 0101-1", "BIRB 0796", "JNK Inhibitor VIII", "681640", "Nutlin-3a (-)", "PD-173074", "ZM-447439", "RO-3306", "MK-2206", "PD-0332991", "BEZ235", "GDC0941", "AZD8055", "PD-0325901", "SB590885", "selumetinib", "CCT007093", "EHT 1864", "Cetuximab", "PF-4708671", "JNJ-26854165", "HG-5-113-01", "HG-5-88-01", "TW 37", "XMD11-85h", "ZG-10", "XMD8-92", "QL-VIII-58", "CCT018159", "AG-014699", "SB 505124", "Tamoxifen", "QL-XII-61", "PFI-1", "IOX2", "YK 4-279", "(5Z)-7-Oxozeaenol", "piperlongumine", "FK866", "Talazoparib", "rTRAIL", "UNC1215", "SGC0946", "XAV939", "Trametinib", "Dabrafenib", "Temozolomide", "Bleomycin (50 uM)", "SN-38", "MLN4924")
    if(!drug %in% possibleDrugs2016) stop(paste("ERROR: the drug specified was not found. Note drug names are case sensitive. Please select from: ", (paste(possibleDrugs2016, collapse=", "))));

    # load("/home/pgeeleher/Dropbox/HDAC_project_Scripts/r_package_files/pRRophetic/data/PANCANCER_IC_Tue_Aug_9_15_28_57_2016.RData") #
    data("PANCANCER_IC_Tue_Aug_9_15_28_57_2016") # contains "drugData2016" the 2016 drug IC50 data, downloaded from (http://www.cancerrxgene.org/translation/drug/download#ic50)

    # load("/home/pgeeleher/Dropbox/HDAC_project_Scripts/r_package_files/pRRophetic/data/cgp2016ExprRma.RData") #
    data("cgp2016ExprRma") # contains "cgp2016ExprRma" the 2016 gene expression data. Data was obtained from "http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip". Cosmic Ids (in the column names) were mapped to cell line names using data from this file: ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Cell_Lines_Details.xlsx

    # get the IC50s and tissue types for cell lines screened with this drug.
    ic50s <-  drugData2016[, "IC50"][drugData2016[, "Drug.name"] == drug]
    names(ic50s) <- drugData2016[ ,"Cell.line.name"][drugData2016[, "Drug.name"] == drug]
    tissue <- drugData2016[, "Tissue"][drugData2016[, "Drug.name"] == drug]
    names(tissue) <- drugData2016[ ,"Cell.line.name"][drugData2016[, "Drug.name"] == drug]

    # if a tissue type has been specified, use only tissues of that type.
    if(tissueType != "all")
    {
      if(tissueType == "allSolidTumors")
      {
        tissueType <- ic50s <- ic50s[!(tissue %in% "blood")]
      }
      else
      {
        ic50s <- ic50s[tissue %in% tissueType]
      }
    }

    # get the ordered subsetted expression and IC50 data
    commonCellLines <- names(ic50s)[names(ic50s) %in% colnames(cgp2016ExprRma)]
    ic50sOrd <- ic50s[commonCellLines]
    trainDataOrd <- cgp2016ExprRma[, commonCellLines]

    return(list(ic50sOrd=ic50sOrd, trainDataOrd=trainDataOrd))
  }
}

#' Predict from the CGP data using a logistic model
#'
#' Predict from the CGP data using a logistic model.
#'
#' @param testMatrix a gene expression matrix with gene names as row ids and sample names as column ids.
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param tissueType specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated). This be one any of "all" (for everything, default option), "allSolidTumors" (everything except for blood), "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive"
#' @param batchCorrect The type of batch correction to be used. Options are "eb", "none", .....
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @param printOutput Set to FALSE to supress output
#' @param numGenesSelected Specifies how genes are selected for "variableSelectionMethod". Options are "tTests", "pearson" and "spearman".
#' @param numSens The number of sensitive cell lines to be fit in the logistic regression model.
#' @param numRes The number of resistant cell lines fit in the logistic regression model.
#' @param minNumSamples The minimum number of test samples, print an error if the number of columns of "testExprData" is below this threshold. A large number of test samples may be necessary to correct for batch effects.
#' @return A predicted probability of sensitive or resistant from the logistic regression model.
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
pRRopheticLogisticPredict <- function(testMatrix, drug, tissueType="all", batchCorrect="eb", minNumSamples=10, selection=-1, printOutput=TRUE, numGenesSelected=1000, numSens=15, numRes=55)
{
  cgpTrainData <- getCGPinfo(drug, tissueType) # get the IC50 and expression data for this drug/tissueType

  predictedPtype <- classifySamples(cgpTrainData$trainDataOrd, cgpTrainData$ic50sOrd, testMatrix, batchCorrect=batchCorrect,minNumSamples=minNumSamples, selection=selection, printOutput=printOutput, numGenesSelected=numGenesSelected, numSens=numSens, numRes=numRes)

  return(predictedPtype[,1])
}


#' Check the distribution of the drug response (IC50) data using a QQ-plot.
#'
#' Visualize the distribution of the transformed IC50 data for a drug of interest using a QQ plot. If the distribution of the IC50 values deviates wildly from a normal distribtion, it is likely not suitalbe for prediction using a linear model (like linear ridge regression). This drug may be more suitable to constructing a model using a logistic or other type of model.
#
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439
#' @importFrom car powerTransform
#' @keywords internal
#' @return pRRopheticQQplot
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
pRRopheticQQplot <- function(drug)
{
  possibleDrugs <- c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")
  if(!drug %in% possibleDrugs) stop(paste("ERROR: the drug specified was not found. Note drug names are case sensitive. Please select from: ", (paste(possibleDrugs, collapse=", "))));

  # load("/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/data/drugAndPhenoCgp.RData") # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)
  data("drugAndPhenoCgp") # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)

  colIc50Name <- paste(drug, "_IC_50", sep="")
  ic50s <- as.numeric(drugSensitivityDataCgp[, colIc50Name])
  names(ic50s) <- drugSensitivityDataCgp[ ,"Cell.Line"]
  whichNas <- which(is.na(ic50s))
  ic50s <- ic50s[-whichNas]

  offset = 0
  if(min(ic50s) < 0) # all numbers must be postive for a powerTranform to work, so make them positive.
  {
    offset <- -min(ic50s) + 1
    ic50s <- ic50s + offset
  }

  transForm <- powerTransform(ic50s)[[6]]
  ic50s <- ic50s^transForm

  qqnorm(ic50s, main=paste("QQplot on power-transformed IC50 values for", drug))
  qqline(ic50s, col="red")
}


#' Calculate PPV and NPV and a cutpoint from the training data.
#'
#' Calculate PPV and NPV and a cutpoint from the training data.
#'
#' @param predResponders a numeric vector of the predictions that were obtained for the known drug responders
#' @param predNonResponders a numeric vector of the predictions for the known drug non-responders
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param tissueType specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated). This be one any of "all" (for everything, default option), "allSolidTumors" (everything except for blood), "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive"#'
#' @return a list with two entries, trainDataOrd the ordered expression data and ic50sOrd the drug sensitivity data.
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
getPPV <- function(predResponders, predNonResponders, drug, tissueType="all")
{
  possibleDrugs <- c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")
  if(!drug %in% possibleDrugs) stop(paste("ERROR: the drug specified was not found. Note drug names are case sensitive. Please select from: ", (paste(possibleDrugs, collapse=", "))));

  # load("/home/pgeeleher/postdoc_stuff/HDAC_project/Scripts/r_package_files/data/drugAndPhenoCgp.RData") # contains drugToCellLineDataCgp (to map cell lines to .CEL file names), gdsc_brainarray_syms (the gene expression data), drugSensitivityDataCgp (the drug sensitivity data)
  cutpoint <- mean(getCGPinfo(drug, tissueType)$ic50sOrd)

  # positive predictive values is the number of true positive / the total number of positive calls.
  ppv <- sum(predResponders < cutpoint) / (sum(predResponders < cutpoint) + sum(predNonResponders < cutpoint))

  # negative predictive values is the number of true negatives / the total number of negative calls.
  npv <- sum(predNonResponders > cutpoint) / (sum(predResponders > cutpoint) + sum(predNonResponders > cutpoint))

  cat(paste("\nPPV: ", round(ppv, 2), "\nNPV: ", round(npv, 2), "\nCutpoint: ", round(cutpoint, 2), "\n\n"))

  return(list(ppv=ppv, npv=npv, cutpoint=cutpoint))
}

#' Take an expression matrix and if duplicate genes are measured, summarize them by their means
#'
#' This function accepts two expression matrices, with gene ids as rownames() and
#' sample ids as colnames(). It will find all duplicate genes and summarize their
#' expression by their mean.
#'
#' @param exprMat a gene expression matrix with gene names as row ids and sample names as column ids.
#' @return a gene expression matrix that does not contain duplicate gene ids
#' @keywords internal
#' @author Paul Geeleher, Nancy Cox, R. Stephanie Huang
summarizeGenesByMean <- function(exprMat)
{
  geneIds <- rownames(exprMat)
  t <- table(geneIds) # how many times is each gene name duplicated
  allNumDups <- unique(t)
  allNumDups <- allNumDups[-which(allNumDups == 1)]

  # create a *new* gene expression matrix with everything in the correct order....
  # start by just adding stuff that isn't duplicated
  exprMatUnique <- exprMat[which(geneIds %in% names(t[t == 1])), ]
  gnamesUnique <- geneIds[which(geneIds %in% names(t[t == 1]))]

  # add all the duplicated genes to the bottom of "exprMatUniqueHuman", summarizing as you go
  for(numDups in allNumDups)
  {
    geneList <- names(which(t == numDups))

    for(i in 1:length(geneList))
    {
      exprMatUnique <- rbind(exprMatUnique, colMeans(exprMat[which(geneIds == geneList[i]), ]))
      gnamesUnique <- c(gnamesUnique, geneList[i])
      # print(i)
    }
  }

  if(class(exprMatUnique) == "numeric")
  {
    exprMatUnique <- matrix(exprMatUnique, ncol=1)
  }

  rownames(exprMatUnique) <- gnamesUnique
  return(exprMatUnique)
}
