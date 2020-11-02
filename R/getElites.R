#' @name getElites
#' @title Get elites for clustering
#' @description This function provides several methods to help selecting elites from input features, which aims to reduce data dimension for multi-omics integrative clustering analysis.
#' @param dat A data.frame of one omics data, can be continuous or binary data.
#' @param surv.info A data.frame with rownames of observations and with at least two columns of `futime` for survival time and `fustat` for survival status (0: censoring; 1: event)
#' @param method A string value to indicate the filtering method for selecting elites. Allowed values contain c('mad', 'sd', 'pca', 'cox', 'freq'). 'mad' means median absolute deviation, 'sd' means standard deviation, 'pca' means principal components analysis, 'cox' means univariate Cox proportional hazards regression which needs surv.info also, 'freq' only works for binary data; "mad" by default.
#' @param na.action A string value to indicate the action for handling NA missing value. Allowed values contain c('rm', 'impute'). 'rm' means removal of all features containing any missing values, 'impute' means imputation for missing values by k-nearest neighbors; "rm" by default.
#' @param doLog2 A logic value to indicate if performing log2 transformation for data before calculating statistics (e.g., sd, mad , pca and cox). FALSE by default.
#' @param lowpct A numeric cutoff for removing low expression values. NULL by default; 0.1 is recommended for continuous data which means features that have no expression in more than 10\% samples will be removed. Otherwise default value of NULL should be kept for binary data.
#' @param p.cutoff A numeric cutoff for nominal p value derived from univariate Cox proportional hazards regression; 0.05 by default.
#' @param elite.pct A numeric cutoff of percentage for selecting elites. NOTE: epite.pct works for all methods except for 'cox', but two scenarios exist. 1) when using method of 'mad' or 'sd', features will be descending sorted by mad or sd, and top elites.pct \* feature size of elites (features) will be selected; 2) when using method of 'freq' for binary data, frequency for value of 1 will be calculated for each feature, and features that have value of 1 in greater than elites.pct \* sample size will be considered elites. This argument will be discarded if elite.num is provided simultaneously. Set this argument with 1 and leave elite.num NULL will return all the features as elites after dealing with NA values.
#' @param elite.num A integer cutoff of exact number for selecting elites. NOTE: elite.num works for all methods except for 'cox', but two scenarios exist. 1) when using method of 'mad' or 'sd', features will be descending sorted by mad or sd, and top elite.num of elites (features) will be selected; 2) when using method of 'freq' for binary data, frequency for value of 1 will be calculated for each feature, and features that have value of 1 in greater than elite.num of sample size will be considered elites.
#' @param pca.ratio A numeric value ranging from 0 to 1 which represents the ratio of principal components is selected; 0.9 by default.
#' @param scaleFlag A logic value to indicate if scaling the data after filtering. FALSE by default.
#' @param centerFlag A logic value to indicate if centering the data after filtering. FALSE by default.
#' @import survival
#' @importFrom impute impute.knn
#' @return A list containing the following components:
#'
#'         \code{elite.dat} a data.frame containing data for selected elites (features).
#'
#'         \code{pca.res} a data.frame containing results for principal components analysis if \code{method == 'pca'}
#'
#'         \code{unicox.res} a data.frame containing results for univariate Cox proportional hazards regression if \code{method == 'cox'}
#' @export
#' @examples # There is no example and please refer to vignette.
getElites <- function(dat        = NULL,
                      surv.info  = NULL,
                      method     = "mad",
                      na.action  = "rm",
                      doLog2     = FALSE,
                      lowpct     = NULL,
                      p.cutoff   = 0.05,
                      elite.pct  = NULL,
                      elite.num  = NULL,
                      pca.ratio  = 0.9,
                      scaleFlag  = FALSE,
                      centerFlag = FALSE) {
  # check argument
  if(!is.element(method, c("mad","sd","pca","cox","freq"))) {
    stop("unsupportive method to filter elite features, allowed values contain c('mad','sd','pca','cox','freq').")}
  if(!is.element(na.action, c("rm","impute"))) {
    stop("unsupportive method to handle NA missing value, allowed values contain c('rm','impute').")}

  # deal with missing value
  if(na.action == "rm") {
    df <- as.data.frame(na.omit(dat))
    if(nrow(df) != nrow(dat)) {message(paste0("--",(nrow(dat) - nrow(df)), " features with NA values are removed."))}
  } else if(na.action == "impute") {
    df <- as.data.frame(impute::impute.knn(as.matrix(dat))$data)
  } else {
    df <- as.data.frame(dat)
  }

  if(doLog2) {df <- as.data.frame(log2(df + 1))}

  # remove low expression
  if(!is.null(lowpct)) {
    n.raw <- nrow(df)
    PASSFlag <- rowSums(df == 0) < ceiling(lowpct * ncol(df))
    names(PASSFlag) <- rownames(df)
    df <- df[PASSFlag,]
    if(nrow(df) != n.raw) {
      message(paste0("--remove ",n.raw-nrow(df)," features with 0 values in more than ", lowpct*100, "% samples."))
    }
  }

  # select elite features
  if(method == "freq") {
    message("--method of 'freq' only supports binary omics data (e.g., somatic mutation matrix), and in this manner, elite.pct and elite.num are used to cut frequency.")
    if(sum(unique(as.vector(as.matrix(df)))) == 1) {
      #message("--it seems the input data correctely contains binary values of 0 and 1 only.")
      statistic <- rowSums(df); names(statistic) <- rownames(df)
    } else {
      stop("it seems the input data contains more values than 0 and 1, please check and choose another proper method!")
    }
  }
  if(method != "freq" & method != "cox" & sum(unique(as.vector(as.matrix(df)))) == 1) {
    stop("it seems the input data is a binary matrix with entries of 0 and 1 only, please check and choose 'freq' or 'cox' methods if appropriate!")
  }

  if(method == "mad") {
    statistic <- apply(df, 1, mad)
    names(statistic) <- rownames(df)
    statistic <- sort(statistic, decreasing = TRUE)
  }

  if(method == "sd") {
    statistic <- apply(df, 1, sd)
    names(statistic) <- rownames(df)
    statistic <- sort(statistic, decreasing = TRUE)
  }

  if(method == "cox" & !is.null(surv.info)) {
    if(all(is.element(c("futime", "fustat"), colnames(surv.info)))) {

      display.progress = function (index, totalN, breakN=20) {
        if ( index %% ceiling(totalN/breakN)  == 0  ) {
          cat(paste(round(index*100/totalN), "% ", sep = ""))
        }
      }

      comsam <- intersect(colnames(df), rownames(surv.info))
      df <- df[,comsam]; surv.info <- surv.info[comsam,,drop = FALSE]
      if(length(comsam) == ncol(df)) {
        message("--all sample matched between omics matrix and survival data.")
      } else {
        message(paste0("--removed ",(ncol(df) - length(comsam))," samples that mismatched between omics matrix and survival data."))
      }

      unicox <- data.frame()
      for(i in 1:nrow(df)){

        display.progress(index = i, totalN = nrow(df))
        gene <- rownames(df)[i]
        tmp <- data.frame(expr = as.numeric(df[i,]),
                          futime = surv.info$futime,
                          fustat = surv.info$fustat,
                          stringsAsFactors = FALSE)
        cox <- survival::coxph(survival::Surv(futime, fustat) ~ expr, data = tmp)
        coxSummary <- summary(cox)
        unicox <- rbind.data.frame(unicox,
                                   data.frame(gene             = gene,
                                              HR               = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                              z                = as.numeric(coxSummary$coefficients[,"z"])[1],
                                              pvalue           = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                              lower            = as.numeric(coxSummary$conf.int[,3][1]),
                                              upper            = as.numeric(coxSummary$conf.int[,4][1]),
                                              stringsAsFactors = FALSE),
                                   stringsAsFactors            = FALSE)
      }
      statistic <- unicox$pvalue
      names(statistic) <- unicox$gene
    } else {stop("survival information must contain columns of 'futime' for survival time and 'fustat' for survival status.")}
  }

  if(method == "freq") {
    if(!is.null(elite.num) & !is.null(elite.pct)) {
      message("elite.num has been provided then discards elite.pct.")
      elite <- names(statistic[statistic > elite.num])
    }
    if(!is.null(elite.num) & is.null(elite.pct)) {
      message("missing elite.pct then use elite.num.")
      elite <- names(statistic[statistic > elite.num])
    }
    if(is.null(elite.num) & !is.null(elite.pct)) {
      message("missing elite.num then use elite.pct")
      elite <- names(statistic[statistic > elite.pct * ncol(df)])
    }
  }

  if(method == "cox") {
    elite <- names(statistic[statistic < p.cutoff])
  }
  if(method %in% c("sd", "mad")) {
    if(!is.null(elite.num) & !is.null(elite.pct)) {
      message("elite.num has been provided then discards elite.pct.")
      if(elite.num > length(statistic)) {
        stop(paste0("number of remaining features (n=", length(statistic), ") less than elite.num of ",elite.num,", please consider to use elite.pct."))
      } else {
        elite <- names(statistic)[1:elite.num]
      }
    }
    if(!is.null(elite.num) & is.null(elite.pct)) {
      message("missing elite.pct then use elite.num.")
      if(elite.num > length(statistic)) {
        stop(paste0("number of remaining features (n=", length(statistic), ") less than elite.num of ",elite.num,", please consider to use elite.pct."))
      } else {
        elite <- names(statistic)[1:elite.num]
      }
    }
    if(is.null(elite.num) & !is.null(elite.pct)) {
      message("missing elite.num then use elite.pct")
      elite <- names(statistic)[1:(elite.pct * nrow(df))]
    }
  }

  if(method == "pca") {
    message(paste0("--the ratio used to select principal component is set as ", pca.ratio))
    tdf <- t(df)
    pca.res <- prcomp(tdf, scale = T, center = T)
    vars <- apply(pca.res$x, 2, var)
    props <- vars / sum(vars)
    cprops <- as.vector(cumsum(props))
    outdat <- (tdf %*% pca.res$ro)[,1:which(cprops > pca.ratio)[1]]
    outdat <- as.data.frame(t(scale(outdat, center = centerFlag, scale = scaleFlag)))
  } else {
    outdat <- as.data.frame(t(scale(t(df[elite, , drop = FALSE]), center = centerFlag, scale = scaleFlag)))
  }

  if(method == "cox") {
    return(list(elite.dat = outdat, unicox.res = unicox))
  } else if(method == "pca") {
    return(list(elite.dat = outdat, pca.res = pca.res))
  } else {return(list(elite.dat = outdat))}
}
