#' @name runPAM
#' @title Run partition around medoids classifier
#' @description Using partition around medoids (PAM) classifier to predict potential subtype label on external cohort and calculate in-group proportions (IGP) statistics.
#' @param train.expr A matrix of normalized expression training data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param test.expr A matrix of normalized expression testing data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.
#' @param gene.subset A string vector to indicate a subset of genes to be used.
#' @return A list with the following components:
#'
#'         \code{IGP}        a named numeric vector storing the in-group proportion (see \link[clusterRepro]{IGP.clusterRepro}).
#'
#'         \code{clust.res}  similar to `clust.res` returned by `getMOIC()` or `get%algorithm_name%` or `getConsensusMOIC()`.
#'
#'         \code{mo.method}  a string value indicating the method used for prediction.
#' @details This function first trains a partition around medoids (PAM) classifier in the discovery (training) cohort
#'  to predict the subtype for patients in the external validation (testing) cohort,
#'  and each sample in the validation cohort was assigned to a subtype label whose centroid had the highest Pearson correlation with the sample.
#'  Finally, the in-group proportion (IGP) statistic will be performed to evaluate the similarity and reproducibility of the acquired subtypes between discovery and validation cohorts.
#' @export
#' @importFrom pamr pamr.train
#' @importFrom clusterRepro IGP.clusterRepro
#' @importFrom grDevices colorRamp
#' @references Tibshirani R, Hastie T, Narasimhan B and Chu G (2002). Diagnosis of multiple cancer types by shrunken centroids of gene expression. Proc Natl Acad Sci, 99,6567–6572.
#'
#' Kapp A V, Tibshirani R. (2007). Are clusters found in one dataset present in another dataset?. Biostatistics, 8(1):9-31.
#' @examples # There is no example and please refer to vignette.
runPAM <- function(train.expr  = NULL,
                   moic.res    = NULL,
                   test.expr   = NULL,
                   gene.subset = NULL) {

  # check sample
  comsam <- intersect(moic.res$clust.res$samID, colnames(train.expr))
  if(length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(moic.res$clust.res)-length(comsam))," samples mismatched from current subtypes."))
  }

  moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  train.expr <- train.expr[,comsam]

  # check if using subset of genes
  if(is.null(gene.subset)) {
    comgene <- intersect(rownames(train.expr), rownames(test.expr))
    message(paste0("--a total of ",length(comgene)," genes shared and used."))

  } else {
    comgene <- intersect(intersect(rownames(train.expr), rownames(test.expr)), gene.subset)
    message(paste0("--a total of ",length(comgene)," genes shared in the gene subset and used."))
  }

  train.expr <- train.expr[comgene,]
  test.expr <- test.expr[comgene,]
  train.subt <- paste0("CS",moic.res$clust.res$clust)

  # check data
  if(max(train.expr) < 25 | (max(train.expr) >= 25 & min(train.expr) < 0)) {
    message("--training expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    train.expr <- train.expr
  }
  if(max(train.expr) >= 25 & min(train.expr) >= 0){
    message("--log2 transformation done for training expression data.")
    train.expr <- log2(train.expr + 1)
  }
  train.expr <- as.data.frame(t(scale(t(train.expr))))

  if(max(test.expr) < 25 | (max(test.expr) >= 25 & min(test.expr) < 0)) {
    message("--testing expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    test.expr <- test.expr
  }
  if(max(test.expr) >= 25 & min(test.expr) >= 0){
    message("--log2 transformation done for testing expression data.")
    test.expr <- log2(test.expr + 1)
  }
  test.expr <- as.data.frame(t(scale(t(test.expr))))

  # train a partition around medoids (PAM) classifier in training dataset
  mylist <- list(x = as.matrix(train.expr), # x should be the expression matrix
                 y = as.vector(train.subt)) # y should be a vector of classification
  pamr.classifier <- quiet(pamr.train(mylist))

  # classify samples using centroids and calculates the in-group proportions
  # pamr.pred.train <- IGP.clusterRepro(Centroids = pamr.classifier$centroids,
  #                                     Data      = as.matrix(train.expr))
  pamr.pred.test  <- IGP.clusterRepro(Centroids = pamr.classifier$centroids,
                                      Data      = as.matrix(test.expr))

  IGP <- pamr.pred.test$IGP; names(IGP) <- paste0("CS", 1:length(IGP))

  ex.moic.res <- data.frame(samID = names(pamr.pred.test$Class),
                            clust = as.character(pamr.pred.test$Class),
                            row.names = names(pamr.pred.test$Class),
                            stringsAsFactors = FALSE)

  return(list(IGP = IGP, clust.res = ex.moic.res, mo.method = "PAM"))
}
