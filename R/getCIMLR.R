#' @name getCIMLR
#' @title Get subtypes from CIMLR
#' @description This function wraps the CIMLR (Cancer Integration via Multikernel Learning) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters.
#' @param cores.ratio Ratio of the number of cores to be used when computing the multi-kernel.
#' @param verbose A logic value to indicate if supressing progression.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @return A list with the following components:
#'
#'         \code{fit}       an object returned by \link[CIMLR]{CIMLR}.
#'
#'         \code{clust.res} a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{feat.res}  the results of features selection process.
#'
#'         \code{mo.method} a string value indicating the method used for multi-omics integrative clustering.
#' @import CIMLR
#' @importFrom SimDesign quiet
#' @export
#' @examples # There is no example and please refer to vignette.
#' @references Ramazzotti D, Lal A, Wang B, Batzoglou S, Sidow A (2018). Multi-omic tumor data reveal diversity of molecular mechanisms that correlate with survival. Nat Commun, 9(1):4453.
getCIMLR <- function(data        = NULL,
                     N.clust     = NULL,
                     type        = rep("gaussian", length(data)),
                     cores.ratio = 0,
                     verbose     = T){

  useless.argument <- type
  if(verbose) {
    fit <- quiet(CIMLR(data,
                       c= N.clust,
                       cores.ratio = cores.ratio))
  } else {
    fit <- CIMLR(data,
                 c= N.clust,
                 cores.ratio = cores.ratio)
  }

  message("clustering done...")
  input_dat <- do.call(rbind,lapply(seq(along = data), function(x){
    ddd <- data[[x]]
    rownames(ddd) <- paste(rownames(ddd), names(data)[x], sep = "+")
    ddd
  }))

  if(verbose) {
    ranks <- quiet(CIMLR_Feature_Ranking(A = fit$S, X = input_dat))
  } else {
    ranks <- CIMLR_Feature_Ranking(A = fit$S, X = input_dat)
  }
  ranks$names <- rownames(input_dat)[ranks$aggR]
  fit$selectfeatures <- ranks
  message("feature selection done...")

  clustres <- data.frame(samID = colnames(data[[1]]),
                         clust = fit$y$cluster,
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = F)
  #clustres <- clustres[order(clustres$clust),]

  f <- sapply(strsplit(ranks$name, "+",fixed = T), "[",1)
  d <- sapply(strsplit(ranks$name, "+",fixed = T), "[",2)

  featres <- data.frame(feature = f,
                        dataset = d,
                        pvalue = ranks$pval,
                        stringsAsFactors = F)
  feat.res <- NULL
  for (d in unique(featres$dataset)) {
    tmp <- featres[which(featres$dataset == d),]
    tmp <- tmp[order(tmp$pvalue, decreasing = F),]
    tmp$rank <- 1:nrow(tmp)
    feat.res <- rbind.data.frame(feat.res,tmp)
  }

  return(list(fit = fit, clust.res = clustres, feat.res = feat.res, mo.method = "CIMLR"))
}
