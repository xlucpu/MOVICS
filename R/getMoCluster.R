#' @name getMoCluster
#' @title Get subtypes from MoCluster
#' @description This function wraps the MoCluster (Multiple omics data integrative clustering) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters.
#' @param ncomp An integer value to indicate the number of components to calculate. To calculate more components requires longer computational time.
#' @param method A string value can be one of CPCA, GCCA and MCIA; CPCA by default.
#' @param option A string value could be one of c('lambda1', 'inertia', 'uniform') to indicate how the different matrices should be normalized.
#' @param k A numeric value to indicate the absolute number (if k >= 1) or the proportion (if 0 < k < 1) of non-zero coefficients for the variable loading vectors. It could be a single value or a vector has the same length as x so the sparsity of individual matrix could be different.
#' @param center A logical value to indicate if the variables should be centered. TRUE by default.
#' @param scale A logical value to indicate if the variables should be scaled. TRUE by default.
#' @param clusterAlg A string value to indicate the cluster algorithm for distance.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @return A list with the following components:
#'
#'         \code{fit}        an object returned by \link[mogsa]{mbpca}.
#'
#'         \code{clust.res}  a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{feat.res}   the results of features selection process.
#'
#'         \code{clust.dend} a dendrogram of sample clustering.
#'
#'         \code{mo.method}  a string value indicating the method used for multi-omics integrative clustering.
#' @export
#' @examples # There is no example and please refer to vignette.
#' @import mogsa
#' @importFrom dplyr %>%
#' @references Meng C, Helm D, Frejno M, Kuster B (2016). moCluster: Identifying Joint Patterns Across Multiple Omics Data Sets. J Proteome Res, 15(3):755-765.
getMoCluster <- function(data       = NULL,
                         N.clust    = NULL,
                         type       = rep("gaussian", length(data)),
                         ncomp      = NULL,
                         method     = "CPCA",
                         option     = "lambda1",
                         k          = 10,
                         center     = TRUE,
                         scale      = TRUE,
                         clusterAlg = "ward.D"){

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  useless.argument <- type
  if(!is.element(method, c("CPCA","GCCA","MCIA"))) {
    stop("method should be one of CPCA [consensus PCA], GCCA [generalized canonical correlation analysis], or MCIA [multiple co-inertia analysis]!")
  }

  if(is.null(ncomp)) {
    ncomp = N.clust
  }

  moas <- data %>% mogsa::mbpca(ncomp      = ncomp,
                                k          = k,
                                method     = switch(method,
                                                    "CPCA" = "globalScore",
                                                    "GCCA" = "blockScore",
                                                    "MCIA" = "blockLoading"),
                                option     = option,
                                center     = center,
                                scale      = scale,
                                moa        = TRUE,
                                svd.solver = "fast",
                                maxiter    = 1000,
                                verbose    = FALSE)

  scrs <- moas %>% moaScore
  dist <- scrs %>% dist
  clust.dend <- hclust(dist, method = clusterAlg)

  clustres <- data.frame(samID = colnames(data[[1]]),
                         clust = cutree(clust.dend,k = N.clust),
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = FALSE)
  #clustres <- clustres[order(clustres$clust),]
  message("clustering done...")


  featres <- moas@loading[which(moas@loading[,1] != 0),]
  f <- sub('_[^_]*$', '', rownames(featres))
  d <- sub('.*_', '', rownames(featres))
  featres <- data.frame(feature = f,
                        dataset = d,
                        load = featres[,1],
                        stringsAsFactors = FALSE)
  feat.res <- NULL
  for (d in unique(featres$dataset)) {
    tmp <- featres[which(featres$dataset == d),]
    feat.res <- rbind.data.frame(feat.res,tmp)
  }
  message("feature selection done...")

  return(list(fit = moas, clust.res = clustres, feat.res = feat.res, clust.dend = clust.dend, mo.method = "MoCluster"))
}
