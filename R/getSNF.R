#' @name getSNF
#' @title Get subtypes from SNF
#' @description This function wraps the SNF (Similarity Network Fusion) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters.
#' @param K An integer value to indicate the number of neighbors in K-nearest neighbors part of the algorithm.
#' @param t An integer value to indicate the number of interations for the diffusion process.
#' @param sigma A numerical value to indicate the variance for local model.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @return A list with the following components:
#'
#'         \code{fit}       an object returned by \link[SNFtool]{SNF}.
#'
#'         \code{clust.res} a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{mo.method} a string value indicating the method used for multi-omics integrative clustering.
#' @export
#' @examples # There is no example and please refer to vignette.
#' @import SNFtool
#' @importFrom dplyr %>%
#' @references Wang B, Mezlini AM, Demir F, et al (2014). Similarity network fusion for aggregating data types on a genomic scale. Nat Methods, 11(3):333-337.
getSNF <- function(data    = NULL,
                   N.clust = NULL,
                   type    = rep("gaussian", length(data)),
                   K       = 30,
                   t       = 20,
                   sigma   = 0.5){

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  useless.argument <- type
  data <- lapply(data, t)

  dat <- lapply(data, function (dd){
    dd <- dd %>% as.matrix
    W <- dd %>% SNFtool::dist2(dd) %>% SNFtool::affinityMatrix(K = K, sigma = sigma)
  })
  W <-  SNFtool::SNF(Wall = dat,
                     K    = K,
                     t    = t)
  clust.SNF = W %>% SNFtool::spectralClustering(N.clust)

  clustres <- data.frame(samID = rownames(data[[1]]),
                         clust = clust.SNF,
                         row.names = rownames(data[[1]]),
                         stringsAsFactors = FALSE)
  #clustres <- clustres[order(clustres$clust),]

  return(list(fit = W, clust.res = clustres, mo.method = "SNF"))
}
