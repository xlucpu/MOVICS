#' @name getNEMO
#' @title Get subtypes from NEMO
#' @description This function wraps the NEMO (Neighborhood based multi-omics clustering) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters.
#' @param num.neighbors The number of neighbors to use for each omic.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @return A list with the following components:
#'
#'         \code{fit}       an object returned by \link[NEMO]{nemo.clustering}.
#'
#'         \code{clust.res} a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{mo.method} a string value indicating the method used for multi-omics integrative clustering.
#' @import NEMO
#' @importFrom SNFtool affinityMatrix dist2
#' @export
#' @references Rappoport N, Shamir R (2019). NEMO: cancer subtyping by integration of partial multi-omic data. Bioinformatics, 35(18):3348-3356.
#' @examples # There is no example and please refer to vignette.
getNEMO <- function(data          = NULL,
                    N.clust       = NULL,
                    type          = rep("gaussian", length(data)),
                    num.neighbors = NA) {

  useless.argument <- type
  fit <- nemo.clustering(omics.list    = data,
                         num.clusters  = N.clust,
                         num.neighbors = num.neighbors)

  clustres <- data.frame(samID = colnames(data[[1]]),
                         clust = as.numeric(fit),
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = F)
  #clustres <- clustres[order(clustres$clust),]

  return(list(fit = fit, clust.res = clustres, mo.method = "NEMO"))
}



