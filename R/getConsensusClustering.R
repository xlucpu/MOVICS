#' @name getConsensusClustering
#' @title Get subtypes from ConsensusClustering
#' @description This function wraps the Consensus Clustering algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters.
#' @param norMethod A string vector indicate the normalization method for consensus clustering.
#' @param reps An integer value to indicate the number of subsamples.
#' @param pItem A numerical value to indicate the proportion of items to sample.
#' @param pFeature A numerical value to indicate the proportion of features to sample.
#' @param clusterAlg A string value to indicate the cluster algorithm.
#' @param innerLinkage A string value to indicate the heirachical linakge method for subsampling.
#' @param finalLinkage A string value to indicate the heirarchical method for consensus matrix.
#' @param distance A string value to indicate the distance function.
#' @param plot A string value to indicate the output format for heatmap.
#' @param writeTable A logical value to indicate if writing output and log to csv.
#' @param title A string value for output directory.
#' @param seed A numerical value to set random seed for reproducible results.
#' @param verbose A logical value to indicate if printing messages to the screen to indicate progress.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @return A list with the following components:
#'
#'         \code{fit}        an object returned by \link[ConsensusClusterPlus]{ConsensusClusterPlus}.
#'
#'         \code{clust.res}  a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{clust.dend} a dendrogram of sample clustering.
#'
#'         \code{mo.method}  a string value indicating the method used for multi-omics integrative clustering.
#' @export
#' @examples # There is no example and please refer to vignette.
#' @import ConsensusClusterPlus
#' @importFrom dplyr %>%
#' @references Monti S, Tamayo P, Mesirov J, et al (2003). Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data. Mach Learn, 52:91-118.
getConsensusClustering <- function(data         = NULL,
                                   N.clust      = NULL,
                                   type         = rep("gaussian", length(data)),
                                   norMethod    = "none",
                                   reps         = 500,
                                   pItem        = 0.8,
                                   pFeature     = 0.8,
                                   clusterAlg   = "hc",
                                   innerLinkage = "ward.D",
                                   finalLinkage = "ward.D",
                                   distance     = "pearson",
                                   plot         = NULL,
                                   writeTable   = F,
                                   title        = file.path(getwd(),"consensuscluster"),
                                   seed         = 123456,
                                   verbose      = F){

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  useless.argument <- type
  if(is.null(norMethod)) {
    d <- do.call(rbind, data)
  } else {
    if(!is.element(norMethod, c("median-centered","mean-centered","z-score","none"))) {
      stop("the normalized method should be one of median-centered, mean-centered, z-score or none!")
    }
    if(norMethod == "median-centered") {
      d <- do.call(rbind, data)
      d <- sweep(d,1, apply(d,1,median,na.rm=T))
    }
    if(norMethod == "mean-centered") {
      d <- do.call(rbind, data)
      d <- sweep(d,1, apply(d,1,mean,na.rm=T))
    }
    if(norMethod == "z-score") {
      d <- do.call(rbind, data)
      d <- t(scale(t(d)))
    }
    if(norMethod == "none") {
      d <- do.call(rbind, data)
    }
  }

    fit <-  ConsensusClusterPlus(d            = as.matrix(d),
                                 maxK         = ifelse(N.clust == 2, 3, N.clust), # cannot set as 2
                                 reps         = reps,
                                 pItem        = pItem,
                                 pFeature     = pFeature,
                                 clusterAlg   = clusterAlg,
                                 innerLinkage = innerLinkage,
                                 finalLinkage = finalLinkage,
                                 distance     = distance,
                                 seed         = seed,
                                 verbose      = verbose,
                                 plot         = plot,
                                 writeTable   = writeTable,
                                 title        = title)
  res <- fit[[N.clust]]

  clustres <- data.frame(samID = colnames(data[[1]]),
                         clust = as.numeric(res$consensusClass),
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = F)
  #clustres <- clustres[order(clustres$clust),]

  return(list(fit = fit, clust.res = clustres, clust.dend = res$consensusTree, mo.method = "ConsensusClustering"))
}
