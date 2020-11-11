#' @name getPINSPlus
#' @title Get subtypes from PINSPlus
#' @description This function wraps the PINSPlus (Perturbation Clustering for data INtegration and disease Subtyping) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters
#' @param clusteringMethod The name of built-in clustering algorithm that PerturbationClustering will use. Currently supported algorithm are kmeans, pam and hclust. Default value is "kmeans".
#' @param iterMin The minimum number of iterations. Default value is 50
#' @param iterMax The maximum number of iterations. Default value is 500.
#' @param norMethod A string vector indicate the normalization method for consensus clustering.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @return A list with the following components:
#'
#'         \code{fit}       an object returned by \link[PINSPlus]{PerturbationClustering}.
#'
#'         \code{clust.res} a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{mo.method} a string value indicating the method used for multi-omics integrative clustering.
#' @export
#' @examples # There is no example and please refer to vignette.
#' @importFrom PINSPlus PerturbationClustering
#' @importFrom dplyr %>%
#' @references Nguyen H, Shrestha S, Draghici S, Nguyen T (2019). PINSPlus: a tool for tumor subtype discovery in integrated genomic data. Bioinformatics, 35(16):2843-2846.
getPINSPlus <- function(data             = NULL,
                        N.clust          = NULL,
                        type             = rep("gaussian", length(data)),
                        norMethod        = "none",
                        clusteringMethod = "kmeans",
                        iterMin          = 50,
                        iterMax          = 500){

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
      d <- sweep(d,1, apply(d, 1, median, na.rm = TRUE))
    }
    if(norMethod == "mean-centered") {
      d <- do.call(rbind, data)
      d <- sweep(d,1, apply(d, 1, mean, na.rm = TRUE))
    }
    if(norMethod == "z-score") {
      d <- do.call(rbind, data)
      d <- t(scale(t(d)))
    }
    if(norMethod == "none") {
      d <- do.call(rbind, data)
    }
  }

  if(!is.element(clusteringMethod, c("kmeans", "hclust", "pam"))) {
    stop("clusteringMethod should be one of kmeans, hclust, or pam!")
  }

  data <- t(d)

  # for multi-omics but cannot determine cluster number
  # fit <- SubtypingOmicsData(data,
  #                           kMin = N.clust,
  #                           kMax = N.clust,
  #                           clusteringMethod = clusteringMethod,
  #                           iterMin = iterMin,
  #                           iterMax = iterMax,
  #                           verbose = T)

  # for one "feature" but can determine cluster number
  fit <- PerturbationClustering(data             = data,
                                kMin             = N.clust,
                                kMax             = N.clust,
                                clusteringMethod = clusteringMethod,
                                iterMin          = iterMin,
                                iterMax          = iterMax,
                                verbose          = TRUE)

  clustres <- data.frame(samID = rownames(data),
                         clust = fit$cluster,
                         row.names = rownames(data),
                         stringsAsFactors = FALSE)
  #clustres <- clustres[order(clustres$clust),]

  return(list(fit = fit, clust.res = clustres, mo.method = "PINSPlus"))
}
