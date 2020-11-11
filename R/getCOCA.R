#' @name getCOCA
#' @title Get subtypes from COCA
#' @description This function wraps the COCA (Cluster-of-Clusters Analysis) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters.
#' @param methods A string vector storing the names of clustering methods to be used to cluster the observations in each subdataset.
#' @param distances A string vector storing the name of distances to be used in the clustering step for each subdataset.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @return A list with the following components:
#'
#'         \code{fit}        an object returned by \link[coca]{coca}.
#'
#'         \code{clust.res}  a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{clust.dend} a dendrogram of sample clustering.
#'
#'         \code{mo.method}  a string value indicating the method used for multi-omics integrative clustering.
#' @import coca
#' @importFrom vegan vegdist
#' @export
#' @examples # There is no example and please refer to vignette.
#' @references Hoadley KA, Yau C, Wolf DM, et al (2014). Multiplatform analysis of 12 cancer types reveals molecular classification within and across tissues of origin. Cell, 158(4):929-944.
getCOCA <- function(data      = NULL,
                    N.clust   = NULL,
                    type      = rep("gaussian", length(data)),
                    methods   = "hclust",
                    distances = "euclidean") {

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

  ### Build matrix of clusters
  outputBuildMOC <- coca::buildMOC(data,
                                   M         = length(data),
                                   K         = N.clust,
                                   methods   = methods,
                                   distances = distances)

  ### Extract matrix of clusters and dataset indicator vector
  moc <- outputBuildMOC$moc
  datasetIndicator <- outputBuildMOC$datasetIndicator

  hcs <- hclust(vegdist(as.matrix(moc), method = "jaccard"), "ward.D")
  coca <- cutree(hcs,N.clust)
  #coca <- coca::coca(moc, K = N.clust)

  clustres <- data.frame(samID = rownames(data[[1]]),
                         clust = as.numeric(coca),
                         row.names = rownames(data[[1]]),
                         stringsAsFactors = FALSE)
  #clustres <- clustres[order(clustres$clust),]

  return(list(fit = outputBuildMOC, clust.res = clustres, clust.dend = hcs, mo.method = "COCA"))
}
