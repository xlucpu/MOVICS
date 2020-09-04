#' @name getLRAcluster
#' @title Get subtypes from LRAcluster
#' @description This function wraps the LRAcluster (Integrated cancer omics data anlsysi by low rank approximation) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters.
#' @param clusterAlg A string value to indicate the cluster algorithm for similarity matrix; 'ward.D' by default.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion; 'gaussian' by default.
#' @return A list with the following components:
#'
#'         \code{fit}        an object returned by \link[LRAcluster]{LRAcluster}.
#'
#'         \code{clust.res}  a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{clust.dend} a dendrogram of sample clustering.
#'
#'         \code{mo.method}  a string value indicating the method used for multi-omics integrative clustering.
#' @examples # There is no example and please refer to vignette.
#' @import LRAcluster
#' @importFrom dplyr %>%
#' @export
#' @references Wu D, Wang D, Zhang MQ, Gu J (2015). Fast dimension reduction and integrative clustering of multi-omics data using low-rank approximation: application to cancer molecular classification. BMC Genomics, 16(1):1022.
getLRAcluster <- function(data       = NULL,
                          N.clust    = NULL,
                          type       = rep("gaussian", length(data)),
                          clusterAlg = "ward.D"){

  data <- lapply(data, as.matrix)

  if(is.element("binomial",type)) {
    bindex <- which(type == "binomial")
    a <- which(rowSums(data[[bindex]]) == 0)
    b <- which(rowSums(data[[bindex]]) == ncol(data[[bindex]]))
    if(length(a) > 0) {
      data[[bindex]] <- data[[bindex]][which(rowSums(data[[bindex]]) != 0),] # remove all zero
    }

    if(length(b) > 0) {
      data[[bindex]] <- data[[bindex]][which(rowSums(data[[bindex]]) != ncol(data[[bindex]])),] # remove all one
    }

    if(length(a) + length(b) > 0) {
      message(paste0("remove a total of ",length(a) + length(b), " features because their categories are not equal to 2!"))
    }
    type[bindex] <- "binary"
  }

  fit <- LRAcluster(data, dimension = N.clust, types = as.list(type))
  dist <- fit$coordinate %>% t %>% dist
  clust.dend <- hclust(dist, method = clusterAlg)

  clustres <- data.frame(samID = colnames(data[[1]]),
                         clust = cutree(clust.dend,k = N.clust),
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = FALSE)
  #clustres <- clustres[order(clustres$clust),]

  return(list(fit = fit, clust.res = clustres, clust.dend = clust.dend, mo.method = "LRAcluster"))
}
