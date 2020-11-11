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
#' @importFrom dplyr %>%
#' @export
#' @references Wu D, Wang D, Zhang MQ, Gu J (2015). Fast dimension reduction and integrative clustering of multi-omics data using low-rank approximation: application to cancer molecular classification. BMC Genomics, 16(1):1022.
getLRAcluster <- function(data       = NULL,
                          N.clust    = NULL,
                          type       = rep("gaussian", length(data)),
                          clusterAlg = "ward.D"){

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  data <- lapply(data, as.matrix)

  if(is.element("binomial",type)) {
    bindex <- which(type == "binomial")
    for (i in bindex) {
      a <- which(rowSums(data[[i]]) == 0)
      b <- which(rowSums(data[[i]]) == ncol(data[[i]]))
      if(length(a) > 0) {
        data[[i]] <- data[[i]][which(rowSums(data[[i]]) != 0),] # remove all zero
      }

      if(length(b) > 0) {
        data[[i]] <- data[[i]][which(rowSums(data[[i]]) != ncol(data[[i]])),] # remove all one
      }

      if(length(a) + length(b) > 0) {
        message(paste0("--", names(data)[i],": a total of ",length(a) + length(b), " features were removed due to the categories were not equal to 2!"))
      }
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
