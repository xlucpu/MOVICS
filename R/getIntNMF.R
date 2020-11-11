#' @name getIntNMF
#' @title Get subtypes from IntNMF
#' @description This function wraps the IntNMF (Integrative Clustering via Non-negative Matrix Factorization) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices.
#' @param N.clust Number of clusters.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @return A list with the following components:
#'
#'         \code{fit}       an object returned by \link[IntNMF]{nmf.mnnals}.
#'
#'         \code{clust.res} a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{mo.method} a string value indicating the method used for multi-omics integrative clustering.
#' @import IntNMF
#' @importFrom dplyr %>%
#' @export
#' @examples # There is no example and please refer to vignette.
#' @references Chalise P, Fridley BL (2017). Integrative clustering of multi-level omic data based on non-negative matrix factorization algorithm. PLoS One, 12(5):e0176278.
getIntNMF <- function(data      = NULL,
                      N.clust   = NULL,
                      type      = rep("gaussian", length(data))){

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  # remove features that made of categories not equal to 2 otherwise Error in svd(X) : a dimension is zero
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
  }

  # In order to make the input data fit non-negativity constraint of intNMF,
  # the values of the data were shifted to positive direction by adding absolute value of the smallest negative number.
  # Further, each data was rescaled by dividing by maximum value of the data to make the magnitudes comparable (between 0 and 1) across the several datasets.
  dat <- lapply(data, function (dd){
    if (!all(dd >= 0)) dd <- pmax(dd + abs(min(dd)), 0) + .Machine$double.eps # .Machine$double.eps as The smallest positive floating-point number x
    dd <- dd/max(dd)
    return(dd %>% as.matrix)
  })

  # The function nmf.mnnals requires the samples to be on rows and variables on columns.
  #dat <- lapply(dat, t)
  dat <- lapply(dat, function(x) t(x) + .Machine$double.eps)

  result.intNMF <- dat %>% IntNMF::nmf.mnnals(k = N.clust)
  clust.intNMF <- result.intNMF$clusters

  clustres <- data.frame(samID = colnames(data[[1]]),
                         clust = as.numeric(clust.intNMF),
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = FALSE)
  #clustres <- clustres[order(clustres$clust),]

  return(list(fit = result.intNMF, clust.res = clustres, mo.method = "IntNMF"))
}
