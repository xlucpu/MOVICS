#' @title affinityMatrix
#' @name affinityMatrix
#' @param Diff Distance Matrix.
#' @param K Number of nearest neighbors.
#' @param sigma Variance for local model.
#' @author Bo Wang, Aziz Mezlini, Feyyaz Demir, Marc Fiume, Zhuowen Tu, Michael Brudno, Benjamin Haibe-Kains, Anna Goldenberg
#' @keywords internal
affinityMatrix <- function(Diff,K=20,sigma=0.5) {
  ###This function constructs similarity networks.
  N = nrow(Diff)

  Diff = (Diff + t(Diff)) / 2
  diag(Diff) = 0;
  sortedColumns = as.matrix(t(apply(Diff,2,sort)))
  finiteMean <- function(x) { mean(x[is.finite(x)]) }
  means = apply(sortedColumns[,1:K+1],1,finiteMean)+.Machine$double.eps;

  avg <- function(x,y) ((x+y)/2)
  Sig = outer(means,means,avg)/3*2 + Diff/3 + .Machine$double.eps;
  Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
  densities = dnorm(Diff,0,sigma*Sig,log = FALSE)

  W = (densities + t(densities)) / 2
  return(W)
}

#' @title dist2
#' @name dist2
#' @param X A data matrix where each row is a different data point.
#' @param C A data matrix where each row is a different data point. If this matrix is the same as X, pairwise distances for all data points are computed.
#' @author Bo Wang, Aziz Mezlini, Feyyaz Demir, Marc Fiume, Zhuowen Tu, Michael Brudno, Benjamin Haibe-Kains, Anna Goldenberg
#' @keywords internal
dist2 <- function(X,C) {
  ndata = nrow(X)
  ncentres = nrow(C)

  sumsqX = rowSums(X^2)
  sumsqC = rowSums(C^2)

  XC = 2 * (X %*% t(C))

  res = matrix(rep(sumsqX,times=ncentres),ndata,ncentres) + t(matrix(rep(sumsqC,times=ndata),ncentres,ndata)) - XC
  res[res < 0] = 0
  return(res)
}

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



