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
#' @import SNFtool
#' @export
#' @references Rappoport N, Shamir R (2019). NEMO: cancer subtyping by integration of partial multi-omic data. Bioinformatics, 35(18):3348-3356.
#' @examples # There is no example and please refer to vignette.
getNEMO <- function(data          = NULL,
                    N.clust       = NULL,
                    type          = rep("gaussian", length(data)),
                    num.neighbors = NA) {

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  useless.argument <- type
  fit <- nemo.clustering(omics.list    = data,
                         num.clusters  = N.clust,
                         num.neighbors = num.neighbors)

  clustres <- data.frame(samID = colnames(data[[1]]),
                         clust = as.numeric(fit),
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = FALSE)
  #clustres <- clustres[order(clustres$clust),]

  return(list(fit = fit, clust.res = clustres, mo.method = "NEMO"))
}

#' @title affinityMatrix
#' @name affinityMatrix
#' @param Diff Distance Matrix.
#' @param K Number of nearest neighbors.
#' @param sigma Variance for local model.
#' @author Bo Wang, Aziz Mezlini, Feyyaz Demir, Marc Fiume, Zhuowen Tu, Michael Brudno, Benjamin Haibe-Kains, Anna Goldenberg
#' @references Wang B, Mezlini AM, Demir F, et al (2014). Similarity network fusion for aggregating data types on a genomic scale. Nat Methods, 11(3):333-337.
#' @keywords internal
#' @return affinityMatrix
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
#' @references Wang B, Mezlini AM, Demir F, et al (2014). Similarity network fusion for aggregating data types on a genomic scale. Nat Methods, 11(3):333-337.
#' @keywords internal
#' @return dist2
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


#' @title Spectral clustering
#' @name spectralClustering
#' @keywords internal
#' @author Bo Wang, Aziz Mezlini, Feyyaz Demir, Marc Fiume, Zhuowen Tu, Michael Brudno, Benjamin Haibe-Kains, Anna Goldenberg
#' @references Wang B, Mezlini AM, Demir F, et al (2014). Similarity network fusion for aggregating data types on a genomic scale. Nat Methods, 11(3):333-337.
#' @return spectralClustering
spectralClustering = SNFtool::spectralClustering

#' @title NEMO num clusters
#' @name nemo.num.clusters
#' @description Estimates the number of clusters in an affinity graph.
#' @param W the affinity graph.
#' @param NUMC possible values for the number of clusters. Defaults to 2:15.
#' @return the estimated number of clusters in the graph.
#' @author Nimrod Rappoport
#' @references Rappoport N, Shamir R (2019). NEMO: cancer subtyping by integration of partial multi-omic data. Bioinformatics, 35(18):3348-3356.
#' @keywords internal
nemo.num.clusters <- function(W, NUMC=2:15) {
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC = NUMC[NUMC > 1]
  }
  W = (W + t(W))/2
  diag(W) = 0
  if (length(NUMC) > 0) {
    degs = rowSums(W)
    degs[degs == 0] = .Machine$double.eps
    D = diag(degs)
    L = D - W
    Di = diag(1/sqrt(degs))
    L = Di %*% L %*% Di
    print(dim(L))
    eigs = eigen(L)
    eigs_order = sort(eigs$values, index.return = TRUE)$ix
    eigs$values = eigs$values[eigs_order]
    eigs$vectors = eigs$vectors[, eigs_order]
    eigengap = abs(diff(eigs$values))
    eigengap = (1:length(eigengap)) * eigengap

    t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = TRUE)$ix
    return(NUMC[t1[1]])
  }
}

#' @title NEMO affinity graph
#' @name nemo.affinity.graph
#' @description Constructs a single affinity graph measuring similarity across different omics.
#' @param raw.data A list of the data to be clustered, where each an entry is a matrix of features x samples.
#' @param k The number of neighbors to use for each omic. It can either be a number, a list of numbers
#' or NA. If it is a number, this is the number of neighbors used for all omics. If this is a list,
#' the number of neighbors are taken for each omic from that list. If it is NA, each omic chooses the
#' number of neighbors to be the number of samples divided by NUM.NEIGHBORS.RATIO.
#' @return A single matrix measuring similarity between the samples across all omics.
#' @author Nimrod Rappoport
#' @references Rappoport N, Shamir R (2019). NEMO: cancer subtyping by integration of partial multi-omic data. Bioinformatics, 35(18):3348-3356.
#' @keywords internal
nemo.affinity.graph <- function(raw.data, k = NA, NUM.NEIGHBORS.RATIO = 6) {
  if (is.na(k)) {
    k = as.numeric(lapply(1:length(raw.data), function(i) round(ncol(raw.data[[i]]) / NUM.NEIGHBORS.RATIO)))
  } else if (length(k) == 1) {
    k = rep(k, length(raw.data))
  }
  sim.data = lapply(1:length(raw.data), function(i) {affinityMatrix(dist2(as.matrix(t(raw.data[[i]])),
                                                                          as.matrix(t(raw.data[[i]]))), k[i], 0.5)})
  affinity.per.omic = lapply(1:length(raw.data), function(i) {
    sim.datum = sim.data[[i]]
    non.sym.knn = apply(sim.datum, 1, function(sim.row) {
      returned.row = sim.row
      threshold = sort(sim.row, decreasing = TRUE)[k[i]]
      returned.row[sim.row < threshold] = 0
      row.sum = sum(returned.row)
      returned.row[sim.row >= threshold] = returned.row[sim.row >= threshold] / row.sum
      return(returned.row)
    })
    sym.knn = non.sym.knn + t(non.sym.knn)
    return(sym.knn)
  })
  patient.names = Reduce(union, lapply(raw.data, colnames))
  num.patients = length(patient.names)
  returned.affinity.matrix = matrix(0, ncol = num.patients, nrow=num.patients)
  rownames(returned.affinity.matrix) = patient.names
  colnames(returned.affinity.matrix) = patient.names

  shared.omic.count = matrix(0, ncol = num.patients, nrow=num.patients)
  rownames(shared.omic.count) = patient.names
  colnames(shared.omic.count) = patient.names

  for (j in 1:length(raw.data)) {
    curr.omic.patients = colnames(raw.data[[j]])
    returned.affinity.matrix[curr.omic.patients, curr.omic.patients] = returned.affinity.matrix[curr.omic.patients, curr.omic.patients] + affinity.per.omic[[j]][curr.omic.patients, curr.omic.patients]
    shared.omic.count[curr.omic.patients, curr.omic.patients] = shared.omic.count[curr.omic.patients, curr.omic.patients] + 1
  }

  final.ret = returned.affinity.matrix / shared.omic.count
  lower.tri.ret = final.ret[lower.tri(final.ret)]
  final.ret[shared.omic.count == 0] = mean(lower.tri.ret[!is.na(lower.tri.ret)])

  return(final.ret)
}

#' @title NEMO clustering
#' @name nemo.clustering
#' @description Performs multi-omic clustering on a datset using the NEMO algorithm.
#' Uses nemo.num.clusters to estimate the number of clusters.
#' @param omics.list A list of the data to be clustered, where each an entry is a matrix of features x samples.
#' @param k The number of neighbors to use for each omic. It can either be a number, a list of numbers
#' or NA. If it is a number, this is the number of neighbors used for all omics. If this is a list,
#' the number of neighbors are taken for each omic from that list. If it is NA, each omic chooses the
#' number of neighbors to be the number of samples divided by NUM.NEIGHBORS.RATIO.
#' @return A single matrix measuring similarity between the samples across all omics.
#' @author Nimrod Rappoport
#' @references Rappoport N, Shamir R (2019). NEMO: cancer subtyping by integration of partial multi-omic data. Bioinformatics, 35(18):3348-3356.
#' @keywords internal
nemo.clustering <- function(omics.list, num.clusters = NULL, num.neighbors = NA) {
  if (is.null(num.clusters)) {
    num.clusters = NA
  }

  graph = nemo.affinity.graph(omics.list, k = num.neighbors)
  if (is.na(num.clusters)) {
    num.clusters = nemo.num.clusters(graph)
  }
  clustering = spectralClustering(graph, num.clusters)
  names(clustering) = colnames(graph)
  return(clustering)
}
