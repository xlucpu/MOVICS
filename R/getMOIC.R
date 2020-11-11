#' @name getMOIC
#' @title Get subtypes from multi-omics integrative clustering
#' @description Using `getMOIC()`, users can choose one out of the ten algorithms embedded in `MOVICS`. Users can implement multi-omics clustering in a simplest way of which the only requirement is to specify and at least specify a list of matrices (argument of `data`), a number of cluster (argument of `N.clust`), and clustering method (argument of `methodslist`) in `getMOIC()`. It is possible to pass various arguments that are specific to each method. Of course, users can also directly call different algorithms by using functions start with `get` and end with the name of the algorithm (e.g., `getSNF`; please refer to `?get%algorithm_name%` for more details about the editable arguments)
#' @param data List of matrices (Maximum number of matrices is 6).
#' @param methodslist A string list specifying one or multiple methods to run (See Details).
#' @param N.clust Number of clusters.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @param ... Additionnal parameters for each method (only works when only one method chosen)
#' @examples # There is no example and please refer to vignette.
#' @export
#' @return A list of results returned by each specified algorithms.
#' @import SNFtool
#' @import IntNMF
#' @import mogsa
#' @import coca
#' @import iClusterPlus
#' @import CIMLR
#' @import PINSPlus
#' @import ConsensusClusterPlus
#' @details
#' Method for integrative clustering will be chosed according to the value of argument 'methodslist':
#'
#' If \code{methodslist == "IntNMF"}, Integrative clustering methods using Non-Negative Matrix Factorization
#'
#' If \code{methodslist == "SNF"}, Similarity network fusion.
#'
#' If \code{methodslist == "LRAcluster"}, Integrated cancer omics data analysis by low rank approximation.
#'
#' If \code{methodslist == "PINSPlus"}, Perturbation Clustering for data integration and disease subtyping
#'
#' If \code{methodslist == "ConsensusClustering"}, Consensus clustering
#'
#' If \code{methodslist == "NEMO"}, Neighborhood based multi-omics clustering
#'
#' If \code{methodslist == "COCA"}, Cluster Of Clusters Analysis
#'
#' If \code{methodslist == "CIMLR"}, Cancer Integration via Multikernel Learning (Support Feature Selection)
#'
#' If \code{methodslist == "MoCluster"}, Identifying joint patterns across multiple omics data sets (Support Feature Selection)
#'
#' If \code{methodslist == "iClusterBayes"}, Integrative clustering of multiple genomic data by fitting a Bayesian latent variable model (Support Feature Selection)
#'
#' @references
#' Pierre-Jean M, Deleuze J F, Le Floch E, et al. Clustering and variable selection evaluation of 13 unsupervised methods for multi-omics data integration[J]. Briefings in Bioinformatics, 2019.
#'
#' intNMF:
#' Chalise P, Fridley BL. Integrative clustering of multi-level omic data based on non-negative matrix factorization algorithm. PLoS One. 2017;12(5):e0176278.
#'
#' iClusterBayes:
#' Mo Q, Shen R, Guo C, Vannucci M, Chan KS, Hilsenbeck SG. A fully Bayesian latent variable model for integrative clustering analysis of multi-type omics data. Biostatistics. 2018;19(1):71-86.
#'
#' SNF:
#' Wang B, Mezlini AM, Demir F, et al. Similarity network fusion for aggregating data types on a genomic scale. Nat Methods. 2014;11(3):333-337.
#'
#' Mocluster:
#' Meng C, Helm D, Frejno M, Kuster B. moCluster: Identifying Joint Patterns Across Multiple Omics Data Sets. J Proteome Res. 2016;15(3):755-765.
#'
#' LRAcluster:
#' Wu D, Wang D, Zhang MQ, Gu J. Fast dimension reduction and integrative clustering of multi-omics data using low-rank approximation: application to cancer molecular classification. BMC Genomics. 2015;16:1022.
#'
#' CIMLR:
#' Ramazzotti D, Lal A, Wang B, Batzoglou S, Sidow A. Multi-omic tumor data reveal diversity of molecular mechanisms that correlate with survival. Nat Commun. 2018;9(1):4453.
#'
#' PINSPlus:
#' Nguyen H, Shrestha S, Draghici S, Nguyen T. PINSPlus: a tool for tumor subtype discovery in integrated genomic data. Bioinformatics. 2019;35(16):2843-2846.
#'
#' ConsensusClustering:
#' Monti S, Tamayo P, Mesirov J, et al. Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data. Machine Learning. 2003;52:91-118.
#'
#' NEMO:
#' Rappoport N, Shamir R. NEMO: cancer subtyping by integration of partial multi-omic data. Bioinformatics. 2019;35(18):3348-3356.
#'
#' COCA:
#' Hoadley KA, Yau C, Wolf DM, et al. Multiplatform analysis of 12 cancer types reveals molecular classification within and across tissues of origin. Cell. 2014;158(4):929-944.
getMOIC <- function(data        = NULL,
                    methodslist = list("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "ConsensusClustering", "IntNMF", "iClusterBayes"),
                    N.clust     = NULL,
                    type        = rep("gaussian", length(data)),
                    ...){

  # check argument
  if (!is.list(data)) {
    stop("data is not a list!")
  }

  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 omics data.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  if(is.null(names(data))){
    names(data) <- sprintf("dat%s", 1:length(data))
  }

  num.methods <- length(unlist(methodslist))

  if(is.vector(methodslist)) {methodslist <- as.list(methodslist)}
  if(!all(is.element(unlist(methodslist), c("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster", "LRAcluster", "ConsensusClustering", "IntNMF", "iClusterBayes")))) {
    stop("current version of MOVICS supports 10 algorithms. Allowed values contain c('SNF', 'CIMLR', 'PINSPlus', 'NEMO', 'COCA', 'MoCluster', 'LRAcluster', 'ConsensusClustering', 'IntNMF', 'iClusterBayes').")
  }

  if(num.methods > 1) {
    message("--you choose more than 1 algorithm and all of them shall be run with parameters by default.")
  }

  # Check dimension
  if(max(sapply(data, dim)[2,]) != min(sapply(data, dim)[2,])){
    message(sprintf("number of samples in dat %s is %s\n", 1:length(data), sapply(data, dim)[2,]))
    stop("data do not contain the same number of samples!")
  }
  reslist <- list()
  for (method in unlist(methodslist)) {
    doMOIC <- switch(method,
                     "IntNMF"              = getIntNMF,
                     "iClusterBayes"       = getiClusterBayes,
                     "SNF"                 = getSNF,
                     "MoCluster"           = getMoCluster,
                     "LRAcluster"          = getLRAcluster,
                     "CIMLR"               = getCIMLR,
                     "PINSPlus"            = getPINSPlus,
                     "ConsensusClustering" = getConsensusClustering,
                     "NEMO"                = getNEMO,
                     "COCA"                = getCOCA
    )
    reslist[[method]] <- doMOIC(data, N.clust, type, ...)
    message(paste0(method," done..."))
  }

  if(num.methods == 1) {
    return(reslist[[1]])
  } else {
    return(reslist)
  }
}
