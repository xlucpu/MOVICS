#' @name getiClusterBayes
#' @title Get subtypes from iClusterBayes
#' @description This function wraps the iClusterBayes (Integrative clustering by Bayesian latent variable model) algorithm and provides standard output for `getMoHeatmap()` and `getConsensusMOIC()`.
#' @param data List of matrices with maximum of 6 subdatasets.
#' @param N.clust Number of clusters.
#' @param type Data type corresponding to the list of matrics, which can be gaussian, binomial or possion.
#' @param n.burnin An integer value to indicate the number of MCMC burnin.
#' @param n.draw An integer value to indicate the number of MCMC draw.
#' @param prior.gamma A numerical vector to indicate the prior probability for the indicator variable gamma of each subdataset.
#' @param sdev A numerical value to indicate the standard deviation of random walk proposal for the latent variable.
#' @param thin A numerical value to thin the MCMC chain in order to reduce autocorrelation.
#' @export
#' @return A list with the following components:
#'
#'         \code{fit}       an object returned by \link[iClusterPlus]{iClusterBayes}.
#'
#'         \code{clust.res} a data.frame storing sample ID and corresponding clusters.
#'
#'         \code{feat.res}  the results of features selection process.
#'
#'         \code{mo.method} a string value indicating the method used for multi-omics integrative clustering.
#' @import iClusterPlus
#' @importFrom dplyr %>%
#' @examples # There is no example and please refer to vignette.
#' @references Mo Q, Shen R, Guo C, Vannucci M, Chan KS, Hilsenbeck SG (2018). A fully Bayesian latent variable model for integrative clustering analysis of multi-type omics data. Biostatistics, 19(1):71-86.
getiClusterBayes <- function(data        = NULL,
                             N.clust     = NULL,
                             type        = rep("gaussian", length(data)),
                             n.burnin    = 18000,
                             n.draw      = 12000,
                             prior.gamma = rep(0.5,length(data)),
                             sdev        = 0.05,
                             thin        = 3) {

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  # remove features that made of categories not equal to 2
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

  data.backup <- data
  data <- lapply(data, t)

  if(n_dat == 1){
    res <- iClusterBayes(dt1         = data[[1]],
                         type        = type,
                         K           = N.clust - 1,
                         n.burnin    = n.burnin,
                         n.draw      = n.draw,
                         prior.gamma = prior.gamma,
                         sdev        = sdev,
                         thin        = thin)
  }
  if(n_dat == 2){
    res <- iClusterBayes(dt1         = data[[1]],
                         dt2         = data[[2]],
                         type        = type,
                         K           = N.clust - 1,
                         n.burnin    = n.burnin,
                         n.draw      = n.draw,
                         prior.gamma = prior.gamma,
                         sdev        = sdev,
                         thin        = thin)
  }
  if(n_dat == 3){
    res <- iClusterBayes(dt1         = data[[1]],
                         dt2         = data[[2]],
                         dt3         = data[[3]],
                         type        = type,
                         K           = N.clust - 1,
                         n.burnin    = n.burnin,
                         n.draw      = n.draw,
                         prior.gamma = prior.gamma,
                         sdev        = sdev,
                         thin        = thin)
  }
  if(n_dat == 4){
    res <- iClusterBayes(dt1         = data[[1]],
                         dt2         = data[[2]],
                         dt3         = data[[3]],
                         dt4         = data[[4]],
                         type        = type,
                         K           = N.clust - 1,
                         n.burnin    = n.burnin,
                         n.draw      = n.draw,
                         prior.gamma = prior.gamma,
                         sdev        = sdev,
                         thin        = thin)
  }
  if(n_dat == 5){
    res <- iClusterBayes(dt1         = data[[1]],
                         dt2         = data[[2]],
                         dt3         = data[[3]],
                         dt4         = data[[4]],
                         dt5         = data[[5]],
                         type        = type,
                         K           = N.clust - 1,
                         n.burnin    = n.burnin,
                         n.draw      = n.draw,
                         prior.gamma = prior.gamma,
                         sdev        = sdev,
                         thin        = thin)
  }
  if(n_dat == 6){
    res <- iClusterBayes(dt1         = data[[1]],
                         dt2         = data[[2]],
                         dt3         = data[[3]],
                         dt4         = data[[4]],
                         dt5         = data[[5]],
                         dt6         = data[[6]],
                         type        = type,
                         K           = N.clust - 1,
                         n.burnin    = n.burnin,
                         n.draw      = n.draw,
                         prior.gamma = prior.gamma,
                         sdev        = sdev,
                         thin        = thin)
  }
  message("clustering done...")
  clustres <- data.frame(samID = rownames(data[[1]]),
                         clust = res$clusters,
                         row.names = rownames(data[[1]]),
                         stringsAsFactors = FALSE)
  #clustres <- clustres[order(clustres$clust),]

  featres <- data.frame(feature = as.character(unlist(lapply(data.backup, rownames))),
                        dataset = rep(names(data),as.numeric(sapply(data.backup, function(x) length(rownames(x))))),
                        post.prob = unlist(res$beta.pp),
                        stringsAsFactors = FALSE)

  feat.res <- NULL
  for (d in unique(featres$dataset)) {
    tmp <- featres[which(featres$dataset == d),]
    tmp <- tmp[order(tmp$post.prob, decreasing = TRUE),]
    tmp$rank <- 1:nrow(tmp)
    feat.res <- rbind.data.frame(feat.res,tmp)
  }
  message("feature selection done...")

  return(list(fit = res, clust.res = clustres, feat.res = feat.res, mo.method = "iClusterBayes"))
}
