#' @name getClustNum
#' @title Get estimation of optimal clustering number
#' @description This function provides two measurements (i.e., clustering prediction index [CPI] and Gap-statistics) and aims to search the optimal number for multi-omics integrative clustering. In short, the peaks reach by the red (CPI) and blue (Gap-statistics) lines should be referred to determine `N.clust`.
#' @param data List of matrices.
#' @param is.binary A logicial vector to indicate if the subdata is binary matrix of 0 and 1 such as mutation.
#' @param try.N.clust A integer vector to indicate possible choices of number of clusters.
#' @param center A logical value to indicate if the variables should be centered. TRUE by default.
#' @param scale A logical value to indicate if the variables should be scaled. FALSE by default.
#' @param fig.path A string value to indicate the output figure path.
#' @param fig.name A string value to indicate the name of the figure.
#' @export
#' @return A figure that helps to choose the optimal clustering number (argument of `N.clust`) for `get%algorithm_name%()` or `getMOIC()`, and a list contains the following components:
#'
#'         \code{CPI}   possible cluster number identified by clustering prediction index
#'
#'         \code{Gapk}  possible cluster number identified by Gap-statistics
#' @import IntNMF
#' @import mogsa
#' @import SNFtool
#' @importFrom ggplot2 alpha
#' @importFrom dplyr %>%
#' @importFrom grDevices dev.copy2pdf
#' @examples # There is no example and please refer to vignette.
#' @references Chalise P, Fridley BL (2017). Integrative clustering of multi-level omic data based on non-negative matrix factorization algorithm. PLoS One, 12(5):e0176278.
#'
#'Tibshirani, R., Walther, G., Hastie, T. (2001). Estimating the number of data clusters via the Gap statistic. J R Stat Soc Series B Stat Methodol, 63(2):411-423.
getClustNum <- function(data        = NULL,
                        is.binary   = rep(FALSE, length(data)),
                        try.N.clust = 2:8,
                        center      = TRUE,
                        scale       = TRUE,
                        fig.path    = getwd(),
                        fig.name    = "optimal_number_cluster") {

  # check data
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  data.backup <- data # save a backup

  #--------------------------------------------#
  # Cluster Prediction Index (CPI) from IntNMF #
  # remove features that made of categories not equal to 2 otherwise Error in svd(X) : a dimension is zero
  if(!all(!is.binary)) {
    bindex <- which(is.binary == TRUE)
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

  #dat <- lapply(dat, t)
  dat <- lapply(dat, function(x) t(x) + .Machine$double.eps)

  message("calculating Cluster Prediction Index...")
  optk1 <- IntNMF::nmf.opt.k(dat      = dat,
                            n.runs    = 5,
                            n.fold    = 5,
                            k.range   = try.N.clust,
                            st.count  = 10,
                            maxiter   = 100,
                            make.plot = FALSE)
  optk1 <- as.data.frame(optk1)

  #-------------------------------#
  # Gap-statistics from MoCluster #
  message("calculating Gap-statistics...")
  moas <- data.backup %>% mogsa::mbpca(ncomp      = 2,
                                       k          = "all",
                                       method     = "globalScore",
                                       center     = center,
                                       scale      = scale,
                                       moa        = TRUE,
                                       svd.solver = "fast",
                                       maxiter    = 1000,
                                       verbose    = FALSE)
  gap <- mogsa::moGap(moas, K.max = max(try.N.clust), cluster = "hclust", plot = FALSE)
  optk2 <- as.data.frame(gap$Tab)[-1,] # remove k=1

  #---------------------#
  # Eigen-gaps from SNF #
  # message("Calculating Eigen-gaps...")
  # data <- lapply(data.backup, t)
  # dat <- lapply(data, function (dd){
  #   dd <- dd %>% as.matrix
  #   W <- dd %>% SNFtool::dist2(dd) %>% SNFtool::affinityMatrix(K = 30, sigma = 0.5)
  # })
  # W <-  SNFtool::SNF(Wall = dat,
  #                    K = 30,
  #                    t = 20)
  # optk3 <- SNFtool::estimateNumberOfClustersGivenGraph(W, NUMC = try.N.clust)

  # nemo.ag <- nemo.affinity.graph(data.backup, k = 20)
  # optk4 <- nemo.num.clusters(nemo.ag, NUMC = try.N.clust)

  #---------------------------#
  # calculate optimal N.clust #
  # N.clust <- as.numeric(names(which(table(c(as.numeric((which.max(apply(optk1, 1, mean)) + 1)),
  #                                           as.numeric((which.max(optk2$gap) + 1)),
  #                                           optk3$`Eigen-gap best`)) >= 2)))

  N.clust <- as.numeric(which.max(apply(optk1, 1, mean) + optk2$gap)) + 1
  if(length(N.clust) == 0) {
    message("--fail to define the optimal cluster number!")
    N.clust <- "null"
  }

  #---------------#
  # Visualization #
  outFig <- paste0(fig.name,".pdf")
  par(bty="o", mgp = c(1.9,.33,0), mar=c(3.1,3.1,2.1,3.1)+.1, las=1, tcl=-.25)
  plot(NULL, NULL,
       xlim = c(min(try.N.clust),max(try.N.clust)),
       #ylim = c(min(optk1), max(optk1)),
       ylim = c(0,1),
       xlab = "Number of Multi-Omics Clusters",ylab = "")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#EAE9E9",border = FALSE)
  grid(col = "white", lty = 1, lwd = 1.5)
  # for (m in 1:n.runs) points(try.N.clust, optk1[, m], pch = 20, cex = 1.5, col = "#224A8D")
  points(try.N.clust, apply(optk1, 1, mean), pch = 19, col = ggplot2::alpha("#224A8D"), cex = 1.5)
  lines(try.N.clust, apply(optk1, 1, mean), col = "#224A8D", lwd = 2, lty = 4)
  mtext("Cluster Prediction Index", side = 2, line = 2, cex = 1.5, col = "#224A8D", las = 3)

  par(new = TRUE, xpd = FALSE)
  plot(NULL,NULL,
       xlim = c(min(try.N.clust),max(try.N.clust)),
       #ylim = c(min(optk2$gap), max(optk2$gap)),
       ylim = c(0,1),
       xlab = "",ylab = "",xaxt = "n",yaxt = "n")
  points(try.N.clust, optk2$gap, pch = 19, col = ggplot2::alpha("#E51718",0.8), cex = 1.5)
  lines(try.N.clust, optk2$gap, col = "#E51718", lwd = 2, lty = 4)
  #axis(side = 4, at = pretty(range(optk2$gap), n = 6))
  axis(side = 4, at = seq(0,1,0.2), labels = c("0.0","0.2","0.4","0.6","0.8","1.0"))
  mtext("Gap-statistics", side = 4, line = 2,las = 3, cex = 1.5, col = "#E51718")

  # abline(v = optk3$`Eigen-gap best`, col = "#008B8A", lwd = 2, lty = 4)
  # text(optk3$`Eigen-gap best`, par("usr")[3] + 0.1,"Eigen-gaps", cex = 1.5, col = "#008B8A", adj = -0.05)
  invisible(dev.copy2pdf(file = file.path(fig.path, outFig), width = 5.5, height = 5))

  message("visualization done...")
  if(N.clust > 1) {
    message(paste0("--the imputed optimal cluster number is ", N.clust, " arbitrarily, but it would be better referring to other priori knowledge."))
  }
  #return(list(N.clust = N.clust, CPI = optk1, Gapk = optk2, Eigen = optk3))
  return(list(N.clust = N.clust, CPI = optk1, Gapk = optk2))
}
