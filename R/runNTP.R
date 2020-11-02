#' @name runNTP
#' @title Run nearest template prediction
#' @description Using Nearest Template Prediction (NTP) based on predefined templates derived from current identified subtypes to assign potential subtype label on external cohort.
#' @param expr A numeric matrix with row features and sample columns; data is recommended to be z-scored.
#' @param templates A data frame with at least two columns; class (coerced to factor) and probe (coerced to character).
#' @param scaleFlag A logic value to indicate if the expression data should be further scaled. TRUE by default.
#' @param centerFlag A logic value to indicate if the expression data should be further centered. TRUE by default.
#' @param nPerm An integer value to indicate the permutations for p-value estimation.
#' @param distance A string value to indicate the distance measurement. Allowed values contain c('cosine', 'pearson', 'spearman', 'kendall'); "cosine" by default.
#' @param seed An integer value for p-value reproducibility.
#' @param verbose A logic value to indicate whether console messages are to be displayed; TRUE by default.
#' @param doPlot A logic value to indicate whether to produce prediction heatmap; FALSE by default.
#' @param fig.path A string value to indicate the output path for storing the nearest template prediction heatmap.
#' @param fig.name A string value to indicate the name of the nearest template prediction heatmap.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @return A figure of predictive heatmap by NTP (.pdf) and a list with the following components:
#'
#'         \code{ntp.res}    a data.frame storing the results of nearest template prediction (see \link[CMScaller]{ntp}).
#'
#'         \code{clust.res}  similar to `clust.res` returned by `getMOIC()` or `get%algorithm_name%` or `getConsensusMOIC()`.
#'
#'         \code{mo.method}  a string value indicating the method used for prediction.
#' @export
#' @importFrom CMScaller ntp subHeatmap
#' @importFrom grDevices dev.copy2pdf
#' @examples # There is no example and please refer to vignette.
#' @references Hoshida, Y. (2010). Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE 5, e15543.
runNTP <- function(expr       = NULL,
                   templates  = NULL,
                   scaleFlag  = TRUE,
                   centerFlag = TRUE,
                   nPerm      = 1000,
                   distance   = "cosine",
                   seed       = 123456,
                   verbose    = TRUE,
                   doPlot     = FALSE,
                   fig.path   = getwd(),
                   fig.name   = "ntpheatmap",
                   width      = 5,
                   height     = 5) {

  # message("Using up- or down-regulated biomarkers (templates) are highly recommended.\n")

  if(!is.element(distance, c("cosine", "pearson", "spearman", "kendall"))) {
    stop("the argument of distance should be one of cosine, pearson, spearman, or kendall.")
  }

  com_feat <- intersect(rownames(expr), templates$probe)
  message(paste0("--original template has ",nrow(templates), " biomarkers and ", length(com_feat)," are matched in external expression profile."))
  expr <- expr[com_feat, , drop = FALSE]
  templates <- templates[which(templates$probe %in% com_feat), , drop = FALSE]

  if(is.element(0,as.numeric(table(templates$class)))) {
    stop("at least one class has no probes/genes matched in template file!")
  }

  emat <- t(scale(t(expr), scale = scaleFlag, center = centerFlag))
  if(doPlot) {
    outFig <- paste0(fig.name,".pdf")
    ntp.res <- ntp(emat      = emat,
                   templates = templates,
                   doPlot    = doPlot,
                   nPerm     = nPerm,
                   distance  = distance,
                   nCores    = 1,
                   seed      = seed,
                   verbose   = verbose)
    invisible(dev.copy2pdf(file = file.path(fig.path, outFig), width = width, height = height))

  } else {
    ntp.res <- ntp(emat      = emat,
                   templates = templates,
                   doPlot    = doPlot,
                   nPerm     = nPerm,
                   distance  = distance,
                   nCores    = 1,
                   seed      = seed,
                   verbose   = verbose)
  }

  ntp.res[,setdiff(colnames(ntp.res),"prediction")] <- round(ntp.res[,setdiff(colnames(ntp.res),"prediction")], 4)
  ex.moic.res <- data.frame(samID = rownames(ntp.res),
                            clust = gsub("CS","",ntp.res$prediction),
                            row.names = rownames(ntp.res),
                            stringsAsFactors = FALSE)

  return(list(ntp.res = ntp.res, clust.res = ex.moic.res, mo.method = "NTP"))
}
