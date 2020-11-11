#' @name getStdiz
#' @title Get standardized omics data
#' @description This function prepare standardized data for generating heatmap. Omics data, especially for expression, should be centered or scaled or z-scored (both centered and scaled). Generally, DNA methylation beta matrix and somatic mutation (0 and 1 binary matrix) should not be normalized. This function also provides an argument of `halfwidth` for continuous omics data; such argument is used to truncate the 'extremum' after normalization; specifically, normalized values that exceed the halfwidth boundaries will be replaced by the halfwidth, which is vary useful to map colors in heatmap.
#' @param data A list of data.frame or matrix storing raw multiple omics data with rows for features and columns for samples.
#' @param halfwidth A numeric vector to assign marginal cutoff for truncating values in data; 1 by default.
#' @param centerFlag A logical vector to indicate if each subdata should be centered; TRUE by default.
#' @param scaleFlag A logical vector to indicate if each subdata should be scaled; TRUE by default.
#' @export
#' @return A standardized data.frame containing multi-omics data.
#' @examples # There is no example and please refer to vignette.
getStdiz <- function(data       = NULL,
                     halfwidth  = rep(1, length(data)),
                     centerFlag = rep(TRUE, length(data)),
                     scaleFlag  = rep(TRUE, length(data))) {

  # check data
  if(is.null(names(data))){
    names(data) <- sprintf("dat%s", 1:length(data))
  }

  n_dat <- length(data)
  if(n_dat > 6){
    stop("current version of MOVICS can support up to 6 datasets.")
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }

  outdata <- list()
  for (i in 1:n_dat) {
    tmp <- t(scale(t(data[[i]]), center = centerFlag[i], scale = scaleFlag[i]))
    if (!is.na(halfwidth[i])) {
      tmp[tmp > halfwidth[i]] = halfwidth[i]
      tmp[tmp < (-halfwidth[i])] = -halfwidth[i]
    }
    outdata[[names(data)[i]]] <- tmp
  }

  return(outdata)
}
