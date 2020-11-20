#' @name getSilhouette
#' @title Get silhoutte plot from consensus clustering
#' @description This function visualizes the silhoutte information.
#' @param sil A sil object returned from `getConsensusMOIC()`.
#' @param clust.col A string vector storing colors for annotating each cluster.
#' @param fig.path A string value to indicate the output path for storing the consensus heatmap.
#' @param fig.name A string value to indicate the name of the consensus heatmap.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @return A silhoutte barplot.
#' @export
#' @import cluster
#' @importFrom grDevices dev.copy2pdf
#' @examples # There is no example and please refer to vignette.
getSilhouette <- function(sil       = NULL,
                          clust.col = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627","#023E8A","#9D4EDD"),
                          fig.path  = getwd(),
                          fig.name  = "silhoutte",
                          width     = 5.5,
                          height    = 5) {

  N.clust <- length(unique(sil[,1]))
  colvec <- clust.col[1:N.clust]

  outFig <- paste0(fig.name,".pdf")
  par(bty="o", mgp = c(2.5,0.33,0), mar=c(5.1,2.1,3.1,2.1)+.1, las=1, tcl=-.25)
  plot(sil,
       border = NA,
       col = colvec)
  dev.copy2pdf(file = file.path(fig.path, outFig), width = width, height = height)
}
