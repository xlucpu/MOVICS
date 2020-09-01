#' @name getConsensusMOIC
#' @title Get subtypes from consensus clustering of multiple multi-omics integrative clustering algorithms
#' @description Since this R package integrates 10 mainstream multi-omics clustering algorithms, we borrow the idea of consensus clustering for later integration of the clustering results derived from different algorithms, so as to improve the clustering robustness. The simplest way to run `getConsensusMOIC()` is to pass a list of object returned by `get\%algorithm_name\%()` or by `getMOIC()` with specific argument of `methodslist`.
#' @param moic.res.list A list of object returned by `getMOIC()`.
#' @param distance A string value of distance measurement for hierarchical clustering; 'euclidean' by default.
#' @param linkage A string value of clustering method for hierarchical clustering; 'ward.D' by default.
#' @param mapcolor A string vector for heatmap mapping color.
#' @param clust.col A string vector storing colors for annotating each cluster at the top of heatmap.
#' @param showID A logic value to indicate if showing the sample ID.
#' @param fig.path A string value to indicate the output path for storing the consensus heatmap.
#' @param fig.name A string value to indicate the name of the consensus heatmap.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @return A list contains the following components:
#'
#'         \code{consensus.hm} an object returned by \link[ComplexHeatmap]{pheatmap}
#'
#'         \code{similarity}   a similary matrix for pair-wise samples with entries ranging from 0 to 1
#'
#'         \code{clust.res}    a data.frame storing sample ID and corresponding clusters
#'
#'         \code{clust.dend}   a dendrogram of sample clustering
#'
#'         \code{mo.method}    a string value indicating the method used for multi-omics integrative clustering
#' @importFrom ClassDiscovery distanceMatrix
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom ComplexHeatmap pheatmap
#' @export
#' @examples # There is no example and please refer to vignette.
#' @references Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics, 32(18):2847-2849.
getConsensusMOIC <- function(moic.res.list = NULL,
                             distance      = "euclidean",
                             linkage       = "ward.D",
                             mapcolor      = c("#000004FF", "#56106EFF", "#BB3754FF", "#F98C0AFF", "#FCFFA4FF"),
                             clust.col     = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627"),
                             showID        = F,
                             fig.path      = getwd(),
                             fig.name      = "consensusheatmap",
                             width         = 5.5,
                             height        = 5) {

  if(!is.list(moic.res.list)) {stop("a list of clust.res returned by getMOIC() should be provided.")}

  if(length(moic.res.list) < 2) {stop("at least two multi-omics clustering methods should be included.")}

  moic.res <- lapply(moic.res.list, function(x) x$clust.res)

  if(!length(unique(sapply(moic.res, function(x) {length(unique(x$clust))}))) == 1) {
    stop("clustering number mismatched across different algorithms!")
  }

  N.clust <- length(unique(moic.res[[1]]$clust))
  colvec <- clust.col[1:N.clust]
  names(colvec) <- paste0("CS",1:N.clust)

  fixed.sam <- moic.res[[1]]$samID

  consensus.matrix <- matrix(0, nrow=length(fixed.sam), ncol=length(fixed.sam))
  for (i in 1:length(moic.res)) {
    res <- moic.res[[i]]
    grpInfo <- res$clust; names(grpInfo) <- res$samID

    ans = as.data.frame(matrix(0, nrow=length(fixed.sam), ncol=length(fixed.sam)))
    rownames(ans) = colnames(ans) <- fixed.sam
    unique.grp = unique(grpInfo)

    for (k in 1:length(unique.grp)){
      grpK.members = names(grpInfo)[grpInfo==unique.grp[k]]
      ans[grpK.members, grpK.members] = 1
    }
    consensus.matrix <- consensus.matrix + as.matrix(ans)
  }

  similarity.matrix <- as.data.frame(consensus.matrix/length(moic.res))
  # hcs <- hclust(vegdist(as.matrix(1-similarity.matrix), method = "jaccard"), "ward.D")
  hcs <- hclust(ClassDiscovery::distanceMatrix(as.matrix(1-similarity.matrix), distance), linkage)
  coca.moic <- cutree(hcs, N.clust)

  clustres <- data.frame(samID = fixed.sam,
                         clust = as.numeric(coca.moic),
                         row.names = fixed.sam,
                         stringsAsFactors = F)
  #clustres <- clustres[order(clustres$clust),]

  annCol <- data.frame("Subtype" = paste0("CS", clustres[fixed.sam,"clust"]), # consensus multi-omics integrative clustering
                       row.names = fixed.sam)
  annColors <- list("Subtype" = colvec)

  # save to pdf
  outFig <- paste0(fig.name,".pdf")
  hm <- ComplexHeatmap::pheatmap(mat               = as.matrix(similarity.matrix),
                                 cluster_rows      = hcs,
                                 cluster_cols      = hcs,
                                 border_color      = NA,
                                 show_colnames     = showID,
                                 show_rownames     = showID,
                                 annotation_col    = annCol,
                                 annotation_colors = annColors,
                                 legend_breaks     = c(0,0.2,0.4,0.6,0.8,1),
                                 legend_labels     = c(0,0.2,0.4,0.6,0.8,1),
                                 color             = grDevices::colorRampPalette(mapcolor)(64))

  pdf(file.path(fig.path, outFig), width = width, height = height)
  draw(hm)
  invisible(dev.off())

  # print to screen
  draw(hm)

  return(list(consensus.hm = hm, similarity.matrix = similarity.matrix, clust.res = clustres, clust.dend = hcs, mo.method = "consensusMOIC"))
}
