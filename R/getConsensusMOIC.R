#' @name getConsensusMOIC
#' @title Get subtypes from consensus clustering of multiple multi-omics integrative clustering algorithms
#' @description Since this R package integrates 10 mainstream multi-omics clustering algorithms, we borrow the idea of consensus clustering for later integration of the clustering results derived from different algorithms, so as to improve the clustering robustness. The simplest way to run `getConsensusMOIC()` is to pass a list of object returned by `get\%algorithm_name\%()` or by `getMOIC()` with specific argument of `methodslist`.
#' @param moic.res.list A list of object returned by `getMOIC()`.
#' @param distance A string value of distance measurement for hierarchical clustering; 'euclidean' by default.
#' @param linkage A string value of clustering method for hierarchical clustering; 'ward.D' by default.
#' @param mapcolor A string vector for heatmap mapping color.
#' @param clust.col A string vector storing colors for annotating each cluster at the top of heatmap.
#' @param showID A logic value to indicate if showing the sample ID; FALSE by default.
#' @param fig.path A string value to indicate the output path for storing the consensus heatmap.
#' @param fig.name A string value to indicate the name of the consensus heatmap.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @return A consensus heatmap  and a list contains the following components:
#'
#'         \code{consensus.hm} an object returned by \link[ComplexHeatmap]{pheatmap}
#'
#'         \code{similarity}   a similary matrix for pair-wise samples with entries ranging from 0 to 1
#'
#'         \code{sil}          a silhouette object that can be further passed to \link[MOVICS]{getSilhouette}
#'
#'         \code{clust.res}    a data.frame storing sample ID and corresponding clusters
#'
#'         \code{clust.dend}   a dendrogram of sample clustering
#'
#'         \code{mo.method}    a string value indicating the method used for multi-omics integrative clustering
#' @importFrom ClassDiscovery distanceMatrix
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom ComplexHeatmap pheatmap
#' @import cluster
#' @export
#' @examples # There is no example and please refer to vignette.
#' @references Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics, 32(18):2847-2849.
getConsensusMOIC <- function(moic.res.list = NULL,
                             distance      = "euclidean",
                             linkage       = "ward.D",
                             mapcolor      = c("#000004FF", "#56106EFF", "#BB3754FF", "#F98C0AFF", "#FCFFA4FF"),
                             clust.col     = c("#2EC4B6", "#E71D36", "#FF9F1C", "#BDD5EA", "#FFA5AB", "#011627","#023E8A","#9D4EDD"),
                             showID        = FALSE,
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
                         stringsAsFactors = FALSE)
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

  # silhoutte
  similarity_matrix <- as.matrix(similarity.matrix)
  similarity_matrix <- (similarity_matrix+t(similarity_matrix))/2
  diag(similarity_matrix) <- 0
  normalize <- function(X) X / rowSums(X)
  similarity_matrix <- normalize(similarity_matrix)

  cluster_id <- 1:N.clust
  sil <- matrix(NA, length(coca.moic), 3, dimnames = list(names(coca.moic), c("cluster","neighbor","sil_width")))
  for(j in 1:N.clust) {
    index <- (coca.moic == cluster_id[j])
    Nj <- sum(index)
    sil[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index, drop = FALSE], 2,
                          function(r) tapply(r, coca.moic[!index], mean)))
    maxC <- apply(dindex, 2, which.max)
    sil[index,"neighbor"] <- cluster_id[-j][maxC]
    s.i <- if(Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj - 1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i) / pmax(b.i, a.i), 0)
    } else {0}
    sil[index,"sil_width"] <- s.i
  }
  attr(sil, "Ordered") <- FALSE
  class(sil) <- "silhouette"

  return(list(consensus.hm = hm, similarity.matrix = similarity.matrix, sil = sil, clust.res = clustres, clust.dend = hcs, mo.method = "consensusMOIC"))
}
