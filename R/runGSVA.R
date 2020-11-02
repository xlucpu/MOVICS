#' @name runGSVA
#' @title Run gene set variation analysis
#' @description Use gene set variation analysis to calculate enrichment score of each sample in each subtype based on given gene set list of interest.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param norm.expr A matrix of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.
#' @param gset.gmt.path A string value to indicate ABSOULUTE PATH/NAME of gene sets of interest stored as GMT format \url{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}.
#' @param annCol A data.frame storing annotation information for samples.
#' @param annColors A list of string vectors for colors matched with annCol.
#' @param clust.col A string vector storing colors for annotating each subtype at the top of heatmap.
#' @param halfwidth A numeric value to assign marginal cutoff for truncating enrichment scores; 1 by default.
#' @param centerFlag A logical vector to indicate if enrichment scores should be centered; TRUE by default.
#' @param scaleFlag A logical vector to indicate if enrichment scores should be scaled; TRUE by default.
#' @param distance A string value of distance measurement for hierarchical clustering; 'euclidean' by default.
#' @param linkage A string value of clustering method for hierarchical clustering; 'ward.D' by default.
#' @param show_rownames A logic value to indicate if showing rownames (feature names) in heatmap; TRUE by default.
#' @param show_colnames A logic value to indicate if showing colnames (sample ID) in heatmap; FALSE by default.
#' @param color A string vector storing colors for heatmap.
#' @param fig.path A string value to indicate the output path for storing the enrichment heatmap.
#' @param fig.name A string value to indicate the name of the enrichment heatmap.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @param ... Additional parameters pass to \link[ComplexHeatmap]{pheatmap}.
#'
#' @return A figure of enrichment heatmap (.pdf) and a list with the following components:
#'
#'         \code{gset.list}  a list storing gene sets information converted from GMT format by \link[clusterProfiler]{read.gmt}.
#'
#'         \code{raw.es}     a data.frame storing raw enrichment score based on given gene sets of interest by using specified \code{gsva.method}.
#'
#'         \code{scaled.es}  a data.frame storing z-scored enrichment score based on given gene sets of interest by using specified \code{gsva.method}.
#'
#' @export
#' @importFrom ClassDiscovery distanceMatrix
#' @importFrom clusterProfiler read.gmt
#' @importFrom GSVA gsva
#' @importFrom ComplexHeatmap pheatmap draw ht_opt
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @references Barbie, D.A. et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(5):108-112.
#'
#' Hänzelmann, S., Castelo, R. and Guinney, J. (2013). GSVA: Gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14(1):7.
#'
#' Lee, E. et al. (2008). Inferring pathway activity toward precise disease classification. PLoS Comp Biol, 4(11):e1000217.
#'
#' Tomfohr, J. et al. (2005). Pathway level analysis of gene expression using singular value decomposition. BMC Bioinformatics, 6(1):1-11.
#'
#' Yu G, Wang L, Han Y, He Q (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5):284-287.
#' @examples # There is no example and please refer to vignette.
runGSVA <- function(moic.res      = NULL,
                    norm.expr     = NULL,
                    gset.gmt.path = NULL,
                    gsva.method   = "gsva",
                    centerFlag    = TRUE,
                    scaleFlag     = TRUE,
                    halfwidth     = 1,
                    annCol        = NULL,
                    annColors     = NULL,
                    clust.col     = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                    distance      = "euclidean",
                    linkage       = "ward.D",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    color         = c("#366A9B", "#4E98DE", "#DDDDDD", "#FBCFA7", "#F79C4A"),
                    fig.path      = getwd(),
                    fig.name      = NULL,
                    width         = 8,
                    height        = 8,
                    ...) {
  
  # standardize function
  standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=TRUE, scaleFlag=TRUE) {
    outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata>halfwidth]=halfwidth
      outdata[outdata<(-halfwidth)]= -halfwidth
    }
    return(outdata)
  }
  
  # check data
  comsam <- intersect(moic.res$clust.res$samID, colnames(norm.expr))
  if(length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(moic.res$clust.res)-length(comsam))," samples mismatched from current subtypes."))
  }
  
  moic.res$clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  norm.expr <- norm.expr[,comsam]
  n.moic <- length(unique(moic.res$clust.res$clust))
  
  # load gene set data and convert gmt to data.frame
  gset <- try(clusterProfiler::read.gmt(gset.gmt.path), silent = TRUE)
  if(class(gset) == "try-error") {stop("please provide correct ABSOLUTE PATH for gene sets of interest.")}
  
  # convert data.frame to list
  term <- unique(gset[,1])
  gset.list <- list()
  for (i in term) {
    gset.list[[i]] <- gset[which(gset[,1] == i),2]
  }
  
  # calculate gene set enrichment scores
  if(max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
  }
  if(max(norm.expr) >= 25 & min(norm.expr) >= 0){
    message("--log2 transformation done for expression data.")
    norm.expr <- log2(norm.expr + 1)
  }
  
  es <- GSVA::gsva(expr          = as.matrix(norm.expr),
                   gset.idx.list = gset.list,
                   method        = gsva.method,
                   parallel.sz   = 1)
  es.backup <- es
  es <- standarize.fun(es, halfwidth = halfwidth, centerFlag = centerFlag, scaleFlag = scaleFlag)
  message(gsva.method," done...")
  
  if(is.null(fig.name)) {
    outFig <- paste0("enrichment_heatmap_using_", gsva.method, ".pdf")
  } else {
    outFig <- paste0(fig.name, "_", gsva.method, ".pdf")
  }
  
  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, decreasing = FALSE), "samID"]
  colvec <- clust.col[1:n.moic]
  names(colvec) <- paste0("CS",1:n.moic)
  if(!is.null(annCol) & !is.null(annColors)) {
    annCol <- annCol[sam.order, , drop = FALSE]
    annCol$Subtype <- paste0("CS",moic.res$clust.res[sam.order,"clust"])
    annColors[["Subtype"]] <- colvec
  } else {
    annCol <- data.frame("Subtype" = paste0("CS",moic.res$clust.res[sam.order,"clust"]),
                         row.names = sam.order,
                         stringsAsFactors = FALSE)
    annColors <- list("Subtype" = colvec)
  }
  
  if(!is.null(annCol) & !is.null(annColors)) {
    for (i in names(annColors)) {
      if(is.function(annColors[[i]])) {
        annColors[[i]] <- annColors[[i]](pretty(range(annCol[,i]),n = 64)) # transformat colorRamp2 function to color vector
      }
    }
  }
  
  ht_opt$message = FALSE
  if(is.null(distance) | is.null(linkage)) {
    hcg <- FALSE
  } else {
    hcg <- hclust(ClassDiscovery::distanceMatrix(t(as.matrix(es[,sam.order])), distance), linkage)
  }
  hm <- ComplexHeatmap::pheatmap(mat               = es[,sam.order],
                                 border_color      = NA,
                                 cluster_cols      = FALSE,
                                 cluster_rows      = hcg,
                                 annotation_col    = annCol,
                                 annotation_colors = annColors,
                                 show_rownames     = show_rownames,
                                 show_colnames     = show_colnames,
                                 color             = grDevices::colorRampPalette(color)(64),
                                 ...)
  
  # save to pdf
  pdf(file.path(fig.path, outFig), width = width, height = height)
  draw(hm)
  invisible(dev.off())
  
  # print to screen
  draw(hm)
  
  return(list(gset.list = gset.list, raw.es = es.backup, scaled.es = es))
}