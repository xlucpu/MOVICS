#' @name runMarker
#' @title Run identification of unique biomarkers
#' @description This function aims to identify uniquely and significantly expressed (overexpressed or downexpressed) biomarkers for each subtype identified by multi-omics integrative clustering algorithms. Top markers will be chosen to generate a template so as to run nearest template prediction for subtype verification.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param dea.method A string value to indicate the algorithm for differential expression analysis. Allowed value contains c('deseq2', 'edger', 'limma').
#' @param dat.path A string value to indicate the path for saving the files of differential expression analysis.
#' @param res.path A string value to indicate the path for saving the results for identifying subtype-specific markers.
#' @param p.cutoff A numeric value to indicate the nominal p value for identifying significant markers; pvalue < 0.05 by default.
#' @param p.adj.cutoff A numeric value to indicate the adjusted p value for identifying significant markers; padj < 0.05 by default.
#' @param dirct A string value to indicate the direction of identifying significant marker. Allowed values contain c('up', 'down'); `up` means up-regulated marker, and `down` means down-regulated marker.
#' @param n.marker A integer value to indicate how many top markers sorted by log2fc should be identified for each subtype; 200 by default.
#' @param norm.expr A matrix of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.
#' @param prefix A string value to indicate the prefix of differential expression file (use for searching files).
#' @param doplot A logic value to indicate if generating heatmap by using subtype-specific markers. TRUE by default.
#' @param annCol A data.frame storing annotation information for samples.
#' @param annColors A list of string vectors for colors matched with annCol.
#' @param clust.col A string vector storing colors for annotating each subtype at the top of heatmap.
#' @param halfwidth A numeric vector to assign marginal cutoff for truncating values in data; 3 by default.
#' @param centerFlag A logical vector to indicate if expression data should be centered; TRUE by default.
#' @param scaleFlag A logical vector to indicate if expression data should be scaled; TRUE by default.
#' @param color A string vector storing colors for heatmap.
#' @param fig.path A string value to indicate the output path for storing the marker heatmap.
#' @param fig.name A string value to indicate the name of the marker heatmap.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @param show_rownames A logic value to indicate if showing rownames (feature names) in heatmap; FALSE by default.
#' @param show_colnames A logic value to indicate if showing colnames (sample ID) in heatmap; FALSE by default.
#' @param ... Additional parameters pass to \link[ComplexHeatmap]{pheatmap}.
#' @return A figure of subtype-specific marker heatmap (.pdf) if \code{doPlot = TRUE} and a list with the following components:
#'
#'         \code{unqlist}   a string vector storing the unique marker across all subtypes.
#'
#'         \code{templates} a data.frame storing the the template information for nearest template prediction, which is used for verification in external cohort.
#'
#'         \code{dirct}     a string value indicating the direction for identifying subtype-specific markers.
#'
#'         \code{heatmap}   a complexheatmap object.
#' @export
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom ComplexHeatmap pheatmap draw
#' @references Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics, 32(18):2847-2849.
#' @examples # There is no example and please refer to vignette.
runMarker <- function(moic.res      = NULL,
                      dea.method    = c("deseq2", "edger", "limma"),
                      prefix        = NULL,
                      dat.path      = getwd(),
                      res.path      = getwd(),
                      p.cutoff      = 0.05,
                      p.adj.cutoff  = 0.05,
                      dirct         = "up",
                      n.marker      = 200,
                      doplot        = TRUE,
                      norm.expr     = NULL,
                      annCol        = NULL,
                      annColors     = NULL,
                      clust.col     = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                      halfwidth     = 3,
                      centerFlag    = TRUE,
                      scaleFlag     = TRUE,
                      show_rownames = FALSE,
                      show_colnames = FALSE,
                      color         = c("#5bc0eb", "black", "#ECE700"),
                      fig.path      = getwd(),
                      fig.name      = NULL,
                      width         = 8,
                      height        = 8,
                      ...) {

  n.moic <- length(unique(moic.res$clust.res$clust))
  mo.method <- moic.res$mo.method
  DEpattern <- paste(mo.method, "_", ifelse(is.null(prefix),"",paste0(prefix,"_")), dea.method, ".*._vs_Others.txt$", sep="")
  DEfiles <- dir(dat.path,pattern = DEpattern)

  if(length(DEfiles) == 0) {stop("no DEfiles!")}
  if(length(DEfiles) != n.moic) {stop("not all the multi-omics clusters have DEfile!")}
  if (!is.element(dirct, c("up", "down"))) {stop( "dirct type error! Allowed value contains c('up', 'down').") }

  if(dirct == "up") { outlabel <- "unique_upexpr_marker.txt" }
  if(dirct == "down") { outlabel <- "unique_downexpr_marker.txt" }

  genelist <- c()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header=TRUE, row.names=NULL, sep="\t", quote="", stringsAsFactors=FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]),]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]

    #rownames(DEres) <- toupper(rownames(DEres))
    if (dirct == "up") {
      genelist <- c( genelist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc > 0, ]) )
    }
    if (dirct == "down") {
      genelist <- c( genelist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc < 0, ]) )
    }
  }
  unqlist <- setdiff(genelist,genelist[duplicated(genelist)])

  marker <- list()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header=TRUE, row.names=NULL, sep="\t", quote="", stringsAsFactors=FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]),]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]

    #rownames(DEres) <- toupper(rownames(DEres))
    if(dirct == "up") {
      outk <- intersect( unqlist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc > 0, ]) )
      outk <- DEres[outk,]
      outk <- outk[order(outk$log2fc, decreasing = TRUE),]

      if(nrow(outk) > n.marker) {
        marker[[filek]] <- outk[1:n.marker,]
      } else {
        marker[[filek]] <- outk
      }
      marker$dirct <- "up"
    }
    if(dirct == "down") {
      outk <- intersect( unqlist, rownames(DEres[!is.na(DEres$padj) & DEres$pvalue < p.cutoff & DEres$padj < p.adj.cutoff & !is.na(DEres$log2fc) & DEres$log2fc < 0, ]) )
      outk <- DEres[outk,]
      outk <- outk[order(outk$log2fc, decreasing = FALSE),]

      if(nrow(outk) > n.marker) {
        marker[[filek]] <- outk[1:n.marker,]
      } else {
        marker[[filek]] <- outk
      }
      marker$dirct <- "down"
    }
    # write file
    write.table(outk, file=file.path(res.path, paste(gsub("_vs_Others.txt","", filek, fixed = TRUE), outlabel, sep = "_")), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }

  # generate templates for nearest template prediction
  templates <- NULL
  for (filek in DEfiles) {
    tmp <- data.frame(probe = rownames(marker[[filek]]),
                      class = sub("_vs_Others.txt","",sub(".*.result.","",filek)),
                      dirct = marker$dirct,
                      stringsAsFactors = FALSE)
    templates <- rbind.data.frame(templates, tmp, stringsAsFactors = FALSE)
  }
  write.table(templates, file=file.path(res.path, paste0(mo.method,"_",dea.method,"_",dirct,"regulated_marker_templates.txt")), row.names = FALSE, sep = "\t", quote = FALSE)

  # generate heatmap with subtype-specific markers
  if(doplot) {
    if(is.null(norm.expr)) {
      stop("please provide a matrix or data.frame of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.")
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

    if(is.null(fig.name)) {
      outFig <- paste0("markerheatmap_using_",dirct,"regulated_genes.pdf")
    } else {
      outFig <- paste0(fig.name, "_using_", dirct,"regulated_genes.pdf")
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

    if(max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 0)) {
      message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
      gset <- norm.expr
    }
    if(max(norm.expr) >= 25 & min(norm.expr) >= 0){
      message("--log2 transformation done for expression data.")
      gset <- log2(norm.expr + 1)
    }

    standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=TRUE, scaleFlag=TRUE) {
      outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
      if (!is.null(halfwidth)) {
        outdata[outdata>halfwidth]=halfwidth
        outdata[outdata<(-halfwidth)]= -halfwidth
      }
      return(outdata)
    }

    plotdata <- standarize.fun(gset[intersect(templates$probe, rownames(gset)), sam.order],
                               halfwidth = halfwidth,
                               centerFlag = centerFlag,
                               scaleFlag = scaleFlag)

    if(!is.null(annCol) & !is.null(annColors)) {
      for (i in names(annColors)) {
        if(is.function(annColors[[i]])) {
          annColors[[i]] <- annColors[[i]](pretty(range(annCol[,i]),n = 64)) # transformat colorRamp2 function to color vector
        }
      }
    }

    hm <- ComplexHeatmap::pheatmap(mat               = plotdata,
                                   border_color      = NA,
                                   cluster_cols      = FALSE,
                                   cluster_rows      = FALSE,
                                   annotation_col    = annCol,
                                   annotation_colors = annColors,
                                   legend_breaks     = pretty(c(-halfwidth,halfwidth)),
                                   legend_labels     = pretty(c(-halfwidth,halfwidth)),
                                   show_rownames     = show_rownames,
                                   show_colnames     = show_colnames,
                                   treeheight_col    = 0,
                                   treeheight_row    = 0,
                                   color             = grDevices::colorRampPalette(color)(64),
                                   ...)

    # save to pdf
    pdf(file.path(fig.path, outFig), width = width, height = height)
    draw(hm,annotation_legend_side = "left",heatmap_legend_side = "left")
    invisible(dev.off())

    # print to screen
    draw(hm,annotation_legend_side = "left",heatmap_legend_side = "left")
    return(list(unqlist = unqlist, templates = templates, dirct = dirct, heatmap = hm))
  } else {
    return(list(unqlist = unqlist, templates = templates, dirct = dirct))
  }
}
