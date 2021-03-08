#' @name runGSEA
#' @title Run identification of unique functional pathways
#' @description Use gene set enrichment analysis to identify subtype-specific (overexpressed or downexpressed) functional pathways for each subtype identified by multi-omics clustering algorithms.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param dea.method A string value to indicate the algorithm for differential expression analysis. Allowed value contains c('deseq2', 'edger', 'limma').
#' @param norm.expr A matrix of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.
#' @param prefix A string value to indicate the prefix of differential expression file (use for searching files).
#' @param dat.path A string value to indicate the path for saving the files of differential expression analysis.
#' @param res.path A string value to indicate the path for saving the results for identifying subtype-specific functional pathways.
#' @param dirct A string value to indicate the direction of identifying significant pathway. Allowed values contain c('up', 'down'); `up` means up-regulated pathway, and `down` means down-regulated pathway; "up" by default.
#' @param n.path A integer value to indicate how many top pathways sorted by NES should be identified for each subtypes; 10 by default.
#' @param msigdb.path A string value to indicate ABSOULUTE PATH/NAME of MSigDB file (GMT file with gene symbols) downloaded from \url{https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H}.
#' @param nPerm A integer value to indicate the number of permutations; 1000 by default and 10000 will be better for reproducibility.
#' @param minGSSize A integer value to indicate minimal size of each geneSet for analyzing; 10 by default.
#' @param maxGSSize A integer value to indicate maximal size of each geneSet for analyzing; 500 by default.
#' @param p.cutoff A numeric value to indicate the nominal p value for identifying significant pathways; pvalue < 0.05 by default.
#' @param p.adj.cutoff A numeric value to indicate the adjusted p value for identifying significant pathways; padj < 0.05 by default.
#' @param gsva.method A string value to indicate the method to employ in the estimation of gene-set enrichment scores per sample. By default this is set to gsva (Hänzelmann et al, 2013) and other options are ssgsea (Barbie et al, 2009), zscore (Lee et al, 2008) or plage (Tomfohr et al, 2005). The latter two standardize first expression profiles into z-scores over the samples and, in the case of zscore, it combines them together as their sum divided by the square-root of the size of the gene set, while in the case of plage they are used to calculate the singular value decomposition (SVD) over the genes in the gene set and use the coefficients of the first right-singular vector as pathway activity profile.
#' @param norm.method A string value to indicate how to calculate subtype-specific pathway enrichment scores. Allowed values contain c('mean', 'median'); mean by default.
#' @param clust.col A string vector storing colors for annotating each subtype at the top of heatmap.
#' @param color A string vector storing colors for heatmap.
#' @param fig.path A string value to indicate the output path for storing the pathway heatmap.
#' @param fig.name A string value to indicate the name of the pathway heatmap.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @return A figure of subtype-specific pathway heatmap (.pdf) and a list with the following components:
#'
#'         \code{gsea.list}  a list storing gsea object returned by \link[clusterProfiler]{GSEA} for each subtype.
#'
#'         \code{raw.es}     a data.frame storing raw enrichment score of identified subtype-specific pathways by using specified \code{gsva.method}.
#'
#'         \code{scaled.es}  a data.frame storing z-scored enrichment score of identified subtype-specific pathways by using specified \code{gsva.method}.
#'
#'         \code{grouped.es} a data.frame storing grouped enrichment score (mean or median value among each subtype) by using specified \code{norm.method}.
#'
#'         \code{heatmap}    a complexheatmap object.
#' @export
#' @importFrom grDevices pdf dev.off colorRampPalette
#' @importFrom ComplexHeatmap pheatmap draw
#' @importFrom GSVA gsva
#' @importFrom clusterProfiler GSEA read.gmt
#' @references Barbie, D.A. et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature, 462(5):108-112.
#'
#' Hänzelmann, S., Castelo, R. and Guinney, J. (2013). GSVA: Gene set variation analysis for microarray and RNA-Seq data. BMC Bioinformatics, 14(1):7.
#'
#' Lee, E. et al. (2008). Inferring pathway activity toward precise disease classification. PLoS Comp Biol, 4(11):e1000217.
#'
#' Tomfohr, J. et al. (2005). Pathway level analysis of gene expression using singular value decomposition. BMC Bioinformatics, 6(1):1-11.
#'
#' Yu G, Wang L, Han Y, He Q (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS, 16(5):284-287.
#'
#' Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics, 32(18):2847–2849.
#' @examples # There is no example and please refer to vignette.
runGSEA <- function(moic.res     = NULL,
                    dea.method   = c("deseq2", "edger", "limma"),
                    norm.expr    = NULL,
                    prefix       = NULL,
                    dat.path     = getwd(),
                    res.path     = getwd(),
                    dirct        = "up",
                    n.path       = 10,
                    msigdb.path  = NULL,
                    nPerm        = 1000,
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    p.cutoff     = 0.05,
                    p.adj.cutoff = 0.05,
                    gsva.method  = "gsva",
                    norm.method  = "mean",
                    clust.col    = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                    color        = NULL,
                    fig.name     = NULL,
                    fig.path     = getwd(),
                    width        = 15,
                    height       = 10) {

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
  mo.method <- moic.res$mo.method
  DEpattern <- paste(mo.method, "_", ifelse(is.null(prefix),"",paste0(prefix,"_")), dea.method, ".*._vs_Others.txt$", sep="")
  DEfiles <- dir(dat.path,pattern = DEpattern)

  if(length(DEfiles) == 0) {stop("no DEfiles!")}
  if(length(DEfiles) != n.moic) {stop("not all multi-omics clusters have DEfile!")}
  if (!is.element(dirct, c("up", "down"))) {stop( "dirct type error! Allowed value contains c('up', 'down').") }
  if (!is.element(gsva.method, c("gsva", "ssgsea", "zscore", "plage"))) {stop( "GSVA method error! Allowed value contains c('gsva', 'ssgsea', 'zscore', 'plage').") }

  if(dirct == "up") { outlabel <- "unique_upexpr_pathway.txt" }
  if(dirct == "down") { outlabel <- "unique_downexpr_pathway.txt" }

  gsea.list <- list()
  gseaidList <- c()
  for (filek in DEfiles) {
    DEres <- read.table(file.path(dat.path, filek), header=TRUE, row.names=NULL, sep="\t", quote="", stringsAsFactors=FALSE)
    DEres <- DEres[!duplicated(DEres[, 1]),]
    DEres <- DEres[!is.na(DEres[, 1]), ]
    rownames(DEres) <- DEres[, 1]
    DEres <- DEres[, -1]

    geneList <- DEres$log2fc; names(geneList) <- rownames(DEres)
    geneList <- sort(geneList,decreasing = TRUE) # ranked gene set

    # run gsea
    msigdb <- try(clusterProfiler::read.gmt(msigdb.path), silent = TRUE)
    if(class(msigdb) == "try-error") {stop("please provide correct ABSOLUTE PATH for MSigDB file.")}

    moic.lab <- paste0("CS",gsub("_vs_Others.txt","",sub(".*CS", "", filek)))
    gsea.list[[moic.lab]] <- suppressWarnings(clusterProfiler::GSEA(geneList     = geneList,
                                                                    TERM2GENE    = msigdb,
                                                                    nPerm        = nPerm,
                                                                    minGSSize    = minGSSize,
                                                                    maxGSSize    = maxGSSize,
                                                                    seed         = TRUE,
                                                                    verbose      = FALSE,
                                                                    pvalueCutoff = 1))
    gsea.dat <- as.data.frame(gsea.list[[moic.lab]])
    write.table(gsea.dat[,setdiff(colnames(gsea.dat),"ID")], file=file.path(res.path, paste(gsub("_vs_Others.txt","", filek, fixed = TRUE), "gsea_all_results.txt", sep = "_")),sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

    if(dirct == "up") {
      gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES > 0 & gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < p.adj.cutoff),]))
    }
    if(dirct == "down") {
      gseaidList <- c(gseaidList, rownames(gsea.dat[which(gsea.dat$NES < 0 & gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < p.adj.cutoff),]))
    }
    unqlist <- setdiff(gseaidList,gseaidList[duplicated(gseaidList)])
  }
  message("GSEA done...")

  GSEApattern <- paste(mo.method, "_", ifelse(is.null(prefix),"",paste0(prefix,"_")), dea.method, ".*._gsea_all_results.txt$", sep="")
  GSEAfiles <- dir(res.path,pattern = GSEApattern)

  pathway <- pathcore <- list()
  pathnum <- c()
  for (filek in GSEAfiles) {
    GSEAres <- read.table(file.path(res.path, filek), header=TRUE, row.names=1, sep="\t", quote="", stringsAsFactors=FALSE)

    if(dirct == "up") {
      outk <- intersect( unqlist, rownames(GSEAres[which(GSEAres$NES > 0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < p.adj.cutoff),]) )
      outk <- GSEAres[outk,]
      outk <- outk[order(outk$NES, decreasing = TRUE),]

      if(nrow(outk) > n.path) {
        pathway[[filek]] <- outk[1:n.path,]
      } else {
        pathway[[filek]] <- outk
      }
      pathnum <- c(pathnum, nrow(pathway[[filek]]))
      pathway$dirct <- "up"

      for (i in rownames(pathway[[filek]])) {
        pathcore[[i]] <- msigdb[which(msigdb[,1] %in% i),"gene"]
      }
    }
    if(dirct=="down") {
      outk <- intersect( unqlist, rownames(GSEAres[which(GSEAres$NES < 0 & GSEAres$pvalue < p.cutoff & GSEAres$p.adjust < p.adj.cutoff),]) )
      outk <- GSEAres[outk,]
      outk <- outk[order(outk$NES, decreasing = FALSE),]

      if(nrow(outk) > n.path) {
        pathway[[filek]] <- outk[1:n.path,]
      } else {
        pathway[[filek]] <- outk
      }
      pathnum <- c(pathnum, nrow(pathway[[filek]]))
      pathway$dirct <- "down"

      for (i in rownames(pathway[[filek]])) {
        pathcore[[i]] <- msigdb[which(msigdb[,1] %in% i),"gene"]
      }
    }
    write.table(outk, file=file.path(res.path, paste(gsub("_gsea_all_results.txt","", filek, fixed = TRUE), outlabel, sep = "_")), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }

  # calculate single sample enrichment scores
  standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=TRUE, scaleFlag=TRUE) {
    outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
    if (!is.null(halfwidth)) {
      outdata[outdata>halfwidth]=halfwidth
      outdata[outdata<(-halfwidth)]= -halfwidth
    }
    return(outdata)
  }

  rowmean <- function(x) {
    return(apply(x, 1, mean))
  }

  rowmedian <- function(x) {
    return(apply(x, 1, median))
  }

  if(max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    gset <- norm.expr
  }
  if(max(norm.expr) >= 25 & min(norm.expr) >= 0){
    message("--log2 transformation done for expression data.")
    gset <- log2(norm.expr + 1)
  }

  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, decreasing = FALSE), "samID"]
  colvec <- clust.col[1:n.moic]
  names(colvec) <- paste0("CS",1:n.moic)
  annCol <- data.frame("Subtype" = paste0("CS",moic.res$clust.res[sam.order,"clust"]),
                       row.names = sam.order,
                       stringsAsFactors = FALSE)
  annColors <- list("Subtype" = colvec)

  es <- GSVA::gsva(expr          = as.matrix(gset[,rownames(annCol), drop = FALSE]),
                   gset.idx.list = pathcore,
                   method        = gsva.method,
                   parallel.sz   = 1)
  es.backup <- es
  es <- standarize.fun(es, halfwidth = 1, centerFlag = TRUE, scaleFlag = TRUE)
  message(gsva.method," done...")

  # calculate subtype-specific pathway enrichment score
  esm <- data.frame(row.names = rownames(es))
  if(norm.method == "mean") {
    for (i in paste0("CS",1:n.moic)) {
      esm <- cbind.data.frame(esm,
                              data.frame(rowmean(es[,rownames(annCol[which(annCol$Subtype == i), , drop = FALSE])])))
    }
  }

  if(norm.method == "median") {
    for (i in paste0("CS",1:n.moic)) {
      esm <- cbind.data.frame(esm,
                              data.frame(rowmedian(es[,rownames(annCol[which(annCol$Subtype == i), , drop = FALSE])])))
    }
  }

  colnames(esm) <- paste0("CS",1:n.moic)

  # generate heatmap
  annRow <- data.frame(Subtype = rep(paste0("CS",1:n.moic), pathnum),
                       row.names = rownames(esm),
                       stringsAsFactors = FALSE)

  if(is.null(color)) {
    mapcolor <- grDevices::colorRampPalette(c("#0000FF", "#8080FF", "#FFFFFF", "#FF8080", "#FF0000"))(64)
  } else {mapcolor <- grDevices::colorRampPalette(color)(64)}
  hm <- ComplexHeatmap::pheatmap(mat                  = as.matrix(esm),
                                 cluster_rows         = FALSE,
                                 cluster_cols         = FALSE,
                                 show_rownames        = TRUE,
                                 show_colnames        = TRUE,
                                 annotation_row       = annRow,
                                 annotation_colors    = annColors,
                                 annotation_names_row = FALSE,
                                 legend               = TRUE,
                                 color                = mapcolor,
                                 border_color         = "black",
                                 legend_breaks        = c(-1,-0.5,0,0.5,1),
                                 legend_labels        = c(-1,-0.5,0,0.5,1),
                                 cellwidth            = 15,
                                 cellheight           = 10)
  # print to screen
  draw(hm, annotation_legend_side = "left",heatmap_legend_side = "left")

  # save to pdf
  if(is.null(fig.name)) {
    outFig <- paste0("gseaheatmap_using_",dirct,"regulated_pathways.pdf")
  } else {
    outFig <- paste0(fig.name, "_using_", dirct,"regulated_pathways.pdf")
  }
  pdf(file.path(fig.path,outFig), width = width, height = height)
  draw(hm, annotation_legend_side = "left",heatmap_legend_side = "left")
  invisible(dev.off())

  message("heatmap done...")

  return(list(gsea.list = gsea.list, raw.es = es.backup, scaled.es = es, grouped.es = esm, heatmap = hm))
}
