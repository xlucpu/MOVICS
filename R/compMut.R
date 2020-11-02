#' @name compMut
#' @title Comparison of mutational frequency
#' @description This function is used to compare mutational frequency among different multi-omics integerative clusters to test the independency between subtypes and mutational status. An oncoprint will be also generated with significant mutations.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param mut.matrix A binary matrix storing binary mutation data with entries of 0 and 1 only.
#' @param test.method A string value to indicate statistical method for independency testing. Allowed values contain c('fisher', 'chisq'); fisher by default.
#' @param freq.cutoff A numeric value to indicate the frequency cutoff for mutation data. Specifically, only features that mutated in over than such proportion would be included in testing; 0.05 by default.
#' @param p.adj.method A string value to indicate the correction method for multiple comparision. Allowed values contain c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'); BH by default.
#' @param tab.name A string value to indicate the name of the output table.
#' @param res.path A string value to indicate the path for saving the table.
#' @param doWord A logic value to indicate if transformating the .txt outfile to a .docx WORD file (.txt file will be also kept); TRUE by default.
#' @param doPlot A logic value to indicate if generating oncoprint; TRUE by default.
#' @param innerclust A logic value to indicate if perform clustering within each subtype; TRUE by default.
#' @param fig.path A string value to indicate the output path for storing the oncoprint.
#' @param fig.name A string value to indicate the name of the oncoprint.
#' @param annCol A data.frame storing annotation information for samples.
#' @param annColors A list of string vectors for colors matched with annCol.
#' @param clust.col A string vector storing colors for annotating each subtype.
#' @param mut.col A string vector to indicate the mutation color for oncoprint.
#' @param bg.col A string vector to indicate the background color for oncoprint.
#' @param p.cutoff A numeric value to indicate the nominal p value cutoff for significant mutations shown in the oncoprint; 0.05 by default.
#' @param p.adj.cutoff A numeric value to indicate the adjusted p value cutoff for significant mutations shown in the oncoprint; 0.05 by default.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @export
#' @return A figure of mutational oncoprint (.pdf) if \code{doPlot = TRUE}, a data.frame storing the difference of mutational frequency among different subtypes and a corresponding table in WORD format if \code{doWord = TRUE}.
#' @importFrom dplyr %>%
#' @importFrom grid grid.rect gpar
#' @importFrom officer read_docx body_add_par body_add_table body_add_par
#' @importFrom ComplexHeatmap HeatmapAnnotation oncoPrint draw
#' @references Gu Z, Eils R, Schlesner M (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics, 32(18):2847–2849.
#' @examples # There is no example and please refer to vignette.
compMut <- function(moic.res     = NULL,
                    mut.matrix   = NULL,
                    freq.cutoff  = 0.05,
                    test.method  = "fisher",
                    p.adj.method = "BH",
                    doWord       = TRUE,
                    doPlot       = TRUE,
                    innerclust   = TRUE,
                    res.path     = getwd(),
                    tab.name     = NULL,
                    fig.path     = getwd(),
                    fig.name     = NULL,
                    annCol       = NULL,
                    annColors    = NULL,
                    mut.col      = "#21498D",
                    bg.col       = "#dcddde",
                    p.cutoff     = 0.05,
                    p.adj.cutoff = 0.05,
                    clust.col    = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                    width        = 8,
                    height       = 4) {

  if(!is.element(test.method, c("fisher","chisq"))) {
    stop("test.method for independency can be one of fisher or chisq.\n")
  }

  if(!is.element(p.adj.method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))) {
    stop("p.adj.method can be one of holm, hochberg, hommel, bonferroni, BH, BY, fdr.\n")
  }

  create.anntrack <- function(samples=NULL, subtype=NULL, typesToPlot=NULL) {
    if (is.null(samples)) {stop("samples can not be NULL!")}
    if (length(samples)!=length(subtype)) {stop("samples and subtype do not have equal length!")}

    if (is.null(typesToPlot)) {
      typesToPlot <- levels(factor(unique(subtype)))
    }
    names(subtype) <- samples
    return(list(subtype=subtype, typesToPlot=typesToPlot))
  }

  createMutSubtype <- function(indata=NULL, samples=NULL, genename=NULL) {
    if (!is.element(genename, rownames(indata))) { stop (paste(genename, "not found in indata!", sep=" ")) }

    comsam <- intersect(colnames(indata), samples)

    out <- rep("Not Available", times=length(samples))
    names(out) <- samples
    out[comsam] <- as.character(indata[genename, comsam])
    out[is.na(out)] <- "Not Available"
    out[out == "0"] <- "Normal"
    out[out == "1"] <- "Mutated"

    return(out)
  }

  clust.res <- moic.res$clust.res
  clust.res <- clust.res[order(clust.res$clust,decreasing = FALSE), , drop = FALSE]
  comsam <- intersect(clust.res$samID, colnames(mut.matrix))
  clust.res <- clust.res[comsam, , drop = FALSE]
  mut.matrix <- mut.matrix[,comsam]
  n.moic <- length(unique(clust.res$clust))

  # check data
  if(length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message("--",paste0((nrow(moic.res$clust.res)-length(comsam))," samples mismatched from current subtypes."))
  }

  ans <- rep(paste0("CS",1:n.moic),as.numeric(table(clust.res$clust))); names(ans) <- clust.res$samID

  genelist <- rownames(mut.matrix[rowSums(mut.matrix) > freq.cutoff * nrow(clust.res),])
  binarymut <- as.data.frame(matrix(0, nrow=nrow(clust.res), ncol=length(genelist)))
  rownames(binarymut) <- names(ans)
  colnames(binarymut) <- genelist
  for (k in 1:length(genelist)) {
    res <- create.anntrack(samples=names(ans), subtype=createMutSubtype(mut.matrix, names(ans), genelist[k]))
    binarymut[, genelist[k]] <- res$subtype
  }

  out <- matrix(0, nrow=length(genelist), ncol=n.moic + 1)
  colnames(out) <- c(levels(factor(ans)), "pvalue")
  rownames(out) <- paste(genelist, "Mutated", sep="_")
  for (k in 1:length(genelist)) {
    genek <- genelist[k]
    x <- ans
    y <- binarymut[names(x), genek, drop = TRUE]; y <- as.character(y); names(y) <- names(x)
    tmp <- setdiff(names(x), names(y)[y=="Not Available"])
    x <- x[tmp]
    y <- y[tmp]
    res <- table(y, x)
    if (!all(colnames(res)==colnames(out)[1:(ncol(out)-1)])) {stop(paste("colnames mismatch for ", k, sep=""))}
    pct <- paste0("(",format(round(res["Mutated", ]/as.numeric(table(clust.res$clust)) * 100, 1), digits = 3),"%)")
    freqpct <- paste(res["Mutated", ], pct, sep = " ")
    out[k, 1:(ncol(out)-1)] <- freqpct

    if(test.method == "fisher") {
      out[k, "pvalue"] <- as.numeric(fisher.test(x, y, workspace = 1e8)$p.value)
    } else {
      out[k, "pvalue"] <- as.numeric(chisq.test(x, y)$p.value)
    }
  }
  out <- as.data.frame(out)
  out$pvalue <- formatC(as.numeric(as.character(out$pvalue)), format = "e", digits = 2)
  out$padj <- formatC(p.adjust(as.numeric(out$pvalue), method = p.adj.method), format = "e", digits = 2)
  if(is.null(tab.name)) {
    outFile <- "Independent test between subtype and mutation.txt"
  } else {
    outFile <- paste0(tab.name,".txt")
  }
  tmb <- rowSums(mut.matrix[genelist,])
  pct <- paste0("(",format(round(tmb/ncol(mut.matrix) * 100, 1), digits = 1),"%)")
  freqpct <- paste(tmb, pct, sep = " ")
  out <- cbind.data.frame(data.frame("Gene (Mutated)" = gsub("_Mutated","",rownames(out)),
                                     "TMB" = freqpct,
                                     check.names = FALSE),
                          out)
  rownames(out) <- NULL
  # message("show heading...")
  # print(head(out))
  write.table(out, file.path(res.path, outFile), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  # generate WORD format
  if(doWord){
    comtable <- out
    title_name <- paste0("Table *. ", gsub(".txt", "", outFile, fixed = TRUE))
    mynote <- "Note: ..."

    my_doc <- officer::read_docx()
    my_doc %>%
      officer::body_add_par(value = title_name, style = "table title") %>%
      officer::body_add_table(value = comtable, style = "table_template") %>%
      officer::body_add_par(value = mynote) %>%
      print(target = file.path(res.path,paste0("TABLE ", gsub(".txt", "", outFile, fixed = TRUE), ".docx")))
  }

  # generate mutation oncoprint
  if(doPlot) {
    if(is.null(fig.name)) {
      outFig <- paste0("oncoprint for mutations with frequency over than " , freq.cutoff * 100, " pct.pdf")
    } else {
      outFig <- paste0(fig.name,".pdf")
    }

    sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, decreasing = FALSE), "samID"]
    colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
    names(colvec) <- paste0("CS",unique(moic.res$clust.res$clust))
    if(!is.null(annCol) & !is.null(annColors)) {
      annCol <- annCol[sam.order, , drop = FALSE]
      annCol$Subtype <- paste0("CS",moic.res$clust.res[sam.order,"clust"])
      annColors[["Subtype"]] <- colvec
    } else {
      annCol <- data.frame("Subtype" = paste0("CS",moic.res$clust.res[sam.order,"clust"]),
                           row.names = sam.order)
      annColors <- list("Subtype" = colvec)
    }
    sig.mut <- as.character(out[which(as.numeric(out$pvalue) < p.cutoff & as.numeric(out$padj) < p.adj.cutoff), "Gene (Mutated)"])
    onco_dat <- t(binarymut[rownames(annCol), sig.mut, drop = FALSE])
    onco_dat[onco_dat == "Normal"] <- ""; onco_dat <- as.data.frame(onco_dat)

    alter_fun = list(
      background = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = bg.col, col = NA))
      },
      Mutated = function(x, y, w, h) {
        grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mut.col, col = NA))
      })
    col = c("Mutated" = mut.col)

    if(innerclust) {
      sam.reorder <- c()
      for (i in 1:n.moic) {
        sam <- moic.res$clust.res[which(moic.res$clust.res$clust == i), "samID"]
        tmp <- quiet(ComplexHeatmap::oncoPrint(onco_dat[,sam], get_type = function(x) x,
                                               alter_fun = alter_fun, col = col,
                                               remove_empty_columns = FALSE,
                                               show_pct = FALSE,
                                               bottom_annotation = NULL,
                                               top_annotation = NULL,
                                               show_heatmap_legend = FALSE))
        sam.reorder <- c(sam.reorder, sam[tmp@column_order])
      }

      my_annotation = ComplexHeatmap::HeatmapAnnotation(df = annCol[sam.reorder, , drop = FALSE], col = annColors)
      p <- quiet(ComplexHeatmap::oncoPrint(onco_dat[,sam.reorder], get_type = function(x) x,
                                           alter_fun = alter_fun, col = col,
                                           remove_empty_columns = FALSE,
                                           column_order = sam.reorder,
                                           show_pct = TRUE,
                                           bottom_annotation = my_annotation,
                                           top_annotation = NULL,
                                           show_heatmap_legend = FALSE))
    } else {
      my_annotation = ComplexHeatmap::HeatmapAnnotation(df = annCol, col = annColors)
      p <- quiet(ComplexHeatmap::oncoPrint(onco_dat, get_type = function(x) x,
                                           alter_fun = alter_fun, col = col,
                                           remove_empty_columns = FALSE,
                                           column_order = colnames(onco_dat),
                                           show_pct = TRUE,
                                           bottom_annotation = my_annotation,
                                           top_annotation = NULL,
                                           show_heatmap_legend = FALSE))
    }

    # save to pdf
    pdf(file.path(fig.path, outFig), width = width, height = height)
    draw(p)
    invisible(dev.off())

    # print to screen
    draw(p)
  }
  return(out)
}
