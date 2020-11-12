#' @name twoclassdeseq2
#' @title Run two class comparison
#' @description Two class differential expression analysis using DESeq2 algorithm.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param countsTable A matrix of RNA-Seq raw count data with rows for genes and columns for samples.
#' @param prefix A string value to indicate the prefix of output file.
#' @param overwt A logic value to indicate if to overwrite existing results; FALSE by default.
#' @param sort.p A logic value to indicate if to sort adjusted p value for output table; TRUE by default.
#' @param verbose A logic value to indicate if to only output id, log2fc, pvalue, and padj; TRUE by default.
#' @param res.path A string value to indicate the path for saving the results.
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol, 15(12):550-558.
#' @export
#' @return Several .txt files storing differential expression analysis results by DESeq2
#' @keywords internal
#' @examples # There is no example and please refer to vignette.
twoclassdeseq2 <- function(moic.res    = NULL,
                           countsTable = NULL,
                           prefix      = NULL,
                           overwt      = FALSE,
                           sort.p      = TRUE,
                           verbose     = TRUE,
                           res.path    = getwd()) {

  # Create comparison list for differential expression analysis between two classes.
  createList  <- function(moic.res = NULL) {

    mo.method <- moic.res$mo.method
    moic.res <- moic.res$clust.res
    tumorsam <- moic.res$samID
    sampleList = list()
    treatsamList =list()
    treatnameList <- c()
    ctrlnameList <- c()

    n.moic <- length(unique(moic.res$clust))

    for (i in 1:n.moic) {
      sampleList[[i]] <- tumorsam
      treatsamList[[i]] = intersect(tumorsam, moic.res[which(moic.res$clust == i), "samID"])
      treatnameList[i] <- paste0("CS",i)
      ctrlnameList[i] <- "Others"
    }

    return(list(sampleList, treatsamList, treatnameList, ctrlnameList, mo.method))
  }

  complist <- createList(moic.res = moic.res)

  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  mo.method <- complist[[5]]
  allsamples <- colnames(countsTable)

  options(warn = 1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]

    compname <- paste(treatname, "_vs_", ctrlname, sep = "")
    tmp <- rep("others", times = length(allsamples))
    names(tmp) <- allsamples
    tmp[samples] <- "control"
    tmp[treatsam] <- "treatment"

    if(!is.null(prefix)) {
      outfile <- file.path(res.path, paste(mo.method, "_", prefix, "_deseq2_test_result.", compname, ".txt", sep = ""))
    } else {
      outfile <- file.path(res.path, paste(mo.method ,"_deseq2_test_result.", compname, ".txt", sep = ""))
    }
    if (file.exists(outfile) & (overwt == FALSE)) {
      cat(paste0("deseq2 of ",compname, " exists and skipped...\n"))
      next
    }

    saminfo <- data.frame("Type" = as.factor(tmp[samples]),
                          "SampleID" = samples,
                          stringsAsFactors = FALSE)

    cts <- countsTable[,samples]
    coldata <- saminfo[samples,]

    dds <- quiet(DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                                colData = coldata,
                                                design = as.formula("~ Type")))

    dds$Type <- relevel(dds$Type,ref = "control")

    dds <- quiet(DESeq2::DESeq(dds))
    res <- DESeq2::results(dds, contrast = c("Type","treatment","control"))

    if(sort.p) {
      resData <- as.data.frame(res[order(res$padj),])
    } else {
      resData <- as.data.frame(res)
    }
    resData$id <- rownames(resData)
    resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    colnames(resData) <- c("id","baseMean","log2fc","lfcSE","stat","pvalue","padj")
    resData$fc <- 2^resData$log2fc
    if(verbose) {
      resData <- resData[,c("id","fc","log2fc","pvalue","padj")]
    } else {
      resData <- resData[,c("id","fc","log2fc","baseMean","lfcSE","stat","pvalue","padj")]
    }

    write.table(resData, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
    cat(paste0("deseq2 of ",compname, " done...\n"))
  }
  options(warn = 0)
}

#' twoclassedger
#' @title Run two class comparison
#' @description Two class differential expression analysis using edgeR algorithm.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param countsTable A matrix of RNA-Seq raw count data with rows for genes and columns for samples.
#' @param prefix A string value to indicate the prefix of output file.
#' @param overwt A logic value to indicate if to overwrite existing results; FALSE by default.
#' @param sort.p A logic value to indicate if to sort adjusted p value for output table; TRUE by default.
#' @param verbose A logic value to indicate if to only output id, log2fc, pvalue, and padj; TRUE by default.
#' @param res.path A string value to indicate the path for saving the results.
#' @references Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1):139-140.
#'
#' McCarthy DJ, Chen Y, Smyth GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Res. 40(10):4288-4297.
#' @export
#' @return Several .txt files storing differential expression analysis results by edgeR
#' @keywords internal
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT topTags
#' @examples # There is no example and please refer to vignette.
twoclassedger <- function(moic.res    = NULL,
                          countsTable = NULL,
                          prefix      = NULL,
                          overwt      = FALSE,
                          sort.p      = TRUE,
                          verbose     = TRUE,
                          res.path    = getwd()) {

  # Create comparison list for differential expression analysis between two classes.
  createList  <- function(moic.res = NULL) {

    mo.method <- moic.res$mo.method
    moic.res <- moic.res$clust.res
    tumorsam <- moic.res$samID
    sampleList <- list()
    treatsamList <-list()
    treatnameList <- c()
    ctrlnameList <- c()

    n.moic <- length(unique(moic.res$clust))

    for (i in 1:n.moic) {
      sampleList[[i]] <- tumorsam
      treatsamList[[i]] <- intersect(tumorsam, moic.res[which(moic.res$clust == i), "samID"])
      treatnameList[i] <- paste0("CS",i)
      ctrlnameList[i] <- "Others"
    }

    return(list(sampleList, treatsamList, treatnameList, ctrlnameList, mo.method))
  }

  complist <- createList(moic.res = moic.res)

  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  mo.method <- complist[[5]]
  allsamples <- colnames(countsTable)

  options(warn = 1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]

    compname <- paste(treatname, "_vs_", ctrlname, sep="")
    tmp <- rep("others", times = length(allsamples))
    names(tmp) <- allsamples
    tmp[samples] <- "control"
    tmp[treatsam] <- "treatment"

    if(!is.null(prefix)) {
      outfile <- file.path(res.path, paste(mo.method, "_", prefix, "_edger_test_result.", compname, ".txt", sep = ""))
    } else {
      outfile <- file.path(res.path, paste(mo.method,"_edger_test_result.", compname, ".txt", sep = ""))
    }
    if (file.exists(outfile) & (overwt==FALSE)) {
      cat(paste0("edger of ",compname, " exists and skipped...\n"))
      next
    }

    saminfo <- data.frame("Type" = tmp[samples],
                          "SampleID" = samples,
                          stringsAsFactors = FALSE)

    group = factor(saminfo$Type,levels = c("control","treatment"))

    design <- model.matrix(~ group)
    rownames(design) <- samples

    y <- edgeR::DGEList(counts = countsTable[,samples],group = saminfo$Type)
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmFit(y, design)
    lrt <- edgeR::glmLRT(fit)
    ordered_tags <- edgeR::topTags(lrt, n = 100000)
    allDiff <- ordered_tags$table
    allDiff <- allDiff[is.na(allDiff$FDR) == FALSE,]
    diff <- allDiff

    diff$id <- rownames(diff)
    resData <- diff[,c("id","logFC","logCPM","LR","PValue","FDR")]
    colnames(resData) <- c("id","log2fc","logCPM","LR","pvalue","padj")
    resData$fc <- 2^resData$log2fc

    if(sort.p) {
      resData <- as.data.frame(resData[order(resData$padj),])
    } else {
      resData <- as.data.frame(resData)
    }
    if(verbose) {
      resData <- resData[,c("id","fc","log2fc","pvalue","padj")]
    } else {
      resData <- resData[,c("id","fc","log2fc","logCPM","LR","pvalue","padj")]
    }
    write.table(resData, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
    cat(paste0("edger of ",compname, " done...\n"))
  }
  options(warn = 0)
}

#' twoclasslimma
#' @title Run two class comparison
#' @description Two class differential expression analysis using limma algorithm.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param norm.expr A matrix of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.
#' @param prefix A string value to indicate the prefix of output file.
#' @param overwt A logic value to indicate if to overwrite existing results; FALSE by default.
#' @param sort.p A logic value to indicate if to sort adjusted p value for output table; TRUE by default.
#' @param verbose A logic value to indicate if to only output id, log2fc, pvalue, and padj; TRUE by default.
#' @param res.path A string value to indicate the path for saving the results.
#' @references Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res, 43(7):e47.
#' @export
#' @return Several .txt files storing differential expression analysis results by limma
#' @keywords internal
#' @importFrom limma lmFit makeContrasts eBayes topTable
#' @examples # There is no example and please refer to vignette.
twoclasslimma <- function(moic.res  = NULL,
                          norm.expr = NULL,
                          prefix    = NULL,
                          overwt    = FALSE,
                          sort.p    = TRUE,
                          verbose   = TRUE,
                          res.path  = getwd()) {

  # Create comparison list for differential expression analysis between two classes.
  createList  <- function(moic.res = NULL) {

    mo.method <- moic.res$mo.method
    moic.res <- moic.res$clust.res
    tumorsam <- moic.res$samID
    sampleList <- list()
    treatsamList <- list()
    treatnameList <- c()
    ctrlnameList <- c()

    n.moic <- length(unique(moic.res$clust))

    for (i in 1:n.moic) {
      sampleList[[i]] <- tumorsam
      treatsamList[[i]] <- intersect(tumorsam, moic.res[which(moic.res$clust == i), "samID"])
      treatnameList[i] <- paste0("CS",i)
      ctrlnameList[i] <- "Others"
    }

    return(list(sampleList, treatsamList, treatnameList, ctrlnameList, mo.method))
  }

  complist <- createList(moic.res = moic.res)

  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  mo.method <- complist[[5]]
  allsamples <- colnames(norm.expr)

  # log transformation
  if(max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    gset <- norm.expr
  }
  if(max(norm.expr) >= 25 & min(norm.expr) >= 0){
    message("--log2 transformation done for expression data.")
    gset <- log2(norm.expr + 1)
  }

  options(warn = 1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]

    compname <- paste(treatname, "_vs_", ctrlname, sep="")
    tmp <- rep("others", times = length(allsamples))
    names(tmp) <- allsamples
    tmp[samples] <- "control"
    tmp[treatsam] <- "treatment"

    if(!is.null(prefix)) {
      outfile <- file.path(res.path, paste(mo.method, "_", prefix, "_limma_test_result.", compname, ".txt", sep = ""))
    } else {
      outfile <- file.path(res.path, paste(mo.method, "_limma_test_result.", compname, ".txt", sep = ""))
    }

    if (file.exists(outfile) & (overwt == FALSE)) {
      cat(paste0("limma of ",compname, " exists and skipped...\n"))
      next
    }

    pd <- data.frame(Samples = names(tmp),
                     Group = as.character(tmp),
                     stringsAsFactors = FALSE)

    design <-model.matrix(~ -1 + factor(pd$Group, levels = c("treatment","control")))
    colnames(design) <- c("treatment","control")

    fit <- limma::lmFit(gset, design = design);
    contrastsMatrix <- limma::makeContrasts(treatment - control, levels = c("treatment", "control"))
    fit2 <- limma::contrasts.fit(fit, contrasts = contrastsMatrix)
    fit2 <- limma::eBayes(fit2, 0.01)
    resData <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", number = 100000)
    resData <- as.data.frame(subset(resData, select=c("logFC","t","B","P.Value","adj.P.Val")))
    resData$id <- rownames(resData)
    colnames(resData) <- c("log2fc","t","B","pvalue","padj","id")
    resData$fc <- 2^resData$log2fc

    if(sort.p) {
      resData <- resData[order(resData$padj),]
    } else {
      resData <- as.data.frame(resData)
    }
    if(verbose) {
      resData <- resData[,c("id","fc","log2fc","pvalue","padj")]
    } else {
      resData <- resData[,c("id","fc","log2fc","t","B","pvalue","padj")]
    }
    write.table(resData, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
    cat(paste0("limma of ",compname, " done...\n"))
  }
  options(warn = 0)
}

#' runDEA
#' @title Run differential expression analysis
#' @description Using choosen algorithm to run differential expression analysis between two classes identified by multi-omics clustering process.
#' @param dea.method A string value to indicate the algorithm for differential expression analysis. Allowed value contains c('deseq2', 'edger', 'limma'). The former two require RNA-Seq raw count data and the last one requires normalized expression data (FPKM or TPM without log2 transformation is recommended).
#' @param expr A matrix of expression data.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param prefix A string value to indicate the prefix of output file.
#' @param overwt A logic value to indicate if to overwrite existing results; FALSE by default.
#' @param sort.p A logic value to indicate if to sort adjusted p value for output table; TRUE by default.
#' @param verbose A logic value to indicate if to only output id, log2fc, pvalue, and padj; TRUE by default.
#' @param res.path A string value to indicate the path for saving the results.
#' @export
#' @return Several .txt files storing differential expression analysis results by specified algorithm
#' @importFrom limma lmFit makeContrasts eBayes topTable
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT topTags
#' @references Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol, 15(12):550-558.
#'
#' Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26(1):139-140.
#'
#' McCarthy DJ, Chen Y, Smyth GK (2012). Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Res. 40(10):4288-4297.
#'
#' Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res, 43(7):e47.
#'
#' @examples # There is no example and please refer to vignette.
runDEA <- function(dea.method = c("deseq2", "edger", "limma"),
                   expr       = NULL,
                   moic.res   = NULL,
                   prefix     = NULL,
                   overwt     = FALSE,
                   sort.p     = TRUE,
                   verbose    = TRUE,
                   res.path   = getwd()) {

  if(!is.element(dea.method, c("deseq2", "edger", "limma"))) {
    stop("unsupported algorithm: dea.method should be one of 'deseq2', 'edger', or 'limma'.")
  }

  method <- dea.method[1] # deseq2 by default

  comsam <- intersect(moic.res$clust.res$samID, colnames(expr))
  # check data
  if(length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(moic.res$clust.res)-length(comsam))," samples mismatched from current subtypes."))
  }

  moic.res$clust.res <- moic.res$clust.res[comsam,,drop = FALSE]
  expr <- expr[,comsam]

  rundea <- switch(method,
                   "deseq2" = twoclassdeseq2,
                   "edger"  = twoclassedger,
                   "limma"  = twoclasslimma)
  if(method %in% c("deseq2", "edger")) {
    message(paste0("--you choose ", method," and please make sure an RNA-Seq count data was provided."))
    rundea(moic.res    = moic.res,
           countsTable = expr,
           prefix      = prefix,
           overwt      = overwt,
           sort.p      = sort.p,
           verbose     = verbose,
           res.path    = res.path)
  } else {
    message(paste0("--you choose ", method," and please make sure a microarray profile or a normalized expression data [FPKM or TPM without log2 transformation is recommended] was provided."))
    rundea(moic.res    = moic.res,
           norm.expr   = expr,
           prefix      = prefix,
           overwt      = overwt,
           sort.p      = sort.p,
           verbose     = verbose,
           res.path    = res.path)
  }
}
