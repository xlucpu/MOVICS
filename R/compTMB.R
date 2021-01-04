#' @name compTMB
#' @title Comparsion of total mutation burden
#' @description This function calculates Total Mutation Burden (TMB) compares them among curent subtypes identified from multi-omics integrative clustering algorithms.
#'
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param maf A data frame of MAF file that has been already read with at least 10 columns as following: c('Tumor_Sample_Barcode', 'Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2')
#' @param rmDup A logical value to indicate if removing repeated variants in a particuar sample, mapped to multiple transcripts of same Gene. TRUE by default.
#' @param rmFLAGS A logical value to indicate if removing possible FLAGS. These FLAGS genes are often non-pathogenic and passengers, but are frequently mutated in most of the public exome studies, some of which are fishy. Examples of such genes include TTN, MUC16, etc. FALSE by default.
#' @param nonSyn A string vector to indicate a list of variant claccifications that should be considered as non-synonymous and the rest will be considered synonymous (silent) variants. Default value of NULL uses Variant Classifications with High/Moderate variant consequences, including c('Frame_Shift_Del', 'Frame_Shift_Ins', 'Splice_Site', 'Translation_Start_Site', 'Nonsense_Mutation', 'Nonstop_Mutation', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation'). See details at \url{http://asia.ensembl.org/Help/Glossary?id=535}
#' @param clust.col A string vector storing colors for annotating each Subtype.
#' @param test.method A string value to indicate the method for statistical testing. Allowed values contain c('nonparametric', 'parametric'); nonparametric means two-sample wilcoxon rank sum test for two subtypes and Kruskal-Wallis rank sum test for multiple subtypes; parametric means two-sample t-test when only two subtypes are identified, and anova for multiple subtypes comparison; "nonparametric" by default.
#' @param show.size A logical value to indicate if showing the sample size within each subtype at the top of the figure. TRUE by default.
#' @param fig.name  A string value to indicate the name of the boxviolin plot.
#' @param fig.path A string value to indicate the output path for storing the boxviolin plot.
#' @param exome.size An integer value to indicate the estimation of exome size. 38 by default (see \url{https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0424-2}).
#' @param width A numeric value to indicate the width of boxviolin plot.
#' @param height A numeric value to indicate the height of boxviolin plot.
#' @return A figure of TMB and TiTv distribution (.pdf) and a list with the following components:
#'
#'         \code{TMB.dat}       a data.frame storing the TMB per sample within each subtype.
#'
#'         \code{TMB.median}    a data.frame storing the median of TMB for each subtype.
#'
#'         \code{titv.dat}      a data.frame storing the fraction contributions of TiTv per sample within each subtype.
#'
#'         \code{maf.nonsilent} a data.frame storing the information for non-synonymous mutations.
#'
#'         \code{maf.silent}    a data.frame storing the information for synonymous mutations.
#'
#'         \code{maf.FLAGS}     a data.frame storing the information for FLAGS mutations if\code{rmFLAGS = TRUE}.
#'
#'         \code{FLAGS.count}   a data.frame storing the summarization per FLAGS if\code{rmFLAGS = TRUE}.
#' @export
#' @importFrom maftools read.maf titv
#' @importFrom dplyr %>% group_by summarise mutate n
#' @importFrom grDevices dev.copy2pdf
#' @examples # There is no example and please refer to vignette.
#' @references Mayakonda A, Lin D, Assenov Y, Plass C, Koeffler PH (2018). Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Res, 28(11): 1747-1756.
#'
#' Shyr C, Tarailo-Graovac M, Gottlieb M, Lee JJ, van Karnebeek C, Wasserman WW. (2014). FLAGS, frequently mutated genes in public exomes. BMC Med Genomics, 7(1): 1-14.
#'
#' Chalmers Z R, Connelly C F, Fabrizio D, et al. (2017). Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden. Genome Med, 9(1):34.
compTMB <- function(moic.res    = NULL,
                    maf         = NULL,
                    rmDup       = TRUE,
                    rmFLAGS     = FALSE,
                    nonSyn      = NULL,
                    exome.size  = 38,
                    clust.col   = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                    test.method = "nonparametric",
                    show.size   = TRUE,
                    fig.path    = getwd(),
                    fig.name    = NULL,
                    width       = 6,
                    height      = 6) {

  label <- c("Tumor_Sample_Barcode",
             "Hugo_Symbol",
             "Chromosome",
             "Start_Position",
             "End_Position",
             "Variant_Classification",
             "Variant_Type",
             "Reference_Allele",
             "Tumor_Seq_Allele1",
             "Tumor_Seq_Allele2")

  maf <- as.data.frame(maf)

  # check arguments
  if(!all(is.element(label, colnames(maf)))) {
    stop(paste0("maf data must have the following columns: \n ", paste(label, collapse = "\n "),"\n\nmissing required fields from maf: ", paste(setdiff(label, colnames(maf)), collapse = "\n")))
  }

  if(!is.element(test.method, c("nonparametric","parametric"))) {
    stop("test.method can be one of nonparametric or parametric.")
  }

  # check data
  comsam <- intersect(moic.res$clust.res$samID, unique(maf$Tumor_Sample_Barcode))
  if(length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(moic.res$clust.res)-length(comsam))," samples mismatched from current subtypes."))
  }

  maf <- maf[which(maf$Tumor_Sample_Barcode %in% comsam),]
  clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  n.moic <- length(unique(clust.res$clust))
  colvec <- clust.col[1:n.moic]
  names(colvec) <- paste0("CS",unique(clust.res$clust))
  col.titv <- c("#E64B35CC", "#4DBBD5CC", "#00A087CC", "#3C5488CC", "#F39B7FCC", "#8491B4CC")
  names(col.titv) <- c("C>T", "T>C", "C>A", "C>G", "T>A", "T>G")

  if(rmFLAGS) {
    FLAGS <- c("TTN",
               "MUC16",
               "OBSCN",
               "AHNAK2",
               "SYNE1",
               "FLG",
               "MUC5B",
               "DNAH17",
               "PLEC",
               "DST",
               "SYNE2",
               "NEB",
               "HSPG2",
               "LAMA5",
               "AHNAK",
               "HMCN1",
               "USH2A",
               "DNAH11",
               "MACF1",
               "MUC17",
               "DNAH5",
               "GPR98",
               "FAT1",
               "PKD1",
               "MDN1",
               "RNF213",
               "RYR1",
               "DNAH2",
               "DNAH3",
               "DNAH8",
               "DNAH1",
               "DNAH9",
               "ABCA13",
               "APOB",
               "SRRM2",
               "CUBN",
               "SPTBN5",
               "PKHD1",
               "LRP2",
               "FBN3",
               "CDH23",
               "DNAH10",
               "FAT4",
               "RYR3",
               "PKHD1L1",
               "FAT2",
               "CSMD1",
               "PCNT",
               "COL6A3",
               "FRAS1",
               "FCGBP",
               "DNAH7",
               "RP1L1",
               "PCLO",
               "ZFHX3",
               "COL7A1",
               "LRP1B",
               "FAT3",
               "EPPK1",
               "VPS13C",
               "HRNR",
               "MKI67",
               "MYO15A",
               "STAB1",
               "ZAN",
               "UBR4",
               "VPS13B",
               "LAMA1",
               "XIRP2",
               "BSN",
               "KMT2C",
               "ALMS1",
               "CELSR1",
               "TG",
               "LAMA3",
               "DYNC2H1",
               "KMT2D",
               "BRCA2",
               "CMYA5",
               "SACS",
               "STAB2",
               "AKAP13",
               "UTRN",
               "VWF",
               "VPS13D",
               "ANK3",
               "FREM2",
               "PKD1L1",
               "LAMA2",
               "ABCA7",
               "LRP1",
               "ASPM",
               "MYOM2",
               "PDE4DIP",
               "TACC2",
               "MUC2",
               "TEP1",
               "HELZ2",
               "HERC2",
               "ABCA4")

    if(sum(is.element(FLAGS, unique(maf$Hugo_Symbol))) > 0) {
      maf.FLAGS <- maf[which(maf$Hugo_Symbol %in% FLAGS),]
      maf.rmFLAGS <- maf[-which(maf$Hugo_Symbol %in% FLAGS),]
      count.flags <- maf.FLAGS %>% group_by(Hugo_Symbol) %>% summarise(count = dplyr::n())
      message("--remove possible FLAGS as below:")
      head(count.flags)

      maf.ob <- maftools::read.maf(maf                      = maf.rmFLAGS,
                                   removeDuplicatedVariants = rmDup,
                                   vc_nonSyn                = nonSyn)
    } else {
      maf.ob <- maftools::read.maf(maf                      = maf,
                                   removeDuplicatedVariants = rmDup,
                                   vc_nonSyn                = nonSyn)
    }
  }   else {
      maf.ob <- maftools::read.maf(maf                      = maf,
                                   removeDuplicatedVariants = rmDup,
                                   vc_nonSyn                = nonSyn)
  }

  # classifies Single Nucleotide Variants into Transitions and Transversions
  titv.dat <- maftools::titv(maf = maf.ob, plot = FALSE, useSyn = FALSE)$fraction.contribution %>%
    as.data.frame() %>%
    mutate(Subtype = paste0("CS",clust.res[as.character(.$Tumor_Sample_Barcode),"clust"]))
  titv.dat.backup <- titv.dat
  titv.dat <- split(titv.dat, f = titv.dat$Subtype)

  # extract silent and nonsilent mutations
  maf.silent    <- as.data.frame(maf.ob@maf.silent)
  maf.nonsilent <- as.data.frame(maf.ob@data)

  # extract total mutation burden
  TMB.dat <- as.data.frame(maf.ob@variants.per.sample)
  TMB.dat <- data.frame(samID            = as.character(TMB.dat$Tumor_Sample_Barcode),
                        variants         = as.character(TMB.dat$Variants),
                        TMB              = as.numeric(TMB.dat$Variants)/exome.size,
                        log10TMB         = log10(as.numeric(TMB.dat$Variants)/exome.size + 1),
                        Subtype          = paste0("CS", clust.res[as.character(TMB.dat$Tumor_Sample_Barcode), "clust"]),
                        stringsAsFactors = FALSE)
  TMB.dat <- TMB.dat[order(TMB.dat$Subtype), , drop = FALSE]
  TMB.med <- TMB.dat %>% group_by(Subtype) %>% summarize(median = median(TMB)) %>% as.data.frame()
  TMB.dat <- TMB.dat[order(TMB.dat$Subtype, TMB.dat$TMB, decreasing = FALSE), ]

  # sample order in bottom panel
  sampleorder <- TMB.dat %>% split(.$Subtype) %>% lapply("[[", 1) %>% lapply(., as.character)
  TMB.plot <- split(TMB.dat, as.factor(TMB.dat$Subtype))
  TMB.plot <- lapply(seq_len(length(TMB.plot)), function(i) {
    x = TMB.plot[[i]]
    x = data.frame(x = seq(i - 1, i, length.out = nrow(x)),
                   TMB = x[, "TMB"],
                   Subtype = x[, "Subtype"])
    return(x)
  })
  names(TMB.plot) <- levels(as.factor(TMB.dat$Subtype))

  # prepare titv data
  titv.dat2 <- lapply(TMB.med$Subtype, function(x){
    tmp <- titv.dat[[x]]
    tmp <- tmp[match(sampleorder[[x]], as.character(tmp$Tumor_Sample_Barcode)), ]
    return(tmp)
    if (!all(tmp$Tumor_Sample_Barcode == sampleorder[[x]])){
      stop("inconsistent sample order...")
    }
  })

  names(titv.dat2) <- TMB.med$Subtype
  titv.dat2 <- lapply(titv.dat2, function(x){
    x <- as.data.frame(x)
    rownames(x) <- x$Tumor_Sample_Barcode
    x <- x[, setdiff(colnames(x), c("Tumor_Sample_Barcode","Subtype"))]
    x <- t(x)
    #delete samples without mutation
    if (length(which(colSums(x) == 0)) > 0) {
      x = x[, -which(colSums(x) == 0), drop = FALSE]
    }
    return(x)
  })

  TMB.med$Median_Mutations_log10 <- log10(TMB.med$median + 1)

  # statistical testing
  if(n.moic == 2 & test.method == "nonparametric") {
    statistic <- "wilcox.test"
    TMB.test  <- wilcox.test(TMB.dat$log10TMB ~ TMB.dat$Subtype)$p.value
    cat(paste0("Wilcoxon rank sum test p value = ", formatC(TMB.test, format = "e", digits = 2)))
  }

  if(n.moic == 2 & test.method == "parametric") {
    statistic <- "t.test"
    TMB.test  <- t.test(TMB.dat$log10TMB ~ TMB.dat$Subtype)$p.value
    cat(paste0("Student's t test p value = ", formatC(TMB.test, format = "e", digits = 2)))
  }

  if(n.moic > 2 & test.method == "nonparametric") {
    statistic <- "kruskal.test"
    TMB.test  <- kruskal.test(TMB.dat$log10TMB ~ TMB.dat$Subtype)$p.value
    pairwise.TMB.test <- pairwise.wilcox.test(TMB.dat$log10TMB,TMB.dat$Subtype,p.adjust.method = "BH")
    # pairwise.TMB.test <- dunnTest(log10TMB ~ as.factor(Subtype),
    #                               data = TMB.dat,
    #                               method = "bh")
    cat(paste0("Kruskal-Wallis rank sum test p value = ", formatC(TMB.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
    print(formatC(pairwise.TMB.test$p.value, format = "e", digits = 2))
  }

  if(n.moic > 2 & test.method == "parametric") {
    statistic <- "anova"
    TMB.test  <- summary(aov(TMB.dat$log10TMB ~ TMB.dat$Subtype))[[1]][["Pr(>F)"]][1]
    pairwise.TMB.test <- pairwise.t.test(TMB.dat$log10TMB,TMB.dat$Subtype,p.adjust.method = "BH")
    cat(paste0("One-way anova test p value = ", formatC(TMB.test, format = "e", digits = 2),"\npost-hoc pairwise Student's t test with Benjamini-Hochberg adjustment presents below:\n"))
    print(formatC(pairwise.TMB.test$p.value, format = "e", digits = 2))
  }

  # start illustration
  if(is.null(fig.name)) {
    outFig <- "distribution of TMB and titv.pdf"
  } else {
    outFig <- paste0(fig.name,".pdf")
  }

  # base layout
  n1 <- seq(from = 0.105, to = 0.97-(0.97-0.05)*0.04, length.out = n.moic + 1)
  n2 <- n1[2:length(n1)]
  n  <- data.frame(n1 = n1[1:n.moic], n2 = n2, n3 = 0.05, n4 = 0.2) %>%
    rbind(c(0.05, 0.97, 0.25, 1), ., c(0, 0.1, 0.05, 0.25)) %>% as.matrix()

  opar <- par(no.readonly = TRUE)
  invisible(suppressWarnings(split.screen(n, erase = TRUE)))
  screen(1, new = TRUE)
  par(xpd = TRUE, mar = c(3, 1, 2, 0), oma = c(0, 0, 0, 0),bty = "o", mgp = c(2, 0.5, 0), tcl=-.25)
  y_lims = range(log10(unlist(lapply(TMB.plot, function(x) max(x[, "TMB"], na.rm = TRUE))) + 1))
  y_lims[1] = 0
  y_lims[2] = ceiling(max(y_lims))
  y_at = y_lims[1]:y_lims[2]
  x_top_label <- as.numeric(unlist(lapply(TMB.plot, nrow)))
  plot(NA, NA, xlim = c(0, length(TMB.plot)), ylim = y_lims,
       xlab = NA, ylab = NA, xaxt="n", yaxt = "n", xaxs = "r", yaxs = "r")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#EAE9E9", border = FALSE)
  grid(col = "white", lty = 1, lwd = 1.5) # add grid

  par(new = TRUE, bty="o")
  plot(NA, NA,
       col = "white",
       xlim = c(0, length(TMB.plot)), ylim = y_lims,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")

  text(x = length(TMB.plot)/2, y = y_lims[2] - 0.1, cex = 1.2,
       label = paste0(statistic, " p = ", formatC(TMB.test, digits = 1, format = "e")))

  # add median TMB
  invisible(lapply(seq_len(nrow(TMB.med)), function(i) {
    segments(x0  = i - 1,
             x1  = i,
             y0  = TMB.med[i, "Median_Mutations_log10"],
             y1  = TMB.med[i, "Median_Mutations_log10"],
             col = "#2b2d42",
             lwd = 2)}))

  # add scatter TMB
  invisible(lapply(seq_len(length(TMB.plot)), function(i){
    tmp = TMB.plot[[i]]
    points(tmp$x, log10(tmp$TMB + 1), pch = 16, cex = 0.5, col = colvec[i])
  }))

  # modify axis
  axis(side = 1, at = seq(0.5, length(TMB.plot) - 0.5, 1), labels = names(TMB.plot),
       las = 1, tick = TRUE, line = 0)
  axis(side = 2, at = y_at, las = 2, line = 0, tick = TRUE, labels = y_at)
  if(show.size) {
    axis(side = 3, at = seq(0.5, length(TMB.plot) - 0.5, 1), labels = paste0("n = ",x_top_label),
         tick = TRUE, line = 0, cex.axis = 0.9)
  }
  mtext(text = bquote("log"[10]~"(TMB + 1)"), side = 2, line = 1.15, cex = 1.1)

  # add titv
  invisible(lapply(seq_len(length(titv.dat2)), function(i){
    screen(i + 1, new = TRUE)
    par(xpd = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), bty = "o")
    tmp <- titv.dat2[[i]]
    barplot(tmp, col = col.titv[rownames(tmp)], names.arg = rep("", ncol(tmp)),
            xaxs = "i", yaxs = "i",
            axes = FALSE, space = 0, border = NA, lwd = 1.2)
    box()
  }))

  # add legend
  screen(n.moic + 2, new = TRUE)
  par(xpd = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), bty = "n")
  plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = NA, ylab = NA)
  legend("center",
         fill   = col.titv,
         legend = names(col.titv),
         border = NA,
         bty    = "n",
         cex    = 0.8)
  close.screen(all.screens = TRUE)
  invisible(dev.copy2pdf(file = file.path(fig.path, outFig), width = width, height = height))

  if(rmFLAGS) {
    return(list(TMB.dat = TMB.dat, TMB.median = TMB.med, titv.dat = titv.dat.backup, maf.nonsilent = maf.nonsilent, maf.silent = maf.silent, maf.FLAGS = maf.FLAGS, FLAGS.count = as.data.frame(count.flags)))
  } else {
    return(list(TMB.dat = TMB.dat, TMB.median = TMB.med, titv.dat = titv.dat.backup, maf.nonsilent = maf.nonsilent, maf.silent = maf.silent))
  }
}
