#' @name compFGA
#' @title Comparison of fraction genome altered
#' @description This function calculates Fraction Genome Altered (FGA), Fraction Genome Gained (FGG), and Fraction Genome Lost (FGL) seperately, and compares them among curent subtypes identified from multi-omics integrative clustering algorithms.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param segment A data frame containing segmented copy number and columns must exactly include the following elements: c('sample','chrom','start','end','value'). Column of `value` should be segments value when \code{iscopynumber = FALSE} but copy-number value when \code{iscopynumber = TRUE}. Copy-number will be converted to segments by log2(copy-number/2).
#' @param iscopynumber A logical value to indicate if the fifth column of segment input is copy-number. If segment file derived from CNV calling provides copy number instead of segment_mean value, this argument must be switched to TRUE. FALSE by default.
#' @param cnathreshold A numeric value to indicate the cutoff for identifying copy-number gain or loss. 0.2 by default.
#' @param test.method A string value to indicate the method for statistical testing. Allowed values contain c('nonparametric', 'parametric'); nonparametric means two-sample wilcoxon rank sum test for two subtypes and Kruskal-Wallis rank sum test for multiple subtypes; parametric means two-sample t-test when only two subtypes are identified, and anova for multiple subtypes comparison; "nonparametric" by default.
#' @param barcolor A string vector to indicate the mapping color for bars of FGA, FGG and FGL.
#' @param clust.col A string vector storing colors for each subtype.
#' @param fig.path A string value to indicate the output path for storing the barplot.
#' @param fig.name A string value to indicate the name of the barplot.
#' @param width A numeric value to indicate the width of barplot.
#' @param height A numeric value to indicate the height of barplot.
#' @return A list contains the following components:
#'
#'         \code{summary}           a table summarizing the measurements of FGA, FGG, and FGL per sample
#'
#'         \code{FGA.p.value}       a nominal p value quantifying the difference of FGA among current subtypes
#'
#'         \code{pairwise.FGA.test} a pairwise BH adjusted p value matrix for multiple comparisons of FGA if more than 2 subtypes were identified
#'
#'         \code{FGG.p.value}       a nominal p value quantifying the difference of FGG among current subtypes
#'
#'         \code{pairwise.FGG.test} a pairwise BH adjusted p value matrix for multiple comparisons of FGG if more than 2 subtypes were identified
#'
#'         \code{FGL.p.value}       a nominal p value quantifying the difference of FGL among current subtypes
#'
#'         \code{pairwise.FGL.test} a pairwise BH adjusted p value matrix for multiple comparisons of FGL if more than 2 subtypes were identified
#'
#'         \code{test.method}       a string value indicating the statistical testing method to calculate p values
#'
#' @export
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr group_by summarize %>%
#' @importFrom ggpubr stat_compare_means
#' @examples # There is no example and please refer to vignette.
#' @references Cerami E, Gao J, Dogrusoz U, et al. (2012). The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data. Cancer Discov, 2(5):401-404.
#'
#' Gao J, Aksoy B A, Dogrusoz U, et al. (2013). Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Sci Signal, 6(269):pl1-pl1.
compFGA <- function(moic.res     = NULL,
                    segment      = NULL,
                    iscopynumber = FALSE,
                    cnathreshold = 0.2,
                    test.method  = "nonparametric",
                    barcolor     = c("#008B8A", "#F2042C", "#21498D"),
                    clust.col    = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                    fig.path     = getwd(),
                    fig.name     = NULL,
                    width        = 8,
                    height       = 4) {

  # check arguments
  if(!all(is.element(c("sample","chrom","start","end","value"), colnames(segment)))) {
    stop("segment data must have the following columns: sample, chrom, start, end, value.")
  }

  if(iscopynumber) {
    segment$value <- log2(segment$value/2) # convert copy-number to segment-mean value
  }

  comsam <- intersect(moic.res$clust.res$samID, unique(segment[,1]))
  # check data
  if(length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(moic.res$clust.res)-length(comsam))," samples mismatched from current subtypes."))
  }

  if(!is.element(test.method, c("nonparametric","parametric"))) {
    stop("test.method can be one of nonparametric or parametric.")
  }

  clust.res <- moic.res$clust.res[comsam, , drop = FALSE]
  segment <- segment[which(segment$sample %in% comsam),]
  n.moic <- length(unique(clust.res$clust))

  # data process
  segment$bases <- segment$end - segment$start

  # calculate FGA, FGG and FGL
  display.progress = function (index, totalN, breakN=20) {
    if ( index %% ceiling(totalN/breakN)  == 0  ) {
      cat(paste(round(index*100/totalN), "% ", sep = ""))
    }
  }
  std <- function(x, na.rm = TRUE) {
    if(na.rm) {
      x <- as.numeric(na.omit(x))
      sd(x)/sqrt(length(x))
    } else {sd(x)/sqrt(length(x))}
  }

  outTab <- data.frame()
  for (i in 1:length(unique(segment$sample))) {
    display.progress(index = i, totalN = length(unique(segment$sample)))

    tmp <- segment[segment$sample == names(table(segment$sample))[i],]

    if (length(tmp[abs(tmp$value) > cnathreshold,"bases"][6]) == 0) {
      FGA = 0
    } else {
      FGA = sum(tmp[abs(tmp$value) > cnathreshold,"bases"]) / sum(tmp[,"bases"])
    }

    if (length(tmp[tmp$value > cnathreshold,"bases"][6]) == 0) {
      FGG = 0
    } else {
      FGG = sum(tmp[tmp$value > cnathreshold,"bases"]) / sum(tmp[,"bases"])
    }

    if (length(tmp[tmp$value < (-cnathreshold),"bases"][6]) == 0) {
      FGL = 0
    } else {
      FGL = sum(tmp[tmp$value < (-cnathreshold),"bases"]) / sum(tmp[,"bases"])
    }

    tmp <- data.frame(samID            = names(table(segment$sample))[i],
                      FGA              = FGA,
                      FGG              = FGG,
                      FGL              = FGL,
                      stringsAsFactors = FALSE)
    outTab <- rbind.data.frame(outTab, tmp, stringsAsFactors = FALSE)
  }
  outTab$Subtype <- paste0("CS", clust.res[outTab$samID, "clust"])

  # calculate mean and se
  summaryFGA  <- outTab %>% group_by(Subtype) %>% dplyr::summarize(mean = mean(FGA, na.rm = TRUE), se = std(FGA, na.rm = TRUE))
  summaryFGG  <- outTab %>% group_by(Subtype) %>% dplyr::summarize(mean = mean(FGG, na.rm = TRUE), se = std(FGG, na.rm = TRUE))
  summaryFGL  <- outTab %>% group_by(Subtype) %>% dplyr::summarize(mean = mean(FGL, na.rm = TRUE), se = std(FGL, na.rm = TRUE))
  summaryFGGL <- data.frame(rbind.data.frame(summaryFGG,summaryFGL),class = rep(c("FGG","FGL"),c(nrow(summaryFGG),nrow(summaryFGL))),stringsAsFactors = FALSE)

  # statistical testing
  # generate boxviolin plot with statistical testing
  if(n.moic == 2 & test.method == "nonparametric") {
    statistic <- "wilcox.test"
    FGA.test  <- wilcox.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test  <- wilcox.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test  <- wilcox.test(outTab$FGL ~ outTab$Subtype)$p.value
  }

  if(n.moic == 2 & test.method == "parametric") {
    statistic <- "t.test"
    FGA.test  <- t.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test  <- t.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test  <- t.test(outTab$FGL ~ outTab$Subtype)$p.value
  }

  if(n.moic > 2 & test.method == "nonparametric") {
    statistic <- "kruskal.test"
    FGA.test  <- kruskal.test(outTab$FGA ~ outTab$Subtype)$p.value
    FGG.test  <- kruskal.test(outTab$FGG ~ outTab$Subtype)$p.value
    FGL.test  <- kruskal.test(outTab$FGL ~ outTab$Subtype)$p.value
    pairwise.FGA.test <- pairwise.wilcox.test(outTab$FGA,outTab$Subtype,p.adjust.method = "BH")
    pairwise.FGG.test <- pairwise.wilcox.test(outTab$FGG,outTab$Subtype,p.adjust.method = "BH")
    pairwise.FGL.test <- pairwise.wilcox.test(outTab$FGL,outTab$Subtype,p.adjust.method = "BH")
  }

  if(n.moic > 2 & test.method == "parametric") {
    statistic <- "anova"
    FGA.test  <- summary(aov(outTab$FGA ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
    FGG.test  <- summary(aov(outTab$FGG ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
    FGL.test  <- summary(aov(outTab$FGL ~ outTab$Subtype))[[1]][["Pr(>F)"]][1]
    pairwise.FGA.test <- pairwise.t.test(outTab$FGA,outTab$Subtype,p.adjust.method = "BH")
    pairwise.FGG.test <- pairwise.t.test(outTab$FGG,outTab$Subtype,p.adjust.method = "BH")
    pairwise.FGL.test <- pairwise.t.test(outTab$FGL,outTab$Subtype,p.adjust.method = "BH")
  }

  # generate barplot
  FGA.col <- barcolor[1]
  FGG.col <- barcolor[2]
  FGL.col <- barcolor[3]

  p1 <- ggplot(summaryFGA, aes(x = Subtype, y = mean,fill=rep("0",nrow(summaryFGA)))) +
    geom_bar(stat = 'identity') +
    geom_errorbar(aes(ymax = mean+se, ymin = mean-se),position = position_dodge(0.9), width = 0.15) +
    annotate(geom="text",
             #hjust = ifelse(n.moic %% 2 == 0, 0.5, 0),
             x = n.moic / 2 + 0.5,
             #y = as.numeric(summaryFGA[which.max(summaryFGA$mean),"mean"] + summaryFGA[which.max(summaryFGA$mean),"se"]),
             y = as.numeric(summaryFGA[which.max(summaryFGA$mean),"mean"]),
             size = 8, angle = 90, fontface = "bold",
             #label = paste0(statistic, " p = ", formatC(FGA.test, digits = 1, format = "e")),
             label = cut(FGA.test,c(0,0.001,0.01,0.05,0.1,1),labels = c("****","***","**","*","."))) +
    scale_x_discrete(name = "",position = "top") +
    theme_bw() +
    theme(axis.line.y = element_line(size = 0.8),
          axis.ticks.y = element_line(size = 0.2),
          axis.text.y = element_blank(),
          axis.title.x = element_text(vjust = -0.3,size = 12),
          axis.text.x = element_text(size = 10, color = "black"),
          plot.margin = unit(c(0.3, -1.7, 0.3, 0.3), "lines"),
          legend.title = element_blank()) +
    coord_flip() +
    scale_fill_manual(values  = FGA.col, breaks = c("0"), labels = c("Copy number-altered genome")) +
    scale_y_reverse(expand = c(0.01,0),
                    name = "FGA (Fraction of Genome Altered)", position = "left")

  p2 <- ggplot(summaryFGGL, aes(x = Subtype, y = ifelse(class == 'FGG',mean,-mean),fill=class)) +
    geom_bar(stat = 'identity') +
    geom_errorbar(data=summaryFGGL[summaryFGGL$class=='FGG',],aes(ymax = mean+se, ymin =mean-se),position = position_dodge(0.9), width = 0.15) +
    geom_errorbar(data=summaryFGGL[summaryFGGL$class=='FGL',],aes(ymax = -mean-se, ymin =-mean+se),position = position_dodge(0.9), width = 0.15) +
    annotate(geom="text",

             x = n.moic / 2 + 0.5,
             y = -as.numeric(summaryFGL[which.max(summaryFGL$mean),"mean"]),
             size = 8, angle = 90, fontface = "bold",
             label = cut(FGL.test,c(0,0.001,0.01,0.05,0.1,1),labels = c("****","***","**","*","."))) +
    annotate(geom="text",
             x = n.moic / 2 + 0.5,
             y = as.numeric(summaryFGG[which.max(summaryFGG$mean),"mean"]),
             size = 8, angle = 90, fontface = "bold",
             label = cut(FGG.test,c(0,0.001,0.01,0.05,0.1,1),labels = c("****","***","**","*","."))) +
    scale_x_discrete(name = "") +
    theme_bw() +
    theme(axis.line.y = element_line(size = 0.8),
          axis.ticks.y = element_line(size = 0.2),
          axis.text.y = element_blank(),
          axis.title.x = element_text(vjust = -0.3, size = 12),
          axis.text.x = element_text(size = 10, color = "black"),
          plot.margin = unit(c(0.3, 0.3, 0.3, -1), "lines"),
          legend.title = element_blank()) +
    coord_flip() +
    scale_fill_manual(values  = c(FGL.col, FGG.col), breaks = c("FGL","FGG"),
                      labels = c("Copy number-lost genome","Copy number-gained genome")) +
    scale_y_continuous(expand = c(0.01,0),
                       name = "FGL or FGG (Fraction of Genome Lost or Gained)")

  pp <- ggplot() +
    # geom_text(data = summaryFGGL,
    #           aes(label = Subtype, x=Subtype), y = 0.5,
    #           size = 0.8*11/.pt, # match font size to theme
    #           hjust = 0.5, vjust = 0.5) +
    geom_label(data = summaryFGGL,
               aes(label = Subtype, x = Subtype, fill = Subtype),
               y = 0.5,
               color = "white",
               size = 0.9*11/.pt, # match font size to theme
               hjust = 0.4, vjust = 0.5) +
    scale_fill_manual(values = clust.col) +
    theme_minimal()+
    theme(axis.line.y =element_blank(),
          axis.ticks.y =element_blank(),
          axis.text.y =element_blank(),
          axis.title.y =element_blank(),
          axis.title.x =element_blank(),
          plot.margin = unit(c(0.3, 0, 0.3, 0), "lines")
    ) +
    guides(fill = FALSE) +
    coord_flip() +
    scale_y_reverse()

  pal <- p1 + pp + p2 +
    plot_layout(widths = c(7,1,7), guides = 'collect') & theme(legend.position = 'top')

  # save to pdf
  if(is.null(fig.name)) {
    outFig <- "barplot of FGA.pdf"
  } else {
    outFig <- paste0(fig.name,".pdf")
  }
  ggsave(file.path(fig.path, outFig), width = width, height = height)

  # print to screen
  print(pal)

  if(n.moic > 2) {
    return(list(summary = outTab,
                FGA.p.value = FGA.test,
                pairwise.FGA.test = pairwise.FGA.test,
                FGG.p.value = FGG.test,
                pairwise.FGG.test = pairwise.FGG.test,
                FGL.p.value = FGL.test,
                pairwise.FGL.test = pairwise.FGL.test,
                test.method = statistic))
  } else {
    return(list(summary = outTab,
                FGA.p.value = FGA.test,
                FGG.p.value = FGG.test,
                FGL.p.value = FGL.test,
                test.method = statistic))
  }
}
