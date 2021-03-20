#' @name compDrugsen
#' @title Comparison of drug sensitivity
#' @description This function estimates the IC50 of specific drugs for each Subtype by developing a ridge regression predictive model based on all/specific cell lines derived from Genomics of Drug Sensitivity in Cancer (GDSC, \url{https://www.cancerrxgene.org/}).
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param norm.expr A matrix of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.
#' @param drugs A string vector to indicate the names of the drugs for which you would like to predict sensitivity, one of Erlotinib, Rapamycin, Sunitinib, PHA-665752, MG-132, Paclitaxel, Cyclopamine, AZ628, Sorafenib, VX-680, Imatinib, TAE684, Crizotinib, Saracatinib, S-Trityl-L-cysteine, Z-LLNle-CHO, Dasatinib, GNF-2, CGP-60474, CGP-082996, A-770041, WH-4-023, WZ-1-84, BI-2536, BMS-536924, BMS-509744, CMK, Pyrimethamine, JW-7-52-1, A-443654, GW843682X, MS-275, Parthenolide, KIN001-135, TGX221, Bortezomib, XMD8-85, Roscovitine, Salubrinal, Lapatinib, GSK269962A, Doxorubicin, Etoposide, Gemcitabine, Mitomycin C, Vinorelbine, NSC-87877, Bicalutamide, QS11, CP466722, Midostaurin, CHIR-99021, AP-24534, AZD6482, JNK-9L, PF-562271, HG-6-64-1, JQ1, JQ12, DMOG, FTI-277, OSU-03012, Shikonin, AKT inhibitor VIII, Embelin, FH535, PAC-1, IPA-3, GSK-650394, BAY 61-3606, 5-Fluorouracil, Thapsigargin, Obatoclax Mesylate, BMS-754807, Lisitinib, Bexarotene, Bleomycin, LFM-A13, GW-2580, AUY922, Phenformin, Bryostatin 1, Pazopanib, LAQ824, Epothilone B, GSK1904529A, BMS345541, Tipifarnib, BMS-708163, Ruxolitinib, AS601245, Ispinesib Mesylate, TL-2-105, AT-7519, TAK-715, BX-912, ZSTK474, AS605240, Genentech Cpd 10, GSK1070916, KIN001-102, LY317615, GSK429286A, FMK, QL-XII-47, CAL-101, UNC0638, XL-184, WZ3105, XMD14-99, AC220, CP724714, JW-7-24-1, NPK76-II-72-1, STF-62247, NG-25, TL-1-85, VX-11e, FR-180204, Tubastatin A, Zibotentan, YM155, NSC-207895, VNLG/124, AR-42, CUDC-101, Belinostat, I-BET-762, CAY10603, Linifanib , BIX02189, CH5424802, EKB-569, GSK2126458, KIN001-236, KIN001-244, KIN001-055, KIN001-260, KIN001-266, Masitinib, MP470, MPS-1-IN-1, BHG712, OSI-930, OSI-027, CX-5461, PHA-793887, PI-103, PIK-93, SB52334, TPCA-1, TG101348, Foretinib, Y-39983, YM201636, Tivozanib, GSK690693, SNX-2112, QL-XI-92, XMD13-2, QL-X-138, XMD15-27; two common chemodrugs (i.e., Cisplatin and Paclitaxel) will be analyzed by default if no indication.
#' @param tissueType A string value to specify if you would like to train the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated); Allowed values contain c("all", "aero_digestive_tract", "blood", "bone", "breast", "digestive_system", "lung", "nervous_system", "skin", "urogenital_system") and "all" by default.
#' @param clust.col A string vector storing colors for annotating each Subtype.
#' @param prefix  A string value to indicate the prefix of the output plot.
#' @param fig.path A string value to indicate the output path for storing the boxviolin plot.
#' @param seed A integer value to indicate the seed for reproducing ridge regression.
#' @param width A numeric value to indicate the width of boxviolin plot.
#' @param height A numeric value to indicate the height of boxviolin plot.
#' @param test.method A string value to indicate the method for statistical testing. Allowed values contain c('nonparametric', 'parametric'); nonparametric means two-sample wilcoxon rank sum test for two subtypes and Kruskal-Wallis rank sum test for multiple subtypes; parametric means two-sample t-test when only two subtypes are identified, and anova for multiple subtypes comparison; "nonparametric" by default.
#' @return Data.frame(s) storing the estimated IC50 of specified drugs per sample within each Subtype.
#' @export
#' @import ggplot2
#' @importFrom ggpubr stat_compare_means
#' @references Geeleher P, Cox N, Huang R S. (2014). pRRophetic: an R package for prediction of clinical chemotherapeutic response from tumor gene expression levels. PLoS One, 9(9):e107468.
#'
#' Geeleher P, Cox N J, Huang R S. (2014). Clinical drug response can be predicted using baseline gene expression levels and in vitro drug sensitivity in cell lines. Genome Biol, 15(3):1-12.
#' @examples # There is no example and please refer to vignette.
compDrugsen <- function(moic.res    = NULL,
                        norm.expr   = NULL,
                        drugs       = c("Cisplatin", "Paclitaxel"),
                        tissueType  = "all",
                        test.method = "nonparametric",
                        clust.col   = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                        prefix      = NULL,
                        seed        = 123456,
                        fig.path    = getwd(),
                        width       = 5,
                        height      = 5) {

  if(!is.element(test.method, c("nonparametric","parametric"))) {
    stop("test.method can be one of nonparametric or parametric.")
  }

  # data processing
  comsam <- intersect(moic.res$clust.res$samID, colnames(norm.expr))
  # check data
  if(length(comsam) == nrow(moic.res$clust.res)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(moic.res$clust.res)-length(comsam))," samples mismatched from current subtypes."))
  }
  moic.res$clust.res <- moic.res$clust.res[comsam,,drop = FALSE]
  norm.expr <- norm.expr[,comsam]

  n.moic <- length(unique(moic.res$clust.res$clust))
  sam.order <- moic.res$clust.res[order(moic.res$clust.res$clust, decreasing = FALSE), "samID"]
  colvec <- clust.col[1:length(unique(moic.res$clust.res$clust))]
  names(colvec) <- paste0("CS",unique(moic.res$clust.res$clust))

  annCol <- data.frame("Subtype" = paste0("CS",moic.res$clust.res[sam.order,"clust"]),
                       samID = sam.order,
                       row.names = sam.order,
                       stringsAsFactors = FALSE)

  if(max(norm.expr) < 25 | (max(norm.expr) >= 25 & min(norm.expr) < 0)) {
    message("--expression profile seems to have veen standardised (z-score or log transformation), no more action will be performed.")
    gset <- norm.expr
  }
  if(max(norm.expr) >= 25 & min(norm.expr) >= 0){
    message("--log2 transformation done for expression data.")
    gset <- log2(norm.expr + 1)
  }

  # drug sensitivity prediction
  predictedPtype <- predictedBoxdat <- list()

  for (drug in drugs) {
    set.seed(seed)

    predictedPtype[[drug]] <- quiet(pRRopheticPredict(testMatrix    = as.matrix(gset[,rownames(annCol)]),
                                                      drug          = drug,
                                                      tissueType    = tissueType,
                                                      dataset       = "cgp2016",
                                                      minNumSamples = 5,
                                                      selection     = 1)) # 1 indicate if multiple genes existed, mean value will be considered

    if(!all(names(predictedPtype[[drug]]) == rownames(annCol))) {stop("name mismatched!\n")}

    predictedBoxdat[[drug]] <- data.frame("Est.IC50"        = predictedPtype[[drug]],
                                          "Subtype"         = as.character(annCol$Subtype),
                                           row.names        = names(predictedPtype[[drug]]),
                                           stringsAsFactors = FALSE)
    message(drug," done...")

    # generate boxviolin plot with statistical testing
    if(n.moic == 2 & test.method == "nonparametric") {
      statistic = "wilcox.test"
      ic50.test  <- wilcox.test(predictedBoxdat[[drug]]$Est.IC50 ~ predictedBoxdat[[drug]]$Subtype)$p.value
      cat(paste0("Wilcoxon rank sum test p value = ", formatC(ic50.test, format = "e", digits = 2), " for ", drug))
    }
    if(n.moic == 2 & test.method == "parametric") {
      statistic = "t.test"
      ic50.test  <- t.test(predictedBoxdat[[drug]]$Est.IC50 ~ predictedBoxdat[[drug]]$Subtype)$p.value
      cat(paste0("Student's t test p value = ", formatC(ic50.test, format = "e", digits = 2), " for ", drug))
    }
    if(n.moic > 2 & test.method == "nonparametric") {
      statistic = "kruskal.test"
      ic50.test  <- kruskal.test(predictedBoxdat[[drug]]$Est.IC50 ~ predictedBoxdat[[drug]]$Subtype)$p.value
      pairwise.ic50.test <- pairwise.wilcox.test(predictedBoxdat[[drug]]$Est.IC50,predictedBoxdat[[drug]]$Subtype,p.adjust.method = "BH")
      cat(paste0(drug,": Kruskal-Wallis rank sum test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:\n"))
      print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))
    }
    if(n.moic > 2 & test.method == "parametric") {
      statistic = "anova"
      ic50.test  <- summary(aov(predictedBoxdat[[drug]]$Est.IC50 ~ predictedBoxdat[[drug]]$Subtype))[[1]][["Pr(>F)"]][1]
      pairwise.ic50.test <- pairwise.t.test(predictedBoxdat[[drug]]$Est.IC50,predictedBoxdat[[drug]]$Subtype,p.adjust.method = "BH")
      cat(paste0(drug,": One-way anova test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise Student's t test with Benjamini-Hochberg adjustment presents below:\n"))
      print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))
    }

    p <- ggplot(data = predictedBoxdat[[drug]],
                aes(x = Subtype, y = Est.IC50, fill = Subtype)) +
      scale_fill_manual(values = colvec) +
      geom_violin(alpha = 0.4, position = position_dodge(width = .75),
                  size = 0.8, color = "black") +
      geom_boxplot(notch = TRUE, outlier.size = -1,
                   color = "black", lwd = 0.8, alpha = 0.7) +
      geom_point(shape = 21, size = 2,
                 position = position_jitterdodge(),
                 color = "black", alpha = 1) +
      theme_classic() +
      ylab(bquote("Estimated IC"[50]~"of"~.(drug))) + xlab("") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.ticks = element_line(size = 0.2, color = "black"),
            axis.ticks.length = unit(0.2, "cm"),
            legend.position = "none",
            axis.title = element_text(size = 15),
            axis.text = element_text(size = 10)) +
      # add statistical inference
      stat_compare_means(method = statistic,
                         hjust = ifelse(n.moic %% 2 == 0, 0.5, 0),
                         label.x = ifelse(n.moic %% 2 == 0, n.moic / 2 + 0.5, n.moic / 2),
                         label.y = min(predictedBoxdat[[drug]]$Est.IC50))

    # save to pdf
    if(is.null(prefix)) {
      outFig <- paste0("boxviolin of estimated ic50 for ",drug,".pdf")
    } else {
      outFig <- paste0(prefix, " for ", drug,".pdf")
    }
    ggsave(file.path(fig.path, outFig), width = width, height = height)
    # print to screen
    print(p)
  }
  return(predictedBoxdat)
}
