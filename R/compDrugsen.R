#' @name compDrugsen
#' @title Comparison of drug sensitivity
#' @description This function estimates the IC50 of specific drugs for each Subtype by developing a ridge regression predictive model based on all/specific cell lines derived from Genomics of Drug Sensitivity in Cancer (GDSC, \url{https://www.cancerrxgene.org/}).
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param norm.expr A matrix of normalized expression data with rows for genes and columns for samples; FPKM or TPM without log2 transformation is recommended.
#' @param drugs A string vector to indicate the names of the drugs for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439; two common chemodrugs (i.e., Cisplatin and Paclitaxel) will be analyzed by default if no indication.
#' @param tissueType A string value to specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated); Allowed values contain c("all", "aero_digestive_tract", "blood", "bone", "breast", "digestive_system", "lung", "nervous_system", "skin", "urogenital_system") and "all" by default.
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
      cat(paste0(drug,": Kruskal-Wallis rank sum test p value = ", formatC(TMB.test, format = "e", digits = 2),"\npost-hoc pairwise wilcoxon rank sum test with Benjamini-Hochberg adjustment presents below:"))
      print(formatC(pairwise.ic50.test$p.value, format = "e", digits = 2))
    }
    if(n.moic > 2 & test.method == "parametric") {
      statistic = "anova"
      ic50.test  <- summary(aov(predictedBoxdat[[drug]]$Est.IC50 ~ predictedBoxdat[[drug]]$Subtype))[[1]][["Pr(>F)"]][1]
      pairwise.ic50.test <- pairwise.t.test(predictedBoxdat[[drug]]$Est.IC50,predictedBoxdat[[drug]]$Subtype,p.adjust.method = "BH")
      cat(paste0(drug,": One-way anova test p value = ", formatC(ic50.test, format = "e", digits = 2),"\npost-hoc pairwise Student's t test with Benjamini-Hochberg adjustment presents below:"))
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
