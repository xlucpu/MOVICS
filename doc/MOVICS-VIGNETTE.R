## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE--------------------------------------------------------------
Sys.setenv(LANGUAGE = "en")

## ---- warning=FALSE, message=FALSE, echo=TRUE, eval=FALSE---------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  if (!require("devtools"))
#      install.packages("devtools")
#  devtools::install_github("xlucpu/MOVICS", host = "https://api.github.com")

## ---- echo=FALSE, eval=TRUE, warning=FALSE, message=FALSE---------------------
library("dplyr")
library("knitr")
library("kableExtra")
#library("devtools")
#load_all()

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
library("MOVICS")

## ---- eval=TRUE---------------------------------------------------------------
# load example data of breast cancer
load(system.file("extdata", "brca.tcga.RData", package = "MOVICS", mustWork = TRUE))
load(system.file("extdata", "brca.yau.RData",  package = "MOVICS", mustWork = TRUE))

## ---- eval=TRUE---------------------------------------------------------------
# print name of example data
names(brca.tcga)
names(brca.yau)

# extract multi-omics data
mo.data   <- brca.tcga[1:4]

# extract raw count data for downstream analyses
count     <- brca.tcga$count

# extract fpkm data for downstream analyses
fpkm      <- brca.tcga$fpkm

# extract maf for downstream analysis
maf       <- brca.tcga$maf

# extract segmented copy number for downstream analyses
segment   <- brca.tcga$segment

# extract survival information
surv.info <- brca.tcga$clin.info

## ---- eval=TRUE---------------------------------------------------------------
# scenario 1: considering we are dealing with an expression data that have 2 rows with NA values
tmp       <- brca.tcga$mRNA.expr # get expression data
dim(tmp) # check data dimension
tmp[1,1]  <- tmp[2,2] <- NA # set 2 rows with NA values
tmp[1:3,1:3] # check data
elite.tmp <- getElites(dat       = tmp,
                       method    = "mad",
                       na.action = "rm", # NA values will be removed
                       elite.pct = 1) # elite.pct equals to 1 means all (100%) features after NA removal will be selected even using mad method
dim(elite.tmp$elite.dat) # check dimension again and see that we have removed 2 rows with NA data

elite.tmp <- getElites(dat       = tmp,
                       method    = "mad",
                       na.action = "impute", # NA values will be imputed
                       elite.pct = 1) 
dim(elite.tmp$elite.dat) # all data kept
elite.tmp$elite.dat[1:3,1:3] # NA values have been imputed 

# scenario 2: considering we are dealing with continuous data and use mad or sd to select elites
tmp       <- brca.tcga$mRNA.expr # get expression data with 500 features
elite.tmp <- getElites(dat       = tmp,
                       method    = "mad",
                       elite.pct = 0.1) # this time only top 10% features with high mad values are kept
dim(elite.tmp$elite.dat) # get 50 elite left

elite.tmp <- getElites(dat       = tmp,
                       method    = "sd",
                       elite.num = 100, # this time only top 100 features with high sd values are kept
                       elite.pct = 0.1) # this time elite.pct argument will be disabled because elite.num has been already indicated.
dim(elite.tmp$elite.dat) # get 100 elites left

# scenario 3: considering we are dealing with data and use cox to select elite
tmp       <- brca.tcga$mRNA.expr # get expression data 
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                       p.cutoff  = 0.05,
                       elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites
dim(elite.tmp$elite.dat) # get 125 elites
table(elite.tmp$unicox$pvalue < 0.05) # 125 genes have nominal pvalue < 0.05 in univariate Cox regression

tmp       <- brca.tcga$mut.status # get mutation data 
elite.tmp <- getElites(dat       = tmp,
                       method    = "cox",
                       surv.info = surv.info, # must provide survival information with 'futime' and 'fustat'
                       p.cutoff  = 0.05,
                       elite.num = 100) # this time elite.num argument will be disabled because cox method refers to p.cutoff to select elites
dim(elite.tmp$elite.dat) # get 3 elites
table(elite.tmp$unicox$pvalue < 0.05) # 3 mutations have nominal pvalue < 0.05

# scenario 4: considering we are dealing with mutation data using freq to select elites
tmp       <- brca.tcga$mut.status # get mutation data 
rowSums(tmp) # check mutation frequency
elite.tmp <- getElites(dat       = tmp,
                       method    = "freq", # must set as 'freq'
                       elite.num = 80, # note: in this scenario elite.num refer to frequency of mutation
                       elite.pct = 0.1) # discard because elite.num has been already indicated
rowSums(elite.tmp$elite.dat) # only genes that are mutated in over than 80 samples are kept as elites

elite.tmp <- getElites(dat       = tmp,
                       method    = "freq", # must set as 'freq'
                       elite.pct = 0.2) # note: in this scenario elite.pct refer to frequency of mutation / sample size
rowSums(elite.tmp$elite.dat) # only genes that are mutated in over than 0.2*643=128.6 samples are kept as elites

# get mo.data list just like below (not run)
# mo.data <- list(omics1 = elite.tmp$elite.dat,
#                 omics2 = ...)

## ---- fig.align="center", fig.width=6, fig.height=6, fig.cap="Figure 1. Identification of optimal cluster number by calculating CPI (blue line) and Gaps-statistics (red line) in TCGA-BRCA cohort.", eval=TRUE----
# identify optimal clustering number (may take a while)
optk.brca <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-BRCA")

## ---- eval=TRUE---------------------------------------------------------------
# perform iClusterBayes (may take a while)
iClusterBayes.res <- getiClusterBayes(data        = mo.data,
                                      N.clust     = 5,
                                      type        = c("gaussian","gaussian","gaussian","binomial"),
                                      n.burnin    = 1800,
                                      n.draw      = 1200,
                                      prior.gamma = c(0.5, 0.5, 0.5, 0.5),
                                      sdev        = 0.05,
                                      thin        = 3)

## ---- eval=FALSE--------------------------------------------------------------
#  iClusterBayes.res <- getMOIC(data        = mo.data,
#                               N.clust     = 5,
#                               methodslist = "iClusterBayes", # specify only ONE algorithm here
#                               type        = c("gaussian","gaussian","gaussian","binomial"), # data type corresponding to the list
#                               n.burnin    = 1800,
#                               n.draw      = 1200,
#                               prior.gamma = c(0.5, 0.5, 0.5, 0.5),
#                               sdev        = 0.05,
#                               thin        = 3)

## ---- fig.show='hide', message=TRUE, eval=TRUE--------------------------------
# perform multi-omics integrative clustering with the rest of 9 algorithms
moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 5,
                         type        = c("gaussian", "gaussian", "gaussian", "binomial"))

# attach iClusterBayes.res as a list using append() to moic.res.list with 9 results already
moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))

# save moic.res.list to local path
save(moic.res.list, file = "moic.res.list.rda")

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  load(file = "iClusterBayes.res.rda")
#  load(file = "moic.res.list.rda")

## ---- fig.align="center", fig.width=7, fig.height=6, fig.cap="Figure 2. Consensus heatmap based on results from 10 multi-omics integrative clustering algorithms with cluster number of 5.", eval=TRUE----
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")

## ---- eval=TRUE---------------------------------------------------------------
# convert beta value to M value for stronger signal
indata <- mo.data
indata$meth.beta <- log2(indata$meth.beta / (1 - indata$meth.beta))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,F)) # no scale for mutation

## ---- eval=TRUE---------------------------------------------------------------
feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA.expr"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA.expr"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "meth.beta"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "mut.status"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4)

## ---- fig.align="center", fig.width=9, fig.height=8.5, fig.cap="Figure 3. Comprehensive heatmap of multi-omics integrative clustering by iClusterBayes with annotation of potential drivers.", eval=TRUE----
# set color for each omics data
# if no color list specified all subheatmaps will be unified to green and red color pattern
mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
meth.col   <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, meth.col, mut.col)

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = iClusterBayes.res$clust.res, # cluster results
             clust.dend    = NULL, # no dendrogram
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = NULL, # no annotation for samples
             annColors     = NULL, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")

## ---- fig.align="center", fig.width=9, fig.height=9, fig.cap="Figure 4. Comprehensive heatmap of multi-omics integrative clustering by COCA with dendrogram for samples.", eval=TRUE----
# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = moic.res.list$COCA$clust.res, # cluster results
             clust.dend    = moic.res.list$COCA$clust.dend, # show dendrogram for samples
             color         = col.list,
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF COCA")

## ---- fig.align="center", fig.width=9, fig.height=9.5, fig.cap="Figure 5. Comprehensive heatmap based on consensus across 10 algorithms with clinicopathological annotation.", eval=TRUE----
# extract PAM50, pathologic stage and age for sample annotation
annCol    <- surv.info[,c("PAM50", "pstage", "age"), drop = FALSE]

# generate corresponding colors for sample annotation
annColors <- list(age    = circlize::colorRamp2(breaks = c(min(annCol$age),
                                                           median(annCol$age),
                                                           max(annCol$age)), 
                                                colors = c("#0000AA", "#555555", "#AAAA00")),
                  PAM50  = c("Basal" = "blue",
                            "Her2"   = "red",
                            "LumA"   = "yellow",
                            "LumB"   = "green",
                            "Normal" = "black"),
                  pstage = c("T1"    = "green",
                             "T2"    = "blue",
                             "T3"    = "red",
                             "T4"    = "yellow", 
                             "TX"    = "black"))

# comprehensive heatmap (may take a while)
getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.FPKM","lncRNA.FPKM","M value","Mutated"),
             clust.res     = cmoic.brca$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F), # show no dendrogram for features
             annRow        = NULL, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")

## ---- fig.align="center", fig.width=6, fig.height=7, fig.cap="Figure 6. Kaplan-Meier survival curve of 5 identified subtypes of breast cancer in TCGA-BRCA cohort.", eval=TRUE----
# survival comparison
surv.brca <- compSurv(moic.res         = cmoic.brca,
                      surv.info        = surv.info,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h",
                      fig.name         = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC") # draw horizontal line at median survival

print(surv.brca)

## ---- eval=TRUE---------------------------------------------------------------
clin.brca <- compClinvar(moic.res      = cmoic.brca,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("PAM50","pstage","fustat"), # features that are considered categorical variables
                         nonnormalVars = "futime", # feature(s) that are considered using nonparametric test
                         exactVars     = "pstage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")


## ---- eval=FALSE--------------------------------------------------------------
#  print(clin.brca$compTab)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
clin.brca$compTab %>%
  kbl(caption = "Table 1. Comparison of clinical features among 5 identified subtype of breast cancer in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.align="center", fig.width=8, fig.height=3, fig.cap="Figure 7. Mutational OncoPrint of 5 identified subtypes of breast cancer in TCGA-BRCA cohort.", eval=TRUE----
# mutational frequency comparison
mut.brca <- compMut(moic.res     = cmoic.brca,
                    mut.matrix   = brca.tcga$mut.status, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 6, 
                    height       = 2,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")

## ---- eval=FALSE--------------------------------------------------------------
#  print(mut.brca)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
mut.brca %>%
  kbl(caption = "Table 2. Comparison of mutational frequency among 5 identified subtype of breast cancer in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- eval=TRUE---------------------------------------------------------------
head(maf)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  head(maf) %>%
#    kbl(caption = "Table . Demo of MAF data with eligible column names.") %>%
#    kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.align="center", fig.width=6, fig.height=6, fig.cap="Figure 8. Comparison of TMB and TiTv among 5 identified subtypes of breast cancer in TCGA-BRCA cohort.", eval=TRUE----
# compare TMB
tmb.brca <- compTMB(moic.res     = cmoic.brca,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

## ---- eval=FALSE--------------------------------------------------------------
#  head(tmb.brca$TMB.dat)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
head(tmb.brca$TMB.dat) %>%
  kbl(caption = "Table 3. Demo of comparison of TMB among 5 identified subtype of breast cancer in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## -----------------------------------------------------------------------------
# change column names of segment data
colnames(segment) <- c("sample","chrom","start","end","value")

## ---- eval=TRUE---------------------------------------------------------------
head(segment)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  head(segment) %>%
#    kbl(caption = "Table . Demo of segmented copy number data with eligible column names.") %>%
#    kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.align="center", fig.width=10, fig.height=2.5, fig.cap="Figure 9. Barplot of fraction genome altered among 5 identified subtypes of breast cancer in TCGA-BRCA cohort.", eval=TRUE----
# compare FGA, FGG, and FGL
fga.brca <- compFGA(moic.res     = cmoic.brca,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA")

## ---- eval=FALSE--------------------------------------------------------------
#  head(fga.brca$summary)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
head(fga.brca$summary) %>%
  kbl(caption = "Table 4. Demo of comparison of fraction genome altered among 5 identified subtype of breast cancer in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.show = "hold", out.width = "50%", fig.align = "default", fig.width=8, fig.height=6, fig.cap="Figure 10. Boxviolins for estimated IC50 of Cisplatin and Paclitaxel among 5 identified subtypes of breast cancer in TCGA-BRCA cohort.", eval=TRUE----
# drug sensitivity comparison
drug.brca <- compDrugsen(moic.res    = cmoic.brca,
                         norm.expr   = fpkm[,cmoic.brca$clust.res$samID], # double guarantee sample order
                         drugs       = c("Cisplatin", "Paclitaxel"), # a vector of names of drug in GDSC
                         tissueType  = "breast", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50") 

## ---- eval=FALSE--------------------------------------------------------------
#  head(drug.brca$Cisplatin)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
head(drug.brca$Cisplatin) %>%
  kbl(caption = "Table 5. Demo of estimated IC50 for Cisplatin among 5 identified subtype of breast cancer in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.align="center", fig.width=8, fig.height=5, fig.cap="Figure 11. Agreement of 5 identified subtypes of breast cancer with PAM50 classification and pathological stage in TCGA-BRCA cohort.", eval=TRUE----
# customize the factor level for pstage
surv.info$pstage <- factor(surv.info$pstage, levels = c("TX","T1","T2","T3","T4"))

# agreement comparison (support up to 6 classifications include current subtype)
agree.brca <- compAgree(moic.res  = cmoic.brca,
                        subt2comp = surv.info[,c("PAM50","pstage")],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE")

## ---- eval=FALSE--------------------------------------------------------------
#  print(agree.brca)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
agree.brca %>%
  kbl(caption = "Table 6. Agreement of 5 identified subtypes with PAM50 classification and pathological stage in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- eval=TRUE---------------------------------------------------------------
# run DEA with edgeR
runDEA(dea.method = "edger",
       expr       = count, # raw count data
       moic.res   = cmoic.brca,
       prefix     = "TCGA-BRCA") # prefix of figure name

# run DEA with DESeq2
runDEA(dea.method = "deseq2",
       expr       = count,
       moic.res   = cmoic.brca,
       prefix     = "TCGA-BRCA")

# run DEA with limma
runDEA(dea.method = "limma",
       expr       = fpkm, # normalized expression data
       moic.res   = cmoic.brca,
       prefix     = "TCGA-BRCA")

## ---- fig.align="center", fig.width=7, fig.height=6, fig.cap="Figure 12. Heatmap of subtype-specific upregulated biomarkers using edgeR for 5 identified subtypes in TCGA-BRCA cohort.", eval=TRUE----
# choose edgeR result to identify subtype-specific up-regulated biomarkers
marker.up <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "edger", # name of DEA method
                       prefix        = "TCGA-BRCA", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = fpkm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")

## ---- eval=FALSE--------------------------------------------------------------
#  # check the upregulated biomarkers
#  head(marker.up$templates)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
head(marker.up$templates) %>%
  kbl(caption = "Table 7. Demo of subtype-specific upregulated biomarkers for 5 identified subtypes of breast cancer in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.align="center", fig.width=7, fig.height=6, fig.cap="Figure 13. Heatmap of subtype-specific downregulated biomarkers using limma for 5 identified subtypes in TCGA-BRCA cohort.", eval=TRUE----
# choose limma result to identify subtype-specific down-regulated biomarkers
marker.dn <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma",
                       prefix        = "TCGA-BRCA",
                       dirct         = "down",
                       n.marker      = 50, # switch to 50
                       doplot        = TRUE,
                       norm.expr     = fpkm,
                       annCol        = annCol,
                       annColors     = annColors,
                       fig.name      = "DOWNREGULATED BIOMARKER HEATMAP")

## ---- eval=TRUE---------------------------------------------------------------
# MUST locate ABSOLUTE path of msigdb file
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

## ---- fig.align="center", fig.width=10, fig.height=8, fig.cap="Figure 14. Heatmap of subtype-specific upregulated pathways using edgeR algorithm for 5 identified subtypes in TCGA-BRCA cohort.", eval=TRUE----
# run GSEA to identify up-regulated GO pathways using results from edgeR
gsea.up <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "edger", # name of DEA method
                   prefix       = "TCGA-BRCA", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.25, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")

## ---- eval=FALSE--------------------------------------------------------------
#  print(gsea.up$gsea.list$CS1[1:6,3:6])

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
gsea.up$gsea.list$CS1[1:6,3:6] %>%
  kbl(caption = "Table 8. Demo of GSEA results for the first cancer subtype (CS1) of breast cancer in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- eval=FALSE--------------------------------------------------------------
#  head(round(gsea.up$grouped.es,3))

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
head(round(gsea.up$grouped.es,3)) %>%
  kbl(caption = "Table 9. Demo of subtype-specific enrichment scores among 5 identified subtypes of breast cancer in TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.align="center", fig.width=10, fig.height=8, fig.cap="Figure 15. Heatmap of subtype-specific downregulated pathways using limma algorithm for 5 identified subtypes in TCGA-BRCA cohort.", eval=TRUE----
# run GSEA to identify down-regulated GO pathways using results from DESeq2
gsea.dn <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "deseq2",
                   prefix       = "TCGA-BRCA",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = fpkm,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "ssgsea", # switch to ssgsea
                   norm.method  = "median", # switch to median
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP") 

## ---- fig.align="center", fig.width=6, fig.height=6, fig.cap="Figure 16. Heatmap of NTP in Yau cohort using subtype-specific upregulated biomarkers identified from TCGA-BRCA cohort", eval=TRUE----
# run NTP in Yau cohort by using up-regulated biomarkers
brca.pred <- runNTP(expr      = brca.yau$mRNA.expr,
                    templates = marker.up$templates, # the template has been already prepared in runMarker()
                    scale     = TRUE, # scale input data
                    center    = TRUE, # center input data
                    doPlot    = TRUE, # to generate heatmap
                    fig.name  = "NTP HEATMAP FOR YAU") 

## ---- eval=FALSE--------------------------------------------------------------
#  head(brca.pred$ntp.res)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
head(brca.pred$ntp.res) %>%
  kbl(caption = "Table 10. Demo of predicted subtypes in Yau cohort by NTP using subtype-specific upregulated biomarkers identified from TCGA-BRCA cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.align="center", fig.width=6, fig.height=7, fig.cap="Figure 17. Kaplan-Meier survival curve of predicted 5 subtypes of breast cancer in Yau cohort.", eval=TRUE----
# compare survival outcome in Yau cohort
surv.yau <- compSurv(moic.res         = brca.pred,
                     surv.info        = brca.yau$clin.info,
                     convt.time       = "m", # switch to year
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

print(surv.yau)

## ---- fig.align="center", fig.width=8, fig.height=5, fig.cap="Figure 18. Agreement of predicted 5 subtypes of breast cancer with PAM50 classification in Yau cohort.", eval=TRUE----
# compare agreement in Yau cohort
agree.yau <- compAgree(moic.res  = brca.pred,
                       subt2comp = brca.yau$clin.info[, "PAM50", drop = FALSE],
                       doPlot    = TRUE,
                       fig.name  = "YAU PREDICTEDMOIC WITH PAM50")


## ---- eval=FALSE--------------------------------------------------------------
#  print(agree.yau)

## ---- echo=FALSE, eval=TRUE---------------------------------------------------
agree.yau %>%
  kbl(caption = "Table 11. Agreement of 5 predicted subtypes of breast cancer with PAM50 classification in Yau cohort.") %>%
  kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- eval=TRUE---------------------------------------------------------------
# include original clinical information as `clust.res` and a string value for `mo.method` to a list
pseudo.moic.res                 <- list("clust.res" = surv.info,
                                        "mo.method" = "PAM50")

# make pseudo samID
pseudo.moic.res$clust.res$samID <- rownames(pseudo.moic.res$clust.res)

# make pseudo clust using a mapping relationship
pseudo.moic.res$clust.res$clust <- sapply(pseudo.moic.res$clust.res$PAM50,
                                          switch,
                                          "Basal"   = 1, # relabel Basal as 1
                                          "Her2"    = 2, # relabel Her2 as 2
                                          "LumA"    = 3, # relabel LumA as 3
                                          "LumB"    = 4, # relabel LumnB as 4
                                          "Normal"  = 5) # relabel Normal as 5

## ---- eval=TRUE---------------------------------------------------------------
head(pseudo.moic.res$clust.res)

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  head(pseudo.moic.res$clust.res) %>%
#    kbl(caption = "Table . Demo of pseudo object for downstream analyses in MOVICS.") %>%
#    kable_classic(full_width = TRUE, html_font = "Calibri")

## ---- fig.align="center", fig.width=6, fig.height=7, fig.cap="Figure 19. Kaplan-Meier survival curve of PAM50 subtypes of breast cancer with pseudo input in TCGA-BRCA cohort.", eval=TRUE----
# survival comparison
pam50.brca <- compSurv(moic.res         = pseudo.moic.res,
                       surv.info        = surv.info,
                       convt.time       = "y", # convert day unit to year
                       surv.median.line = "h", # draw horizontal line at median survival
                       fig.name         = "KAPLAN-MEIER CURVE OF PAM50 BY PSEUDO")

## ---- echo=TRUE, eval=TRUE----------------------------------------------------
sessionInfo()

## ---- echo=FALSE, eval=FALSE--------------------------------------------------
#  save.image("MOVICS.RData")

