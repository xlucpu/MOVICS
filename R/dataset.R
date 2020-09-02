#' GDSC 2016 gene expression data
#'
#' A matrix frame storing the GDSC 2016 gene expression data.
#'
#' @format A matrix with 17419 rows and 1018 variables.
#' \describe{contains "cgp2016ExprRma" the 2016 gene expression data. Data was obtained from \url{http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/preprocessed/Cell_line_RMA_proc_basalExp.txt.zip}. Cosmic Ids (in the column names) were mapped to cell line names using data from this file: \url{ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/Cell_Lines_Details.xlsx}}
#' @source \url{https://osf.io/5xvsg/}
"cgp2016ExprRma"

#' GDSC 2016 drug data
#'
#' A data frame storing the GDSC 2016 drug data.
#'
#' @format A data frame with 224510 rows and 13 variables.
#' \describe{contains "drugData2016" the 2016 drug IC50 data, downloaded from \url{http://www.cancerrxgene.org/translation/drug/download#ic50)}}
#' @source \url{https://osf.io/5xvsg/}
"drugData2016"

