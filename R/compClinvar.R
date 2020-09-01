#' @name compClinvar
#' @title Comparison of clinical variables
#' @description Create a table summarizing all baseline variables (both continuous and categorical) stratifying by current identified Subtypes and performing statistical tests. The object gives a table that is easy to use in medical research papers.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param var2comp A data.frame of clinical variables that need to compare among current subtypes with rownames for samples and columns for variable names.
#' @param strata A string value to indicate the stratifying (subtype) variable; 'Subtype' by default.
#' @param factorVars A string vectors to indicate the categorical variables. If omitted, only factors are considered categorical variables.
#' @param nonnormalVars A string vector to specify the variables for which the p-values should be those of nonparametric tests. By default all p-values are from normal assumption-based tests (oneway.test)., Default: NULL
#' @param exactVars A string vector to specify the variables for which the p-values should be those of exact tests. By default all p-values are from large sample approximation tests (chisq.test)., Default: NULL
#' @param includeNA A logic value to indicate if NA should be handled as a regular factor level rather than missing. NA is shown as the last factor level in the table if \code{includeNA = T}. Only effective for categorical variables., Default: FALSE
#' @param doWord A logic value to indicate if transformating the .txt outfile to a .docx WORD file (.txt file will be also kept).
#' @param tab.name A string value to indicate the name of the output table.
#' @param res.path A string value to indicate the path for saving the table.
#' @param ... Additionnal parameters pass to jstable::CreateTableOne2
#' @importFrom jstable CreateTableOne2
#' @importFrom dplyr %>%
#' @importFrom officer read_docx body_add_par body_add_table body_add_par
#' @export
#' @return A summarizing table with stastitical testing results.
#' @examples # There is no example and please refer to vignette.
compClinvar <- function(moic.res      = NULL,
                        var2comp      = NULL,
                        strata        = "Subtype",
                        factorVars    = NULL,
                        nonnormalVars = NULL,
                        exactVars     = NULL,
                        includeNA     = F,
                        doWord        = T,
                        tab.name      = NULL,
                        res.path      = getwd(),
                        ...){

  dat <- moic.res$clust.res
  colnames(dat)[2] <- "Subtype"
  dat$Subtype <- paste0("CS", dat$Subtype)
  com_sam <- intersect(dat$samID,rownames(var2comp))

  # check data
  if(length(com_sam) == nrow(dat)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(dat)-length(com_sam))," samples mismatched from current subtypes."))
  }
  dat <- cbind.data.frame("Subtype" = dat[com_sam, "Subtype", drop = F], var2comp[com_sam, , drop = F])

  # summarizing
  stabl <- jstable::CreateTableOne2(vars = setdiff(colnames(dat), strata),
                                    strata = strata,
                                    data = dat,
                                    factorVars = factorVars,
                                    nonnormal = nonnormalVars,
                                    exact = exactVars,
                                    includeNA = includeNA,
                                    showAllLevels = T,
                                    ...)
  comtable <- as.data.frame(stabl)
  comtable <- cbind.data.frame(var = rownames(stabl), comtable)
  rownames(comtable) <- NULL; colnames(comtable)[1] <- " "
  comtable[is.na(comtable)] <- ""
  #print(comtable)

  if(is.null(tab.name)) {
    outFile <- "summarization of clinical variables stratified by current subtype.txt"
  } else {
    outFile <- paste0(tab.name,".txt")
  }

  write.table(stabl, file.path(res.path,outFile), sep = "\t", quote = F)

  # generate WORD format
  if(doWord){

    table_subtitle <- colnames(comtable)

    title_name <- paste0("Table *. ", gsub(".txt", "", outFile, fixed = T))
    mynote <- "Note: ..."

    my_doc <- officer::read_docx()
    my_doc %>%
      officer::body_add_par(value = title_name, style = "table title") %>%
      officer::body_add_table(value = comtable, style = "table_template") %>%
      officer::body_add_par(value = mynote) %>%
      print(target = file.path(res.path,paste0("TABLE ", gsub(".txt", "", outFile, fixed = T), ".docx")))
  }
  return(list(compTab = comtable))
}
