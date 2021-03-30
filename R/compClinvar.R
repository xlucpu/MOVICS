#' @name compClinvar
#' @title Comparison of clinical variables
#' @description Create a table summarizing all baseline variables (both continuous and categorical) stratifying by current identified Subtypes and performing statistical tests. The object gives a table that is easy to use in medical research papers.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param var2comp A data.frame of clinical variables that need to compare among current subtypes with rownames for samples and columns for variable names.
#' @param strata A string value to indicate the stratifying variable. This function will generate an internal 'Subtype' variable which concatenates a string of 'CS' and values from 'clust' column of 'clust.res' in argument of `moic.res`. This argument is set as NULL by default and in this case using 'Subtype' variable as strata.
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
                        strata        = NULL,
                        factorVars    = NULL,
                        nonnormalVars = NULL,
                        exactVars     = NULL,
                        includeNA     = FALSE,
                        doWord        = TRUE,
                        tab.name      = NULL,
                        res.path      = getwd(),
                        ...){

  dat <- moic.res$clust.res
  colnames(dat)[which(colnames(dat) == "clust")] <- "Subtype"
  dat$Subtype <- paste0("CS", dat$Subtype)
  com_sam <- intersect(dat$samID,rownames(var2comp))

  # check data
  if(length(com_sam) == nrow(dat)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(dat)-length(com_sam))," samples mismatched from current subtypes."))
  }
  dat <- cbind.data.frame("Subtype" = dat[com_sam, "Subtype", drop = FALSE], var2comp[com_sam, , drop = FALSE])

  # summarizing
  if(is.null(strata)) {
    strata <- "Subtype"
  }
  if(!is.element(strata, colnames(dat))) {
    stop("fail to find this strata in var2comp. Consider using NULL by default.")
  }

  warn <- NULL
  tryCatch(stabl <-
             jstable::CreateTableOne2(vars = setdiff(colnames(dat), strata),
                                     strata = strata,
                                     data = dat,
                                     factorVars = factorVars,
                                     nonnormal = nonnormalVars,
                                     exact = exactVars,
                                     includeNA = includeNA,
                                     showAllLevels = TRUE,
                                     ...),
           warning=function(w) {
             warn <<- append(warn, conditionMessage(w))}
           )
  if(grepl("NA", warn, fixed = TRUE)) { # if get warning, probably due to the failure of computing exact p value in large data
    set.seed(19991018)
    stabl <- jstable::CreateTableOne2(vars = setdiff(colnames(dat), strata),
                                      strata = strata,
                                      data = dat,
                                      factorVars = factorVars,
                                      nonnormal = nonnormalVars,
                                      exact = exactVars,
                                      includeNA = includeNA,
                                      showAllLevels = TRUE,
                                      argsExact = list(simulate.p.value = T), # use simulated P values then
                                      ...)
  } else {
    stabl <- jstable::CreateTableOne2(vars = setdiff(colnames(dat), strata),
                                      strata = strata,
                                      data = dat,
                                      factorVars = factorVars,
                                      nonnormal = nonnormalVars,
                                      exact = exactVars,
                                      includeNA = includeNA,
                                      showAllLevels = TRUE,
                                      ...)
  }

  # trim output
  comtable <- as.data.frame(stabl)
  comtable <- cbind.data.frame(var = rownames(stabl), comtable)
  rownames(comtable) <- NULL; colnames(comtable)[1] <- " "
  comtable[is.na(comtable)] <- ""
  comtable <- comtable[,setdiff(colnames(comtable),"sig")] # remove significance
  #print(comtable)

  if(is.null(tab.name)) {
    outFile <- "summarization of clinical variables stratified by current subtype.txt"
  } else {
    outFile <- paste0(tab.name,".txt")
  }

  write.table(stabl, file.path(res.path,outFile), sep = "\t", quote = FALSE)

  # generate WORD format
  if(doWord){
    table_subtitle <- colnames(comtable)

    title_name <- paste0("Table *. ", gsub(".txt", "", outFile, fixed = TRUE))
    mynote <- "Note: ..."

    my_doc <- officer::read_docx()
    my_doc %>%
      officer::body_add_par(value = title_name, style = "table title") %>%
      officer::body_add_table(value = comtable, style = "table_template") %>%
      officer::body_add_par(value = mynote) %>%
      print(target = file.path(res.path,paste0("TABLE ", gsub(".txt", "", outFile, fixed = TRUE), ".docx")))
  }
  return(list(compTab = comtable))
}
