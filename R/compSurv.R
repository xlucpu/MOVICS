#' @name compSurv
#' @title Comparison of prognosis by Kaplan-Meier survival curve
#' @description This function calculates Kaplan-meier estimator and generate survival curves with log-rank test to detect prognostic difference among identified subtypes. If more than 2 subtypes are identified, pair-wise comparisons will be performed with an additional table printed on the survival curve.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param surv.info A data.frame with rownames of observations and with at least two columns of `futime` for survival time and `fustat` for survival status (0: censoring; 1: event)
#' @param clust.col A string vector storing colors for each subtype.
#' @param p.adjust.method A string value for indicating method for adjusting p values (see \link[stats]{p.adjust}). Allowed values include one of c(`holm`, `hochberg`, `hommel`, `bonferroni`, `BH`, `BY`, `fdr`, `none`); "BH" by default.
#' @param fig.name A string value to indicate the output path for storing the kaplan-meier curve.
#' @param fig.path A string value to indicate the name of the kaplan-meier curve.
#' @param convt.time A string value to indicate how to convert the survival time; value of `d` for days, `m` for months and `y` for years; "d" by default.
#' @param surv.median.line A string value for drawing a horizontal/vertical line at median survival. Allowed values include one of c(`none`, `hv`, `h`, `v`). v: vertical, h:horizontal; "none" by default.
#' @param xyrs.est An integer vector to estimate probability of surviving beyond a certain number (x) of years (Estimating x-year survival); NULL by default.
#' @param surv.cut A numeric value to indicate the x-axis cutoff for showing the maximal survival time. NULL by default (show 0-maximum survival time range).
#' @return A figure of multi-omics Kaplan-Meier curve (.pdf) and a list with the following components:
#'
#'         \code{fitd}       an object returned by \link[survival]{survdiff}.
#'
#'         \code{fid}        an object returned by \link[survival]{survfit}.
#'
#'         \code{xyrs.est}   x-year probability of survival and the associated lower and upper bounds of the 95% confidence interval are also displayed if argument of `xyrs.est` was set by users.
#'
#'         \code{overall.p}  a nominal p.value calculated by Kaplain-Meier estimator with log-rank test.
#'
#'         \code{pairwise.p} an object of class "pairwise.htest" which is a list containing the p values (see \link[survminer]{pairwise_survdiff}); (\emph{only returned when more than 2 subtypes are identified}).
#' @import survival
#' @import survminer
#' @import ggplot2
#' @importFrom grDevices pdf dev.off pdf.options
#' @importFrom ggpp geom_table
#' @importFrom tibble tibble
#' @export
#' @examples # There is no example and please refer to vignette.
compSurv <- function(moic.res         = NULL,
                     surv.info        = NULL,
                     convt.time       = "d",
                     surv.cut         = NULL,
                     xyrs.est         = NULL,
                     clust.col        = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                     p.adjust.method  = "BH",
                     surv.median.line = "none",
                     fig.name         = NULL,
                     fig.path         = getwd()){

  if(!all(is.element(c("futime","fustat"), colnames(surv.info)))) {
    stop("fail to find variables of futime and fustat.")
  }

  if(!all(is.element(convt.time, c("d","m","y")))) {
    stop("unsupported time conversion. Allowed values contain c('d', 'm', 'y').")
  }

  # get common samples
  comsam <- intersect(rownames(surv.info),rownames(moic.res$clust.res))
  mosurv.res <- cbind.data.frame(surv.info[comsam,c("futime","fustat")],
                                 moic.res$clust.res[comsam, "clust", drop = FALSE])
  message(paste0("--a total of ",length(comsam), " samples are identified."))

  # remove missing data if possible
  if(sum(c(is.na(mosurv.res$futime), is.na(mosurv.res$fustat))) > 0) {
    message("--removed missing values.")
    mosurv.res <- as.data.frame(na.omit(mosurv.res))
    message(paste0("--leaving ",nrow(mosurv.res), " observations."))
  }

  # check time unit
  if(max(mosurv.res$futime) < 365) {
    warning("it seems the 'futime' might not in [day] unit, please make sure you provide the correct survival information.")
  }

  # create new variable of Subtype
  mosurv.res$Subtype <- paste0("CS", mosurv.res$clust)
  mosurv.res <- mosurv.res[order(mosurv.res$Subtype),]

  # estimate x-year survival probability
  if(!is.null(xyrs.est)) {
    if(max(xyrs.est) * 365 < max(mosurv.res$futime)) {
      xyrs <- summary(survfit(Surv(futime, fustat) ~ Subtype, data = mosurv.res), times = xyrs.est * 365)
    } else {
      stop("the maximal survival time is less than the time point you want to estimate!")
    }
  } else {
    xyrs <- "[Not Available]: argument of xyrs.est was not specified."
  }

  # convert survival time
  mosurv.res$futime <- switch(convt.time,
                              "d" = mosurv.res$futime,
                              "m" = round(mosurv.res$futime/30.5,4),
                              "y" = round(mosurv.res$futime/365,4))

  date.lab <- switch(convt.time,
                     "d" = "Days",
                     "m" = "Months",
                     "y" = "Years")

  # truncate survival time
  if(date.lab == "Days" & max(mosurv.res$futime) > 3650) {
    #xlim = c(0, 3650)
    brk = 365
  }
  if(date.lab == "Days" & max(mosurv.res$futime) <= 3650) {
    #xlim = c(0, max(mosurv.res$futime))
    brk = floor(max(mosurv.res$futime)/10)
  }
  if(date.lab == "Months" & max(mosurv.res$futime) > 120) {
    #xlim = c(0, 120)
    brk = 12
  }
  if(date.lab == "Months" & max(mosurv.res$futime) <= 120) {
    #xlim = c(0, max(mosurv.res$futime))
    brk = floor(max(mosurv.res$futime)/10)
  }
  if(date.lab == "Years" & max(mosurv.res$futime) > 10) {
    #xlim = c(0, 10)
    brk = 1
  }
  if(date.lab == "Years" & max(mosurv.res$futime) <= 10) {
    #xlim = c(0, max(mosurv.res$futime))
    brk = 1
  }

  if(is.null(surv.cut)) {
    xlim = c(0, max(mosurv.res$futime))
  } else {
    message(paste0("--cut survival curve up to ",surv.cut," ",tolower(date.lab)))
    xlim = c(0, surv.cut)
  }


  # if(!is.null(surv.cut)) {
  #   xlim = c(0, surv.cut)
  #   brk = surv.cut/10
  # } else {message("--cut survival curve up to 10 years.")}

  n.moic <- length(unique(mosurv.res$Subtype))

  # basic survival analysis
  fitd <- survdiff(Surv(futime, fustat) ~ Subtype,
                   data      = mosurv.res,
                   na.action = na.exclude)
  p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
  fit <- survfit(Surv(futime, fustat)~ Subtype,
                 data      = mosurv.res,
                 type      = "kaplan-meier",
                 error     = "greenwood",
                 conf.type = "plain",
                 na.action = na.exclude)

  # hack strata for better survival curve
  names(fit$strata) <- gsub("Subtype=", "", names(fit$strata))

  # kaplan-meier curve
  p <- suppressWarnings(ggsurvplot(fit               = fit,
                                   conf.int          = FALSE,
                                   risk.table        = TRUE,
                                   risk.table.col    = "strata",
                                   palette           = clust.col[1:n.moic],
                                   data              = mosurv.res,
                                   size              = 1,
                                   xlim              = xlim,
                                   break.time.by     = brk,
                                   legend.title      = "",
                                   surv.median.line  = surv.median.line,
                                   xlab              = paste0("Time (",date.lab,")"),
                                   ylab              = "Survival probability (%)",
                                   risk.table.y.text = FALSE))

  # make survival time as percentage
  p$plot <- quiet(p$plot + scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0,100,25)))

  if(n.moic > 2) {

    # add nominal pvalue for log-rank test
    p.lab <- paste0("Overall P",
                    ifelse(p.val < 0.001, " < 0.001",
                           paste0(" = ",round(p.val, 3))))

    p$plot <- p$plot + annotate("text",
                                x = 0, y = 0.55,
                                hjust = 0,
                                fontface = 4,
                                label = p.lab)

    # calculate pair-wise survival comparison
    ps <- pairwise_survdiff(Surv(futime, fustat)~ Subtype,
                            data            = mosurv.res,
                            p.adjust.method = p.adjust.method)

    # add pair-wise comparison table
    # options(stringsAsFactors = FALSE)
    addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
                                             round(ps$p.value, 3))))
    addTab[is.na(addTab)] <- "-"
    # options(stringsAsFactors = TRUE)

    df <- tibble(x = 0, y = 0, tb = list(addTab))
    p$plot <- p$plot + geom_table(data = df, aes(x = x, y = y, label = tb), table.rownames = TRUE)

  } else {
    # add nominal pvalue for log-rank test
    p.lab <- paste0("P",
                    ifelse(p.val < 0.001, " < 0.001",
                           paste0(" = ",round(p.val, 3))))

    p$plot <- p$plot + annotate("text",
                                x = 0, y = 0.55,
                                hjust = 0,
                                fontface = 4,
                                label = p.lab)
  }

  # output curve to pdf
  if(!is.null(fig.name)) {
    outFile <- file.path(fig.path,paste0(fig.name,".pdf"))
  } else {
    outFile <- file.path(fig.path,paste0("km_curve_",moic.res$mo.method,".pdf"))
  }
  pdf.options(reset = TRUE, onefile = FALSE)
  pdf(outFile, width = 6, height = 7)
  # ggsave(outFile, width = 6, height = 7)
  print(p)
  dev.off()

  # output curve to screen
  print(p)

  if(n.moic > 2) {
    return(list(fitd = fitd, fit = fit, xyrs.est = xyrs, overall.p = p.val, pairwise.p = ps))
  } else {
    return(list(fitd = fitd, fit = fit, xyrs.est = xyrs, overall.p = p.val))
  }
}
