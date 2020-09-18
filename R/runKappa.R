#' @name runKappa
#' @title Run consistency evaluation using Kappa statistic
#' @description Calculate Kappa statistic to measure the consistency between two appraisements
#' @param subt1 A numeric vector to indicate the first appraisement. Order should be exactly the same with subt2 for each sample.
#' @param subt2 A numeric vector to indicate the second appraisement. Order should be exactly the same with subt1 for each sample.
#' @param subt1.lab A string value to indicate the label of the first subtype.
#' @param subt2.lab A string value to indicate the label of the second subtype.
#' @param fig.path A string value to indicate the output path for storing the consistency heatmap.
#' @param fig.name A string value to indicate the name of the consistency heatmap.
#' @param width A numeric value to indicate the width of output figure.
#' @param height A numeric value to indicate the height of output figure.
#' @details This function evaluates the consistency between two appraisements that targets to the same cohort.
#' For example, the NTP-predicted subtype amd PAM-predicted subtype of external cohort, or the current subtype and predicted subtype of discovery cohort.
#' Therefore, the arguments `subt1` and `subt2` can be the `clust` column of `clust.res` derived from `getMOIC()` with one specified algorithm
#' or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms or `runNTP()` or `runPAM()`.
#' However, subtypes identified from different algorithm (i.e., `get\%algorithm_name1\%` and `get\%algorithm_name2\%`) can not be evaluated
#' because the subtype 1 identified from the first algorithm may not be the same subtype 1 from the second algorithm.
#' @return A figure of consistency heatmap (.pdf).
#' @export
#' @importFrom grDevices dev.copy2pdf
#' @examples # There is no example and please refer to vignette.
runKappa <- function(subt1     = NULL,
                     subt2     = NULL,
                     subt1.lab = NULL,
                     subt2.lab = NULL,
                     fig.path  = getwd(),
                     fig.name  = "constheatmap",
                     width     = 5,
                     height    = 5) {

  subt1 <- as.vector(as.character(subt1))
  subt2 <- as.vector(as.character(subt2))

  # check sample
  if(length(subt1) != length(subt2)) {
    stop("subtypes identified from different cohorts.")
  }

  # check subtype
  if(!identical(sort(unique(subt1)),sort(unique(subt2)))) {
    stop("subtypes fail to matched from two appraisements.")
  }

  # check argument
  if(is.null(subt1.lab) | is.null(subt2.lab)) {
    stop("label for subtype1 and subtype2 must be both indicated.")
  }

  # calculate consistency
  comb.subt <- data.frame(subt1 = paste0("CS", subt1),
                          subt2 = paste0("CS", subt2),
                          stringsAsFactors = F)

  tab_classify <- as.data.frame.array(table(comb.subt$subt1,comb.subt$subt2))

  # n.moic <- length(unique(moic.res$clust.res$clust))
  # if(n.moic < 5) {
  #   p.fisher <- fisher.test(tab_classify, workspace = 1e9)$p.value
  # } else {
  #   message("--using simulated p value in fisher's exact test.")
  #   p.fisher <- fisher.test(tab_classify, simulate.p.value = T)$p.value
  # }

  # calculate kappa statistic
  x <- table(comb.subt$subt1,comb.subt$subt2)
  nr <- nrow(x); nc <- ncol(x); N <- sum(x)
  Po <- sum(diag(x))/N; Pe <- sum(rowSums(x) * colSums(x)/N)/N
  kappa <- (Po - Pe)/(1 - Pe)
  seK0 <- sqrt(Pe/(N * (1 - Pe)))
  p.v <- 1 - pnorm(kappa/seK0)
  p.lab <- ifelse(p.v < 0.001, "P < 0.001", paste0("P = ", format(round(p.v,3), nsmall = 3)))

  # generate consistency table
  blue   <- "#204F8D"
  lblue  <- "#498EB9"
  dwhite <- "#B6D1E8"
  white  <- "#E6EAF7"

  par(bty="n", mgp = c(2,0.5,0), mar = c(4.1,4.1,4.1,2.1),tcl=-.25, font.main=3)
  par(xpd=NA)
  plot(c(0,ncol(tab_classify)),c(0,nrow(tab_classify)),
       col = "white",
       xlab = "",xaxt = "n",
       ylab = "",yaxt = "n")
  title(paste0("Consistency between ",subt1.lab," and ",subt2.lab,"\nKappa = ", format(round(kappa,3), nsmall = 3),
               "\n", p.lab),
        adj = 0, line = 0)

  # add y-axis
  axis(2, at = 0.5:(nrow(tab_classify)-0.5), labels = FALSE)
  text(y = 0.5:(nrow(tab_classify)-0.5),
       par("usr")[1],
       labels = rownames(tab_classify)[nrow(tab_classify):1],
       srt = 0, pos = 2, xpd = TRUE)
  mtext(paste0("Subtypes derived from ", subt1.lab), side = 2, line = 3)

  # add x-axis
  axis(1, at = 0.5:(ncol(tab_classify)-0.5), labels = FALSE)
  text(x = 0.5:(ncol(tab_classify)-0.5),
       par("usr")[1] - 0.2,
       labels = colnames(tab_classify),
       srt = 45, pos = 1, xpd = TRUE)
  mtext(paste0("Subtypes derived from ", subt2.lab), side = 1, line = 3)

  # generate colors
  input_matrix <- as.matrix(tab_classify)
  mat.max = max(input_matrix)
  unq.value <- unique(sort(as.vector(input_matrix)))
  rbPal <- colorRampPalette(c(white,dwhite,lblue,blue))
  col.vec <- rbPal(max(unq.value) + 1)
  col.mat <- matrix(NA,byrow = T,ncol = ncol(input_matrix),nrow = nrow(input_matrix))

  # fill matrix
  for (i in 1:nrow(col.mat)) {
    for (j in 1:ncol(col.mat)) {
      col.mat[i,j] <- col.vec[input_matrix[i,j] + 1]
    }
  }

  # generate heatmap
  x_size <- ncol(input_matrix)
  y_size <- nrow(input_matrix)

  my_xleft = rep(c(0:(x_size-1)),each = x_size)
  my_xright = my_xleft + 1
  my_ybottom = rep(c((y_size-1):0),y_size)
  my_ytop = my_ybottom + 1
  rect(xleft = my_xleft,
       ybottom = my_ybottom,
       xright = my_xright,
       ytop = my_ytop,
       col=col.mat,
       border = F)

  # fill count
  text(my_xleft + 0.5,my_ybottom + 0.5,input_matrix, cex = 1.3)

  # output to pdf
  outFig <- paste0(fig.name,".pdf")
  invisible(dev.copy2pdf(file = file.path(fig.path, outFig), width = width, height = height))
}
