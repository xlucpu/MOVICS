#' @name compAgree
#' @title Comparison of agreement between two subtypes
#' @description Compute the Rand Index, Jaccard Index, Fowlkes-Mallows, and Normalized Mutual Information for agreement of two partitions, and generate alluvial diagrams for visualization.
#' @param moic.res An object returned by `getMOIC()` with one specified algorithm or `get\%algorithm_name\%` or `getConsensusMOIC()` with a list of multiple algorithms.
#' @param subt2comp A data.frame of subtypes that need to compare with current subtype with rownames for samples and columns for other subtypes.
#' @param doPlot A logic value to indicate if generating alluvial diagram to show the agreement of different subtypes; TRUE by default.
#' @param box.width A numeric valur to indicate the width for box in alluvial diagram.
#' @param clust.col A string vector storing colors for each cluster.
#' @param width A numeric value to indicate the width of alluvial diagram.
#' @param height A numeric value to indicate the height of alluvial diagram.
#' @param fig.path A string value to indicate the output path for storing the figure.
#' @param fig.name A string value to indicate the name of the figure.
#' @return A figure of agreement (.pdf) if \code{doPlot = TRUE} and a data.frame storing four agreement measurements, including Rand Index (RI), Adjusted Mutual Information (AMI), Jaccard Index (JI), and Fowlkes-Mallows (FM).
#' @export
#' @import ggplot2
#' @import ggalluvial
#' @importFrom ggalluvial StatStratum
#' @importFrom dplyr group_by tally %>%
#' @importFrom cowplot plot_grid
#' @importFrom flexclust comPart
#' @importFrom aricode AMI ARI
#' @importFrom reshape2 melt
#' @examples # There is no example and please refer to vignette.
compAgree <- function(moic.res  = NULL,
                      subt2comp = NULL,
                      doPlot    = TRUE,
                      clust.col = c("#2EC4B6","#E71D36","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                      box.width = 0.1,
                      fig.name  = NULL,
                      fig.path  = getwd(),
                      width     = 6,
                      height    = 5) {
  dat <- moic.res$clust.res
  colnames(dat)[2] <- "Subtype"
  dat$Subtype <- paste0("CS", dat$Subtype)
  comsam <- intersect(dat$samID,rownames(subt2comp))

  # check data
  if(length(comsam) == nrow(dat)) {
    message("--all samples matched.")
  } else {
    message(paste0("--",(nrow(dat)-length(comsam))," samples mismatched from current subtypes."))
  }
  dat <- cbind.data.frame("Subtype" = dat[comsam, "Subtype", drop = FALSE], subt2comp[comsam, , drop = FALSE])
  dat <- as.data.frame(na.omit(dat))
  if(nrow(dat) != nrow(moic.res$clust.res)) {message("--removed NA values in subt2comp.")}

  var <- colnames(dat)
  n.var <- length(var)
  if(n.var > 6) {stop("please indicate less than 6 subtypes (including current subtypes) that need to compare.")}

  # generate comparsion table
  outTab <- NULL
  c1 <- as.vector(as.numeric(factor(dat[,1])))
  for (i in 2:ncol(dat)) {
    #cat(paste0("Compare ", colnames(dat)[1]," with ", colnames(dat)[i],".\n"))

    c2 <- as.vector(as.numeric(factor(dat[,i])))

    # calculate Rand Index
    RI <- flexclust::comPart(c1, c2, type = c('RI'))

    # calculate Adjusted Mutual Information
    AMI <- aricode::AMI(c1, c2)

    # calculate Jaccard Index
    JI <- flexclust::comPart(c1, c2, type = c('J'))

    # calculate Fowlkes-Mallows
    FM <- flexclust::comPart(c1, c2, type = c('FM'))

    outTab <- rbind.data.frame(outTab,data.frame(current.subtype = colnames(dat)[1],
                                                 other.subtype = colnames(dat)[i],
                                                 RI = as.numeric(RI),
                                                 AMI = as.numeric(AMI),
                                                 JI = as.numeric(JI),
                                                 FM = as.numeric(FM),
                                                 stringsAsFactors = FALSE),
                               stringsAsFactors = FALSE)
  }

  assign("StatStratum", ggalluvial::StatStratum, envir=globalenv())

  if(doPlot) {

    if(is.null(fig.name)) {
      outFig <- "Agreement between current subtype and other classifications.pdf"
    } else {
      outFig <- paste0(fig.name,".pdf")
    }

    # generate barplot for agreement
    agreement <- reshape2::melt(outTab[,2:ncol(outTab)], id.vars = "other.subtype", variable.name = "Method")
    b <- ggplot(data = agreement, aes(x = Method, y = value, fill = other.subtype)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_fill_brewer(palette = "Set1") + ggplot2::labs(x = "", y = "Scalar") +
      scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(color = "black", size = 12, face = "bold", vjust = -5),
            axis.text.y = element_text(color = "black", size = 12, face = "bold"),
            axis.title.x = element_text(color = "black", size = 12, face = "bold"),
            axis.title.y = element_text(color = "black", size = 12, face = "bold")) +
      ggtitle("")

    # generate alluvial diagram
    col = clust.col[1:length(unique(dat$Subtype))]
    var1 <- var[1]

    # 1 subtype
    if(n.var == 2) {
      var2 <- var[2]
      subdf <- dat[,1:n.var]
      colnames(subdf) <- c("Subtype",paste0("Subtype",1:(n.var-1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1) %>% tally(name = "Freq") %>% as.data.frame()

      p <- ggplot(subdf,
                  aes(y = Freq,
                      axis1 = Subtype,
                      axis2 = Subtype1)) +
        scale_fill_manual(values = col) +
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) +
        geom_stratum(width = 1/8, reverse = TRUE) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = TRUE) +
        scale_x_continuous(breaks = 1:n.var, labels = c(var1, var2)) +
        theme_bw() +
        theme(legend.position = "top",
              legend.title = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
        ggtitle("")
    }

    # 2 subtypes
    if(n.var == 3) {
      var2 <- var[2]
      var3 <- var[3]
      subdf <- dat[,1:n.var]
      colnames(subdf) <- c("Subtype",paste0("Subtype",1:(n.var-1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1, Subtype2) %>% tally(name = "Freq") %>% as.data.frame()

      p <- ggplot(subdf,
                  aes(y = Freq,
                      axis1 = Subtype,
                      axis2 = Subtype1,
                      axis3 = Subtype2)) +
        scale_fill_manual(values = col) +
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) +
        geom_stratum(width = box.width, reverse = TRUE) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = TRUE) +
        scale_x_continuous(breaks = 1:n.var, labels = c(var1, var2, var3)) +
        theme_bw() +
        theme(legend.position = "top",
              legend.title = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
        ggtitle("")
    }

    # 3 subtypes
    if(n.var == 4) {
      var2 <- var[2]
      var3 <- var[3]
      var4 <- var[4]
      subdf <- dat[,1:n.var]
      colnames(subdf) <- c("Subtype",paste0("Subtype",1:(n.var-1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1, Subtype2, Subtype3) %>% tally(name = "Freq") %>% as.data.frame()

      p <- ggplot(subdf,
                  aes(y = Freq,
                      axis1 = Subtype,
                      axis2 = Subtype1,
                      axis3 = Subtype2,
                      axis4 = Subtype3)) +
        scale_fill_manual(values = col) +
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) +
        geom_stratum(width = 1/8, reverse = TRUE) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = TRUE) +
        scale_x_continuous(breaks = 1:n.var, labels = c(var1, var2, var3, var4)) +
        theme_bw() +
        theme(legend.position = "top",
              legend.title = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
        ggtitle("")
    }

    # 4 subtypes
    if(n.var == 5) {
      var2 <- var[2]
      var3 <- var[3]
      var4 <- var[4]
      var5 <- var[5]
      subdf <- dat[,1:n.var]
      colnames(subdf) <- c("Subtype",paste0("Subtype",1:(n.var-1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1, Subtype2, Subtype3, Subtype4) %>% tally(name = "Freq") %>% as.data.frame()

      p <- ggplot(subdf,
                  aes(y = Freq,
                      axis1 = Subtype,
                      axis2 = Subtype1,
                      axis3 = Subtype2,
                      axis4 = Subtype3,
                      axis5 = Subtype4)) +
        scale_fill_manual(values = col) +
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) +
        geom_stratum(width = 1/8, reverse = TRUE) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = TRUE) +
        scale_x_continuous(breaks = 1:n.var, labels = c(var1, var2, var3, var4, var5)) +
        theme_bw() +
        theme(legend.position = "top",
              legend.title = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
        ggtitle("")
    }

    # 5 subtypes
    if(n.var == 6) {
      var2 <- var[2]
      var3 <- var[3]
      var4 <- var[4]
      var5 <- var[5]
      var6 <- var[6]
      subdf <- dat[,1:n.var]
      colnames(subdf) <- c("Subtype",paste0("Subtype",1:(n.var-1)))
      subdf <- subdf %>% group_by(Subtype, Subtype1, Subtype2, Subtype3, Subtype4, Subtype5) %>% tally(name = "Freq") %>% as.data.frame()

      p <- ggplot(subdf,
                  aes(y = Freq,
                      axis1 = Subtype,
                      axis2 = Subtype1,
                      axis3 = Subtype2,
                      axis4 = Subtype3,
                      axis5 = Subtype4,
                      axis6 = Subtype5)) +
        scale_fill_manual(values = col) +
        geom_flow(stat = "alluvium", width = 1/8, aes(fill = Subtype)) +
        geom_stratum(width = 1/8, reverse = TRUE) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = TRUE) +
        scale_x_continuous(breaks = 1:n.var, labels = c(var1, var2, var3, var4, var5, var6)) +
        theme_bw() +
        theme(legend.position = "top",
              legend.title = element_blank(),
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
        ggtitle("")
    }
    agree <- list(b,p)
    bp <- plot_grid(plotlist = agree, ncol = 2)

    # save to pdf
    ggsave(file.path(fig.path,outFig), width = width, height = height)

    # output to screen
    print(bp)
  }
  return(outTab)
}
