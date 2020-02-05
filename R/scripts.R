#' Calculate value N std dev away from mean
#'
#' This function allows you to identify genes with loss of expression
#' @param mn Mean of distribution
#' @param sd std dev of distribution
#' @param n number of std dev below mean to calculate
#' @return the value `n` standard deviations below the mean `mn`
#' @keywords internal

nSdBelowMean <- function(mn, sd, n){
    nsdblw <- mn - n*sd
    return(nsdblw)
}

#' Identify genes with expression loss
#'
#' This function allows you to identify genes with loss of expression
#' @param exprMatrNml A matrix of expression values from normal tissue. Each row
#'   is a gene, and each column is a patient or sample. Genes should be in same
#'   order as exprMatrTum.
#' @param exprMatrTum A matrix of expression values from tumor tissue. Each row
#'   is a gene, and each column is a patient or sample. Genes should be in same
#'   order as exprMatrNml.
#' @param nSdBelow The number of SD below the mean of the adjacent normal tissue
#'   to set the boundary between tumor subgroups.
#' @param minPropPerGroup A value between 0-1 that represents the minimum
#' proportion of samples that should be present in each of the two subgroups
#' (defined by the boundary set by nSdBelow) for a particular gene.
#' @return a nx7 matrix, with n equaling the number of genes. The columns are
#' as follows:
#' \itemize{
#'   \item geneNm - the gene name
#'   \item lowerBound - the lower bound, or the value `nSdBelow` the mean of
#'   the normal tissue expression data.
#'   \item propTumLessThBound - the proportion of tumor samples with
#'   expression levels less than `lowerBound`
#'   \item muAb - "mu above", the mean expression value of tumors greater than
#'   (ie above) the `lowerBound`.
#'   \item `muBl` - "mu below", the mean expression value of tumors less than
#'   (ie below) the `lowerBound`.
#'   \item `deltaMu` - the difference between `muAb` and `muBl`.
#'   \item meetsMinPropPerGrp - a logical indicating whether the proportion
#'   of samples in each group is greater than that set by `minPropPerGroup`.
#' }
#' @keywords gene expression, subgroups
#' @importFrom magrittr "%>%"
#' @import tidyr tidyverse stats
#' @export
#' @examples
#' exprMatrNml <- matrix(abs(rnorm(100, mean=2)), nrow=10)
#' exprMatrTum <- matrix(abs(rnorm(100)), nrow=10)
#' geneNames <- paste0(letters[1:nrow(exprMatrNml)], 1:nrow(exprMatrNml))
#' rownames(exprMatrNml) <- rownames(exprMatrTum) <- geneNames
#' nSdBelow <- 2
#' minPropPerGroup <- .2
#' rl <- receptLoss(exprMatrNml, exprMatrTum, nSdBelow, minPropPerGroup)
#' head(rl)

receptLoss <- function(exprMatrNml, exprMatrTum, nSdBelow, minPropPerGroup){
    ## Remove warnings from devtools::check()
    propTumLessThBound <- NULL
    meetsMinPropPerGrp <- NULL

    ## Identify boundary nSdBelow the mean of adj. normal tissue
    boundAll <- apply(exprMatrNml, 1, function(x){
        nSdBelowMean(mean(x), sd(x), nSdBelow)
        }
    )

    propTumLessThanBound <- rowSums(exprMatrTum < boundAll) / ncol(exprMatrTum)

    binIndxMatrGrTh <- ifelse(exprMatrTum > boundAll, 1, NA)
    binIndxMatrLsTh <- ifelse(exprMatrTum < boundAll, 1, NA)
    meansBelow <- rowMeans(exprMatrTum*binIndxMatrLsTh, na.rm=TRUE)
    meansAbove <- rowMeans(exprMatrTum*binIndxMatrGrTh, na.rm=TRUE)
    deltaMu <- meansAbove - meansBelow

    boundAllDf <- dplyr::tibble("geneNm"=rownames(exprMatrNml),
                                "lowerBound"=boundAll,
                                "propTumLessThBound"=propTumLessThanBound,
                                "muAb"=meansAbove,
                                "muBl"=meansBelow,
                                "deltaMu"=deltaMu) %>%
    dplyr::mutate(meetsMinPropPerGrp=propTumLessThBound>minPropPerGroup &
                                propTumLessThBound<{1-minPropPerGroup}) %>%
    dplyr::arrange(-meetsMinPropPerGrp, -deltaMu)
    return(boundAllDf)
}



#' Plot histogram of genes with expression loss
#'
#' This function allows you to plot histograms of tumor and adj normal data
#' @param exprMatrNml A matrix of expression values from normal tissue.
#' Each row is a gene, and each column is a patient or sample. Genes should
#' be in same order as exprMatrTum.
#' @param exprMatrTum A matrix of expression values from tumor tissue. Each row
#'   is a gene, and each column is a patient or sample. Genes should be in same
#'   order as exprMatrNml.
#' @param rldf The dataframe output from running the receptLoss function
#' @param geneName The name of the gene to plot. The name of the gene should
#'   correspond to a row name in both exprMatrNml and exprMatrTum matrices.
#' @param addToTitle A string that can be added to the title, which includes
#'   the gene name.
#' @param clrs Vector of length 2 containing colors to use for plot
#' @return returns an object of class `ggplot`
#' @keywords gene expression, subgroups, visualization
#' @importFrom magrittr "%>%"
#' @import tidyr ggplot2 tidyverse stats
#' @export
#' @examples
#' exprMatrNml <- matrix(abs(rnorm(100, mean=2)), nrow=10)
#' exprMatrTum <- matrix(abs(rnorm(100)), nrow=10)
#' geneNames <- paste0(letters[1:nrow(exprMatrNml)], 1:nrow(exprMatrNml))
#' rownames(exprMatrNml) <- rownames(exprMatrTum) <- geneNames
#' nSdBelow <- 2
#' minPropPerGroup <- .2
#' rl <- receptLoss(exprMatrNml, exprMatrTum, nSdBelow, minPropPerGroup)
#' clrs <- c("#E78AC3", "#8DA0CB")
#' plotReceptLoss(exprMatrNml, exprMatrTum, rl, geneName="g7", clrs=clrs)

plotReceptLoss <- function(exprMatrNml, exprMatrTum, rldf,
                            geneName, addToTitle="", clrs) {
    type <- NULL
    propTumLessThBound <- NULL
    meetsMinPropPerGrp <- NULL
    ..density.. <- NULL

    normal <- exprMatrNml[geneName,] %>%
        data.frame() %>%
        dplyr::mutate(type="normal")
    tumor <- exprMatrTum[geneName,] %>%
        data.frame() %>%
        dplyr::mutate(type="tumor")
    tidyDf <- rbind(normal, tumor) %>%
        dplyr::as_tibble() %>%
        dplyr::rename(expr=".")
    rldf.sub <- rldf[rldf[, "geneNm"] == geneName,]
    p1 <- ggplot(tidyDf, aes(x=expr, color=as.factor(type),
                            fill=as.factor(type), group=as.factor(type)
        )) +
    theme_classic() +
    geom_histogram(
        data=subset(tidyDf, type == "tumor"), alpha=0.5,
        aes(y=..density..), binwidth=0.2,
        color=clrs[2], fill=clrs[2]
    ) +
    scale_fill_manual(values=clrs) +
    theme(
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.position="none",
        text=element_text(size=20, color="black"),
        plot.title=element_text(face="italic")
    ) +
    xlim(c(-0.1, 0.2 + max(exprMatrTum, exprMatrNml))) +
    ggtitle(paste0(geneName, addToTitle)) +
    xlab(expression(Log[2] * "(TPM Reads)")) +
    geom_vline(size=2, alpha=0.8,
        color=clrs[1], xintercept=rldf.sub$lowerBound
    ) +
    stat_function(
        fun="dnorm", colour=clrs[1],
        args=list(mean=mean(normal$.),
                    sd=sd(normal$.)),
        linetype="dashed", size=1.5
    )
    print(p1)
}
