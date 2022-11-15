#' Calculate U test between two dataframes of counts
#'
#' @param seg segmentation data

U.test <- function(seg1,seg2){
    pvals <- NULL
    mean.fold.changes <- NULL
    means1 <- NULL
    medians1 <- NULL
    medians2 <- NULL
    means2 <- NULL
    for (i in 1:ncol(seg1)){
        average1 <- mean(seg1[,i])
        average2 <- mean(seg2[,i])
        median1 <- median(seg1[,i])
        median2 <- median(seg2[,i])
        Utest <- wilcox.test(seg1[,i],seg2[,i], alternative="greater", paired=FALSE)
        mean.fold <-  mean(seg1[,i]) / mean(seg2[,i]) - 1
        medians1 <- c(medians1, median1)
        medians2 <- c(medians2, median2)
        means1 <- c(means1, average1)
        means2 <- c(means2, average2)
        pvals <- c(pvals, Utest$p.value)
        mean.fold.changes <- c(mean.fold.changes, mean.fold)
    }

    window.difference<- data.frame(MbSizes=colnames(seg1),
                                   p.values=pvals,
                                   fold.change=mean.fold.changes,
                                   Mean.events.1 = means1,
                                   Mean.events.2 = means2,
                                   Median.events.1 = medians1,
                                   Median.events.2 = medians2)
    return(window.difference)
}
