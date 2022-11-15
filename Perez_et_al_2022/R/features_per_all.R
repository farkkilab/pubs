#' Extract features per samples
#'
#' @param seg segmentation data

######### Function to extract the number of LOH events of a particular size ########
#Start allways MbSizes with 0
#As output will give default output will give you the LOH events of sizes 0 to 1MB, 1Mb to 4Mb, 7Mb to 10Mb, and so on... according to the windows sizes selected
features.LOH <- function(segs,  MbSizes = c(0,1,4,7,10,13,16,19,22,25,28,31,34,37,40), A_cn=7, B_cn=8){
    #For each of  the samples
    totalCountLOHs <- NULL
    for(sample in unique(segs[,1])){
        sample.seg <- segs[segs[,1] %in% sample,]
        #Selection of LOH events
        segLOH <- sample.seg[sample.seg[,B_cn] == 0 & sample.seg[,A_cn] != 0,,drop=F]
        sizesLOH <- (segLOH[,4] - segLOH[,3])/1e6
        countLoHSample <- NULL
        for (i in 2:length(MbSizes)){
              countLOHs <- length(which(sizesLOH >= MbSizes[i-1] & sizesLOH <= MbSizes[i]))
              countLoHSample <- c(countLoHSample, countLOHs)
        }
        countLOHs <- length(which(sizesLOH > MbSizes[i]))
        countLoHSample <- c(countLoHSample, countLOHs)
        totalCountLOHs <- rbind(totalCountLOHs, countLoHSample)
    }
    row.names(totalCountLOHs) <- unique(segs[,1])
    colnames(totalCountLOHs) <- c(paste("LoH", paste(MbSizes, "Mb", sep="" ), sep="_"))
    #colnames(totalCountLOHs) <- c(paste(MbSizes, "Mb", sep="" ))
    totalCountLOHs <- as.data.frame(totalCountLOHs)
    return(totalCountLOHs)
}


###Segments by proportions2, select for those bins greater than the actual bin but smaller than the next bin (segment size in MbSizes)
segmentSizesProportions <- function(Segment_sizes, MbSizes = c(0,1,4,7,10,13,16,19,22,25,28,31,34,37,40)){
    number_events <- NULL
    for (i in 1:length(MbSizes)){
        if (i == length(MbSizes)){ #If its the last one
            x <- length(which(Segment_sizes >= MbSizes[i])) #Get the number of events greater or equal to the size bin
        }else{
            x <- length(which(Segment_sizes >= MbSizes[i] & Segment_sizes < MbSizes[i+1])) #Get the number of events greater or equal to the size bin
        }
        number_events <- c(number_events, x)
    }
    names(number_events) <- paste(MbSizes,"Mb",sep="")
    proportion_events <- number_events/length(Segment_sizes)
    return(proportion_events)
}


################ Set an dispertion accuracy #########################
#Take as input two vectors, and his status
#Get a cutoff that separe the categories using desition tree
#Then check the number of values that are above the cuts, but should be below
#Calculate ACC using this information

get_dispertion_acc <- function(values1, values2, status1="HRD", status2="HRP"){

        x <- data.frame(LSTS=values1, anyvalue=rep(2,length(values1)), status=rep(status1,length(values1)))
        y <- data.frame(LSTS=values2, anyvalue=rep(2,length(values2)), status=rep(status2,length(values2)))
        z <- rbind(x,y)

        #Select a cutoff using desition tree
        fit <- rpart(status ~ ., data=z, method='class', control=rpart.control(minsplit = 2, minbucket = 1, cp=0.001))

        #Select
        value <- rpart.subrules.table(fit)[1,5]
        cutoff <- as.numeric(value)

        TP <- length(which(values1 >= cutoff))
        FN <- length(which(values1 < cutoff))
        FP <- length(which(values2 >= cutoff))
        TN <- length(which(values2 < cutoff))
        TPR <- TP/(TP+FN)
        TNR <- TN/(TN+FP)
        BA <- (TPR + TNR)/2
        return(BA)
}
