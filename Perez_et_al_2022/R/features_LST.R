######## Extraction of segments sizes ########
#### Also detect the space between segments
segmentSizes <- function(segs){
  #Get segment sizes of segments by Mb
  Segment_sizes <- (segs[,4] - segs[,3])
  #Get size of space between segments
  samples <- unique(segs[,1])
  breaks <- NULL
  for(j in samples){
    sample.seg <- segs[segs[,1] %in% j,]
    chroms <- unique(sample.seg[,2])
    chroms <- chroms[!chroms %in% c(23,24,'chr23','chr24','chrX','chrx','chrY','chry')]
    for(chr in chroms){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% chr,,drop=F]
      if(nrow(sample.chrom.seg) < 2) {next}
      for (id in 2:nrow(sample.chrom.seg)){
        space <- (sample.chrom.seg[id,3] - sample.chrom.seg[id-1,4])
        breaks <- c(breaks, space)
      }
    }
  }
  Segment_sizes <- c(Segment_sizes, breaks)
  return(Segment_sizes)
}



######## Extraction of events sizes by sample ########
###And also detect the proportion of segments equal or greather than determined sizes
### As output a dataframe with proportion of segments equal or greather than the MbSizes
segmentBySamples <- function(segs, MbSizes = c(0,1,4,7,10,13,16,19,22,25,28,31,34,37,40)){
  Proportions <- NULL
  #Get size of space between segments
  samples <- unique(segs[,1])
  breaks <- NULL
  for(j in samples){
    sample.seg <- segs[segs[,1] %in% j,]
    #Get segment sizes of segments by Mb
    Segment_sizes <- (sample.seg[,4] - sample.seg[,3])
    chroms <- unique(sample.seg[,2])
    chroms <- chroms[!chroms %in% c(23,24,'chr23','chr24','chrX','chrx','chrY','chry')]
    for(chr in chroms){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% chr,,drop=F]
      if(nrow(sample.chrom.seg) < 2) {next}
      for (id in 2:nrow(sample.chrom.seg)){
        space <- (sample.chrom.seg[id,3] - sample.chrom.seg[id-1,4])
        breaks <- c(breaks, space)
      }
    }
    Segment_sizes <- Segment_sizes/1e6 #Convert it to Mb notation
    Segment_proportions_Bysizes <- segmentSizesProportions(Segment_sizes, MbSizes)
    Proportions <- cbind(Proportions, Segment_proportions_Bysizes)
  }
  colnames(Proportions) <- samples
  row.names(Proportions) <- MbSizes
  return(Proportions)
}



###Count state transitions by Segment size
#Segsize is the minimun size of two Allelic Imbalances, continuos
LSTs <- function(seg, chrominfo = chrominfo_grch38, segsizes=10e6, mindistance=3e6, tandemelements=2, getpositions = FALSE){
  samples <- unique(seg[,1])
  output <- setNames(rep(0,length(samples)), samples)
  sample.positions.lst <- NULL
  #Run by each sample
  for(j in samples){
    #print(j)
    sample.seg <- seg[seg[,1] %in% j,]
    sample.lst <- 0

    #Run by each chromosome
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]

      if(nrow(sample.chrom.seg) < tandemelements) {next} #If less than the tandem elements, then next
      sample.chrom.seg.new <- sample.chrom.seg
      # split into chromosome arms
      p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,3] <= chrominfo[i,2],,drop=F]
      q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,4] >= chrominfo[i,3],,drop=F]
      p.arm[nrow(p.arm),4] <- chrominfo[i,2] #Set that the last p segment ends at the begining of the centromer
      q.arm[1,3] <- chrominfo[i,3] #Set that first q segment, start at the end of the centromer

      p.arm<- shrink.seg.ai(p.arm)
      q.arm<- shrink.seg.ai(q.arm)

      #Removing segments less than mindistance for the p.arm, smoothing the resultan values acording to their CNvalues and mindistance
      n.3mb <- which((p.arm[,4] - p.arm[,3]) < mindistance)
      while(length(n.3mb) > 0){
        p.arm <- p.arm[-(n.3mb[1]),]
        p.arm <- shrink.seg.ai(p.arm)
        n.3mb <- which((p.arm[,4] - p.arm[,3]) < mindistance)
      }

      #Removing segments less than mindistance for the q.arm, smoothing the resultan values acording to their CNvalues and mindistance
      n.3mb <- which((q.arm[,4] - q.arm[,3]) < mindistance)
      while(length(n.3mb) > 0){
        q.arm <- q.arm[-(n.3mb[1]),]
        q.arm <- shrink.seg.ai(q.arm)
        n.3mb <- which((q.arm[,4] - q.arm[,3]) < mindistance)
      }

      #Count the number of segments in tandem for the q.arm
      if(nrow(q.arm) >= 2){
        q.arm <- cbind(q.arm[,1:8], c(0,1)[match((q.arm[,4]-q.arm[,3]) >= segsizes, c('FALSE','TRUE'))])
        selected_segments <- which(q.arm[,9] == 1) #Proceed if at least one element is bigger than the segment sizes required
        if (length(selected_segments) >= tandemelements){
          counter <- 0        #This will count the number of consecutive events
          positions <- NULL   #Store the row number of the elements in tandem
          for (l in 1:length(q.arm[,9])){
            if (q.arm[l,9] == 1){
              counter <- counter + 1
              positions <- c(positions,l)
              if (counter == tandemelements){
                  counter <- counter - 1 #Rest one in case of tandem events
                  #Now we check all the spaces between segments are less than mindistance
                  condition <- 1
                  for (m in positions[-1]){
                    if ((q.arm[m,3] - q.arm[m-1,4]) >= mindistance){
                      condition <- 0
                    }
                  }
                  positions <- positions[-1] #Remove the first element in case of tandem events
                  if (condition == 1){
                      sample.lst <- sample.lst + 1
                      events <- q.arm[positions,c(1:8)]
                      sample.positions.lst <- rbind(sample.positions.lst, events)
                  }
              }
            }else{
              counter = 0
              positions <- NULL
            }
          }
        }
      }

      #Count the number of segments in tandem for the p.arm
      if(nrow(p.arm) >= 2){
        p.arm <- cbind(p.arm[,1:8], c(0,1)[match((p.arm[,4] - p.arm[,3]) >= segsizes, c('FALSE','TRUE'))])
        selected_segments <- which(p.arm[,9] == 1) #Proceed if at least one element is bigger than the segment sizes required
        if (length(selected_segments) >= tandemelements){
          counter <- 0        #This will count the number of consecutive events
          positions <- NULL   #Store the row number of the elements in tandem
          for (l in 1:length(p.arm[,9])){
            if (p.arm[l,9] == 1){
              counter <- counter + 1
              positions <- c(positions,l)
              if (counter == tandemelements){
                counter <- counter - 1 #Rest one in case of tandem events
                #Now we check all the spaces between segments are less than mindistance
                condition <- 1
                for (m in positions[-1]){
                  if ((p.arm[m,3] - p.arm[m-1,4]) >= mindistance){
                    condition <- 0
                  }
                }
                positions <- positions[-1] #Remove the first element in case of tandem events
                if (condition == 1){
                  sample.lst <- sample.lst + 1
                  events <- p.arm[positions, c(1:8)]
                  sample.positions.lst <- rbind(sample.positions.lst, events)
                }
              }
            }else{
              counter = 0
              positions <- NULL
            }
          }
        }
      }

    } #End of for by chr

    output[j] <- sample.lst
  }
  if (getpositions){
    return(sample.positions.lst)
  }else{
    return(output)
  }

}
