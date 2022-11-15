#' unite neighbouring segments if possible
#'
#' @param seg segmentation data

#Shrink each sample by Chromosome
shrink.seg.ai.wrapper <- function(seg){
  new.seg <- seg[1,]
  #For each of  the samples
  for(j in unique(seg[,1])){
    sample.seg <- seg[seg[,1] %in% j,]
    new.sample.seg <- seg[1,]
    #For each of the chromosomes
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]

      #Just make shrink for each chromosomes, with more than two segments.
      if(nrow(sample.chrom.seg) > 1){
        sample.chrom.seg <- shrink.seg.ai(sample.chrom.seg)
      }
      #Concatenate the rows, of each chromosome
      new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
    }
    #Concatenate the rows of each sample
    new.seg <- rbind(new.seg, new.sample.seg[-1,])
  }
  seg <- new.seg[-1,]
  return(seg)
}


#Shrink each sample by p.arm and q.arm of each chromosome
shrink.seg.ai.wrapper2 <- function(seg){
  new.seg <- seg[1,]
  #For each of  the samples
  for(j in unique(seg[,1])){
    sample.seg <- seg[seg[,1] %in% j,]
    new.sample.seg <- seg[1,]
    #For each of the chromosomes
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]

      # split into chromosome arms
      p.arm <- sample.chrom.seg[sample.chrom.seg[,3] <= chrominfo[i,2],]
      q.arm <- sample.chrom.seg[sample.chrom.seg[,4] >= chrominfo[i,3],]
      p.arm[nrow(p.arm),4] <- chrominfo[i,2] #Set that the last p segment ends at the begining of the centromer
      q.arm[1,3] <- chrominfo[i,3] #Set that first q segment, start at the end of the centromer
      p.arm <- shrink.seg.ai(p.arm)
      q.arm <- shrink.seg.ai(q.arm)
      sample.chrom.seg <- rbind(p.arm, q.arm)

      #Concatenate the rows, of each chromosome
      new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
    }
    #Concatenate the rows of each sample
    new.seg <- rbind(new.seg, new.sample.seg[-1,])
  }
  seg <- new.seg[-1,]
  return(seg)
}
