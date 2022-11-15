#' Remove chrm deletions from a segmentation dataset
#'
#' @param seg segmentation data


rm.chr.deletions <- function(seg, B_cn=8){
  newSeg <- NULL
  for(sample in unique(seg[,1])){
      SegSamp <- seg[seg[,1] %in% sample,,drop=F]
      chrDel <- NULL
      for(j in unique(SegSamp[,2])){
        if(all(SegSamp[SegSamp[,2] == j,B_cn] == 0)) {
          chrDel <- c(chrDel, j)
        }
      }
      SegSamp <- SegSamp[!SegSamp[,2] %in% chrDel,,drop=F]
      newSeg <- rbind(newSeg, SegSamp)
  }
  return(newSeg)
}

#Remove chromosome with no deletions inside at all
rm.chr.normals <- function(seg, B_cn=8,  A_cn=7){
  newSeg <- NULL
  for(sample in unique(seg[,1])){
    SegSamp <- seg[seg[,1] %in% sample,,drop=F]
    chrDel <- NULL
    for(j in unique(SegSamp[,2])){
      if(all(SegSamp[SegSamp[,2] == j,B_cn] == SegSamp[SegSamp[,2] == j,A_cn])) {
        chrDel <- c(chrDel, j)
      }
    }
    SegSamp <- SegSamp[!SegSamp[,2] %in% chrDel,,drop=F]
    newSeg <- rbind(newSeg, SegSamp)
  }
  return(newSeg)
}
