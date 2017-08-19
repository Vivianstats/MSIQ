
# strandmode = 0: unstranded
# strandmode = 1: secondstrand
# strandmode = 2, firststrand

library(GenomicAlignments)
library(rbamtools)



get_reads = function(d, cgene, gene_exons, bam_path, strandmode = 0){
  
  genestart = gene_exons$exon_starts[1]
  geneend = gene_exons$exon_ends[cgene$nExons]
  
  # get reads' coordinates from RNA-seq data
  region = GRanges(paste("chr", cgene$chr, sep = ""), 
                   IRanges(genestart, geneend), strand = cgene$str)
  param <- ScanBamParam(which = region)
  readpairs = readGAlignmentPairs(bam_path[d], param = param, strandMode = strandmode)
  if (strandmode !=0 ){
    readpairs = readpairs[strand(readpairs) == cgene$str]}
  
  reads4dim = data.frame(start(GenomicAlignments::first(readpairs)), end(GenomicAlignments::first(readpairs)),
                         start(GenomicAlignments::second(readpairs)), end(GenomicAlignments::second(readpairs)))
  
  # filter out reads mapped to introns
  indNoIntron = sapply(1:4, function(k){
    indexStarts <- findInterval(reads4dim[,k], gene_exons$exon_starts)
    indexEnds <- findInterval(reads4dim[,k], gene_exons$exon_ends + 1)
    
    ind = (indexStarts - 1 == indexEnds)
  })
  if (class(indNoIntron) == "matrix"){
    indNoIntron = rowSums(indNoIntron) == 4
  }else {return(NULL)}
  
  if (sum(indNoIntron) == 0 ) {return(NULL)}
  reads4dim = reads4dim[indNoIntron, , drop = FALSE]
  
  
  # filter out reads with skip on the same exon
  # reads4bin = sapply(1:4, function(k){
  #   findInterval(reads4dim[, k], gene_exons$exon_starts)
  # })
  
  exon_Len = gene_exons$exon_Len
  readsTxcoords = sapply(1:4, function(k){
    indexs = findInterval(reads4dim[, k], gene_exons$exon_starts)
    position = sapply(indexs, function(x) sum(exon_Len[1:x])) +
      reads4dim[, k] - (gene_exons$exon_ends)[indexs]
    return(position)
  })
  if (nrow(reads4dim) == 1){
    readsTxcoords  = matrix(readsTxcoords, nrow = 1)
  }
  
  maxReadLen = max(c(qwidth(first(readpairs)), qwidth(second(readpairs))))
  minReadLen = min(c(qwidth(first(readpairs)), qwidth(second(readpairs))))
  indNoskip = (readsTxcoords[,2] - readsTxcoords[,1] + 1 <= maxReadLen) &
    (readsTxcoords[,4] - readsTxcoords[,3] + 1 <= maxReadLen) &
    (readsTxcoords[,2] - readsTxcoords[,1] + 1 >= minReadLen) &
    (readsTxcoords[,4] - readsTxcoords[,3] + 1 >= minReadLen) 

  if (sum(indNoskip) == 0 ) {return(NULL)}
  readsTxcoords = readsTxcoords[indNoskip, , drop = FALSE]

  ind = which(readsTxcoords[,1] > readsTxcoords[,3])
  if (length(ind) > 1){
    tp = readsTxcoords[ind, 1:2]
    readsTxcoords[ind, 1:2] = readsTxcoords[ind, 3:4]
    readsTxcoords[ind, 3:4] = tp
  }

  return(readsTxcoords)
}