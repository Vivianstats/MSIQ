# folder = "MSIQ_result/"
# dir.create(folder, showWarnings = TRUE, recursive = FALSE)
# gtf_path <- "~/Dropbox/Iso\ Discovery/Codes/gene_models/human/gencode.v24.nano.gtf"
# save_path = paste(folder, "gene_models.RData", sep = "")

gtf_to_gene_models = function(gtf_path, save_path, ncores){
  
  txdb <- makeTxDbFromGFF(gtf_path, format="gtf")
  #dropSeqlevels(txdb, "chrM")
  keepStandardChromosomes(txdb)
  
  genes_list = genes(txdb)
  gene_names <- names(genes_list)
  chrs = as.character(seqnames(genes_list))
  strands = as.character(strand(genes_list))
  exons_list_by_gene <- exonsBy(txdb, by="gene")
  exons_list_by_tx = exonsBy(txdb, by="tx", use.names = TRUE)
  txs_list_by_gene = transcriptsBy(txdb, "gene")
  
  #geneid = 1
  
  gene_models <- mclapply(1:length(gene_names), function(geneid){
    #gene_models <- lapply(1:2, function(geneid){
    exons = disjoin(exons_list_by_gene[[geneid]])
    exon_starts = start(exons)
    exon_ends = end(exons)
    exon_lens = width(exons)
    
    txs = txs_list_by_gene[[geneid]]
    txs_list = lapply(1:length(txs), function(ii){
      tx_name = txs$tx_name[ii]
      exons_in_tx = exons_list_by_tx[[tx_name]]
      index = sort(subjectHits(findOverlaps(exons_in_tx, exons)))
      unique(index)
    })
    names(txs_list) = txs$tx_name
    if (geneid %% 1000 == 0) {print(geneid); gc()}
    
    return(list(chr = chrs[geneid], exonStarts = exon_starts,
                exonEnds = exon_ends, exonLens = exon_lens,
                txs = txs_list, exonNum = length(exon_starts), 
                txNum = length(txs_list),
                str = strands[geneid]))
    
  }, mc.cores = ncores)
  
  
  
  chr = sapply(gene_models, `[[`, "chr")
  chr = sapply(strsplit(chr, "r"), `[[`, 2)
  nExons =  sapply(gene_models, `[[`, "exonNum")
  exon_starts = sapply(gene_models, `[[`, "exonStarts")
  exon_starts = sapply(exon_starts, function(x) paste(x, collapse = ", "))
  exon_ends =  sapply(gene_models, `[[`, "exonEnds")
  exon_ends = sapply(exon_ends, function(x) paste(x, collapse = ", "))
  transName = lapply(gene_models, function(x) {
    txnames = names(x$txs)
    return(sapply(strsplit(txnames, "\\."), `[[`, 1))
    })

  transName = sapply(transName, function(x) paste(x, collapse = ", "))
  isoforms = lapply(gene_models, function(x) {
    lapply(x$txs, function(xx) paste(xx, collapse = ", "))
  })
  isoforms = sapply(isoforms, function(x){
    temp = paste(unlist(x), collapse = "), (")
    paste("[(", temp, ")]")
  })
  transNum = sapply(gene_models, `[[`, "txNum")
  str = sapply(gene_models, `[[`, "str")
  gene_names = sapply(strsplit(gene_names, "\\."), `[[`, 1)
  gene_df = data.frame(geneName = gene_names, chr = chr, nExons = nExons,
                       exon_starts = exon_starts, exon_ends = exon_ends,
                       transName = transName, isoforms = isoforms,
                       transNum = transNum, str= str, stringsAsFactors = FALSE)
  gene_models = gene_df
  "get annotation from gtf, finished!"
  save(gene_models, file = save_path)
  
}





