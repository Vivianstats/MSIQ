# updated from process_reads_from_bam_20150902.R
# instead of bin_by_bam -->> bin_by_gene
# save an RData for each gene in each bam



#####################################################
#### estimate fragment length
#### from single-iso genes
#####################################################
get_frag_len = function(gene_models, D, bam_path, strandmode = 0, folder, ncores){

  rowID = which(gene_models[,"transNum"] == 1 & gene_models[,"nExons"] > 1)
  gene_exons_list <- lapply(1:nrow(gene_models),function(i){
    exon_starts <-  as.numeric(strsplit(gene_models[i,"exon_starts"],split=",")[[1]])
    exon_ends <-  as.numeric(strsplit(gene_models[i,"exon_ends"],split=",")[[1]])
    exon_Len <- exon_ends - exon_starts + 1
    list(exon_starts=exon_starts, exon_ends=exon_ends, exon_Len=exon_Len)
  })
  names(gene_exons_list) <- gene_models[,"geneName"]

  res = sapply(1:D, function(d){
    #print(paste("replicate", d))
    fglensList = mclapply(rowID, function(rowid){
      #if (rowid %% 200 == 0) print(rowid)
      cgene = gene_models[rowid, ]
      gene_exons = gene_exons_list[[rowid]]
      readTxcoords = get_reads(d, cgene, gene_exons, bam_path, strandmode = 0)

      #nrow(readTxcoords)

      if (is.null(readTxcoords)){ return(NULL)
      }else{
        return(readTxcoords[,4] - readTxcoords[,1] + 1)
      }
    }, mc.cores = ncores)

    fglensList = fglensList[!sapply(fglensList, is.null )]
    fglens = unlist(fglensList)
    fglens = fglens[fglens < quantile(fglens, 0.99)]

    if (length(fglens) > 1000){
      fgmean = mean(fglens)
      fgsd = sd(fglens)
    }else{
      fgmean = 200
      fgsd = 80
    }
    return(c(fgmean, fgsd))
  })

  mu_D = res[1, ]
  sigma_D = res[2, ]


  save(mu_D, file=paste(folder, "mu_D.RData", sep=""))
  save(sigma_D, file=paste(folder, "sigma_D.RData", sep=""))
  return(0)
}



##############################################################
######### get bin and FPKM info
##############################################################
# load(paste(folder,"gene_models_chr21.Rdata", sep=""))

get_FPKM_and_num <- function(gene_models, D, bam_path, bai_path, num_cores){
  gene_exons_list <- lapply(1:nrow(gene_models),function(i){
    exon_starts <-  as.numeric(strsplit(gene_models[i,"exon_starts"],split=",")[[1]])
    exon_ends <-  as.numeric(strsplit(gene_models[i,"exon_ends"],split=",")[[1]])
    exon_Len <- exon_ends - exon_starts + 1
    list(exon_starts=exon_starts, exon_ends=exon_ends, exon_Len=exon_Len)
  })
  names(gene_exons_list) <- gene_models[,"geneName"]


  res <- mclapply(1:D, function(d){

    print(paste("get FPKMs for replicate", d))
    bam <- bam_path[d]
    idx <- bai_path[d]

    reader <- bamReader(bam, indexname=idx, idx=TRUE)

    #loadIndex(reader, idx)
    #indexInitialized(reader)
    #getRefData(reader)

    #print("start getting FPKMs...")
    bam_count <- bamCountAll(reader, verbose=FALSE)
    N_FPKM <- sum(bam_count[,"nAligns"])
    #N_by_bam[d] <- N_FPKM

    RefData <- getRefData(reader)

    ### calculate FPKM

    FPKM_num_list <- lapply(1:nrow(gene_models), function(geneid){
      #if (geneid%%100 == 0)
      #print(paste("geneid", geneid))
      gene = gene_models[geneid, ]
      l_FPKM = sum(gene_exons_list[[geneid]]$exon_Len) # length of gene
      start = gene_exons_list[[gene[1,"geneName"]]]$exon_starts[1]
      end = gene_exons_list[[gene[1,"geneName"]]]$exon_ends[gene[1,"nExons"]]
      coords <- c(RefData[RefData$SN == paste("chr",gene[1,"chr"], sep =""), "ID"],
                  start, end)
      #coords <- c(refid[as.character(gene[1,"chr"])], start, end)
      names(coords) <- c("refid","start","stop")
      range <- bamRange(reader,coords)

      rangedf <- as.data.frame(range)
      if_pair <- rangedf[, "name"] %in% rangedf[, "name"][duplicated(rangedf[, "name"])]

      num <- 0.5 * sum(if_pair)
      FPKM <- num * 10^9 / (l_FPKM * N_FPKM)
      return(list(FPKM=FPKM, num=num))
    })

    FPKM <- lapply(FPKM_num_list, function(x) x$FPKM)
    num_reads <- lapply(FPKM_num_list, function(x) x$num)
    names(FPKM) <- gene_models[, "geneName"]
    names(num_reads) <- gene_models[, "geneName"]

    bamClose(reader)
    return(list(FPKM = FPKM,
                num = num_reads,
                N = N_FPKM))
  }, mc.cores = num_cores)

  FPKM_by_bam <- lapply(res, function(x) x$FPKM)
  N_by_bam <- sapply(res, function(x) x$N) # total number of reads in each bam
  num_by_bam <- lapply(res, function(x) x$num) # total number of reads of each gene in each bam

  FPKM_by_gene <- lapply(1:nrow(gene_models), function(r){
    sapply(FPKM_by_bam, function(x) x[[r]])
  })
  names(FPKM_by_gene) <- gene_models[ ,"geneName"]

  # num_by_gene <- lapply(1:nrow(gene_models), function(r){
  #   sapply(num_by_bam, function(x) x[[r]])
  # })
  # names(num_by_gene) <- gene_models[ ,"geneName"]

  return(list(FPKM_by_gene = FPKM_by_gene,
              #num_by_gene = num_by_gene,
              N_by_bam = N_by_bam))
}

##############################################################
######### get bin and FPKM info
##############################################################


get_bin <- function(gene_models, D, bam_path, num_cores, save_folder,
                    start_id = 1){
  #print("extract exons' information")
  gene_exons_list <- lapply(1:nrow(gene_models),function(i){
    exon_starts <-  as.numeric(strsplit(gene_models[i,"exon_starts"],split=",")[[1]])
    exon_ends <-  as.numeric(strsplit(gene_models[i,"exon_ends"],split=",")[[1]])
    exon_Len <- exon_ends - exon_starts + 1
    list(exon_starts=exon_starts, exon_ends=exon_ends, exon_Len=exon_Len)
  })
  names(gene_exons_list) <- gene_models[,"geneName"]

  num_by_bam <- mclapply(1:D, function(d){

    num <- sapply(start_id:nrow(gene_models), function(id){
      # print(paste("geneid", id))
      if (id%%1000 == 0) print(paste("geneid", id))

      cgene = gene_models[id, ]

      gene_exons = gene_exons_list[[id]]
      readTxcoords = get_reads(d, cgene, gene_exons, bam_path, strandmode = 0)


      if (is.null(readTxcoords)){
        #save(NULL, file=paste(save_folder, "rep", d, "/",id,".RData" ,sep=""))
        #num_by_gene[id, d] = 0
        return(0)
      }else{
        #num_by_gene[id, d] = nrow(readTxcoords)

        cumlen = c(0, cumsum(gene_exons$exon_Len)) + 1
        readIdx = sapply(1:4, function(k){
          findInterval(readTxcoords[, k], cumlen)
        })
        if (nrow(readTxcoords) == 1){
          readIdx  = matrix(readIdx, nrow = 1)
        }
        bin = lapply(1:nrow(readIdx), function(i){
          list(readIdx[i,1]:readIdx[i,2], readIdx[i,3]:readIdx[i,4])
        })
        coordinateOnExon = sapply(c(1,4), function(k){
          readTxcoords[,k] - (cumlen[readIdx[,k] ] - 1)
        })

        res <- list(bin=bin, coordinateOnExon=coordinateOnExon)
        save(res, file=paste(save_folder, "rep", d, "/",id,".RData" ,sep=""))
      }

      return( nrow(readTxcoords) )
    })

    gc()
    return(num)
  }, mc.cores = num_cores)

  num_by_gene = lapply(1:nrow(gene_models), function(id){
    sapply(num_by_bam, "[[", id)
  })
  names(num_by_gene) = gene_models$geneName
  save(num_by_gene, file=paste(save_folder, "num_by_gene.RData", sep=""))

  return( 0 )
}


