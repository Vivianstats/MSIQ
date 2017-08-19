#' use MSIQ to quantify transcript expression from multiple samples
#'
#' @param D An integer specifying the number of RNA-seq sample supplied
#' @param gtf_path A character specifying the full path of the .gtf file
#' @param bam_path A character vector. Each element gives the full path of one .bam file.
#' @param A An integer specifying the first parameter used the Beta prior. Default is 5.
#' @param B An integer specifying the second parameter used the Beta prior. Default is 2.
#' @param ncores An integer denoting the number of cores used for parallel computation.
#' @return Save the estimation results to a text file transcript_summary.txt.
#' @export
#' @import parallel
#' @import gtools
#' @import rbamtools
#' @import rPython
#' @import GenomicFeatures
#' @import GenomicAlignments
msiq = function(D, gtf_path, bam_path, A = 5, B = 2, ncores){
  # supfolder = "~/Dropbox/MSIQ_package/MSIQ/R/"
  # python.load(paste(supfolder, "Gibbs_functions.py", sep=""))
  py_path = system.file("exec", "Gibbs_functions.py", package = "MSIQ")
  print(py_path)
  python.load(py_path)

  ###############################################
  ### get gene_models from annotation (gtf file)
  ###############################################
  folder = "MSIQ_result/"
  dir.create(folder, showWarnings = TRUE, recursive = FALSE)
  gene_path = paste(folder, "gene_models.RData", sep = "")
  gene_models = gtf_to_gene_models(gtf_path, save_path = gene_path, ncores = ncores)


  load(file=gene_path)


  ###############################################
  ### get statistics from bam file
  ###############################################


  bai_path = paste(bam_path, ".bai", sep = "")

  print("extarct basic statistics from bam files ...")

  part1 <- get_FPKM_and_num(gene_models, D, bam_path, bai_path, num_cores = ncores)
  FPKM_by_gene <- part1$FPKM_by_gene
  N_by_bam <- part1$N_by_bam
  save(FPKM_by_gene, file=paste(folder, "FPKM_by_gene.RData", sep=""))
  save(N_by_bam, file=paste(folder, "N_by_bam.RData", sep=""))

  frag_len = get_frag_len(gene_models, D, bam_path, strandmode = 0, folder, ncores = ncores)

  ###############################################
  ### process reads in bam files
  ###############################################
  print("extarct reads from bam files ...")
  for (i in 1:D){
    dir.create(paste(folder, "rep", i, sep = ""), showWarnings = TRUE, recursive = FALSE)
  }

  bin_by_gene <- get_bin(gene_models, D, bam_path, num_cores = ncores,
                      save_folder = folder, start_id = 1)

  ###############################################
  ### MSIQ estimation
  ###############################################
  print("start isoform quantification ...")

  load(paste(folder, "FPKM_by_gene.RData", sep=""))
  load(paste(folder, "num_by_gene.RData", sep=""))
  load(paste(folder, "N_by_bam.RData", sep=""))
  load(paste(folder, "mu_D.RData", sep=""))
  load(paste(folder, "sigma_D.RData", sep=""))


  gene_exons_list <- lapply(1:nrow(gene_models),function(i){
    exon_starts <-  as.numeric(strsplit(gene_models[i,"exon_starts"],split=",")[[1]])
    exon_ends <-  as.numeric(strsplit(gene_models[i,"exon_ends"],split=",")[[1]])
    exon_Len <- exon_ends - exon_starts + 1
    list(exon_starts=exon_starts, exon_ends=exon_ends, exon_Len=exon_Len)
  })
  names(gene_exons_list) <- gene_models[,"geneName"]


  rep_thre = ceiling(D/2) # genes should have reads in at least rep_thre samples

  dir.create(paste(folder, "estimation", sep = ""), showWarnings = TRUE, recursive = FALSE)
  Result <- mclapply(1:nrow(gene_models), function(geneid){
  #Result <- mclapply(geneIDs, function(geneid){
    if(geneid %% 50 == 0) print(paste("gene", geneid))
    gene_name = gene_models[geneid, "geneName"]
    FPKM = FPKM_by_gene[[gene_name]]
    N_gene  <- num_by_gene[[gene_name]]
    Fresult <- try(estFun(geneid, gene_models, gene_exons_list,
                          gene_name, FPKM, N_gene, N_by_bam, D, mu_D, sigma_D,
                          A = A, B = B, rep_thre = rep_thre),
                   silent=TRUE)
    gene_res_name <- paste(folder, "estimation/", geneid, "est.RData", sep="")
    save(Fresult, file=gene_res_name)
    return(geneid)
  }, mc.cores = ncores)

  reslist <- list.files(paste(folder, "estimation/", sep=""))
  reslist <- paste(folder, "estimation/", reslist, sep = "")

  resdata <- lapply(1:length(reslist), function(i){
    # print(i)
    load(reslist[i])
    if (class(Fresult) == "list"){
      mat =  matrix(nrow = nrow(Fresult$iso_FPKM), ncol = 3 + 2*D)
      mat[,1] = rownames(Fresult$iso_FPKM)
      mat[,2] = Fresult$gene_name
      mat[,3] = Fresult$tau_result
      mat[,4:(3+D)] = Fresult$iso_FPKM
      mat[, (4+D):(3+2*D)] = matrix(rep(Fresult$estE.D, nrow(Fresult$iso_FPKM)),
                                     nrow = nrow(Fresult$iso_FPKM), byrow = TRUE)
      #mat[1, 4+2*D] = Fresult$estGamma
      return(mat)
    }else{
      return(NULL)
    }
  })

  resdata = resdata[!sapply(resdata, is.null)]
  resdata = Reduce(rbind, resdata)

  colnames(resdata) = c("transcriptID", "geneID", "frac",
                        paste("FPKM_", 1:D, sep = ""),
                        paste("sample", 1:D, sep = ""))
                      #  "gamma")
  write.table(resdata, paste(folder, "transcript_summary.txt", sep = ""),
              quote = FALSE, row.names = FALSE)

  # subdata = resdata[, c(2, (4+D):(3+2*D))]
  # subdata = subdata[complete.cases(subdata), ]
  # write.table(subdata, paste(folder, "sample_summary.txt", sep = ""),
  #             quote = FALSE, row.names = FALSE)

  # save(resdata, file = paste(folder, "resdata_test.RData", sep = ""))

  print("calculation finished!")
  return(0)
}
