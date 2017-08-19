estFun <- function(geneid, gene_models,
                    gene_exons_list, gene_name,
                    FPKM, N_gene, N_by_bam,
                    D, mu_D, sigma_D, A = 7, B = 2, rep_thre){
  if (geneid %% 100 == 0) print(paste("gene", geneid))
  bamInd = (N_gene > 0)

  if (sum(bamInd) < rep_thre){ return(0) }
  bin_list <- list()
  bamIDs <- which(bamInd)
  folder = "MSIQ_result/"
  for (id in 1:length(bamIDs)){
    load(paste(folder, "rep", bamIDs[id], "/", geneid,".RData" ,sep=""))
    bin_list[[id]] <- res
  }


  J <- gene_models[geneid, "transNum"]
  if (J==1){
    #print("one-isoform gene")
    tau_JSLIDE = 1
    tau_result <- matrix(rep(1,3), nrow=1)
    rownames(tau_result) <- strsplit(gene_models[geneid, "transName"], split=", ")[[1]]
    #estGamma = 1
    estE.D <- rep(1, D)

    transcript.01mat <- getIsoform(x=gene_models[geneid, ],
                        nExons=gene_models[geneid, "nExons"], J=J)
  }else{
    readsInfo <- generateReads(geneid, gene_models, gene_exons_list, bin_list,
                               D = length(bamIDs), mu_D = mu_D[bamIDs],
                               sigma_D = sigma_D[bamIDs])

    n.D <- readsInfo$n.D
    transcript.01mat <- readsInfo$transcript.01mat

    trueLogPrReads.jMinus1.list <- readsInfo$trueLogPrReads.list
    readsInfo.byRep.list <- readsInfo$readsInfo.byRep.list
    logPrReads.list <- readsInfo$logPrReads.list

    #print("JSLIDE-begin")
    ##########################################################
    # JSLIDE
    ##########################################################
    if (sum(n.D)>9000){
      #print("complex gene")
      maxIter=3000
      burnIn=500
    }else{
      maxIter=6000
      burnIn=2000
    }

    python.assign("maxIter", maxIter)
    python.assign("burnIn", burnIn)
    python.assign("J", J)
    #python.assign("D", D)
    python.assign("D", length(bamIDs))
    python.assign("n_D", n.D)
    python.assign("A", A)
    python.assign("B", B)
    python.exec("lambda_J = np.repeat(np.float(1)/J, J)")
    python.exec("initE_D = np.repeat(1, D)")
    #temp <- rep(0, D)
    temp <- rep(0, length(bamIDs))
    python.assign("LogPrReads_jMinus1_list", temp)
    #for (d in 1:D){
    for (d in 1:length(bamIDs)){
      mat_d <- trueLogPrReads.jMinus1.list[[d]]
      colnames(mat_d) <- NULL
      python.assign("mat_d", mat_d)
      python.assign("d", d)
      python.exec("LogPrReads_jMinus1_list[d-1] = np.matrix(mat_d)")
    }

    python.exec("result = Gibbs_Model(maxIter=maxIter, burnIn=burnIn,
                J=J, D=D, n_D = n_D,
                A=A, B=B, lambda_J=lambda_J,
                initE_D=initE_D,
                logPrReads_jMinus1_list=LogPrReads_jMinus1_list,
                seed=1
    )")


    #python.exec("a=result['estGamma'][0]")
    #estGamma <- python.get("a")
    tau_JSLIDE <- python.get("list(result['tau_JSLIDE'])")
    estE.D <- python.get("list(result['estE_D'])")

    names(tau_JSLIDE) = strsplit(gene_models[geneid, "transName"], split=", ")[[1]]
    }


  if (J==1){
    iso_FPKM = matrix(FPKM, nrow=1)
    rownames(iso_FPKM) <- gene_models$transName[geneid]
  }else{
    exon_Len <- (gene_exons_list[[gene_name]])$exon_Len
    iso_Len = as.numeric(transcript.01mat %*% matrix(exon_Len, ncol=1))
    tau_JSLIDE = round(tau_JSLIDE, digits = 4)
    iso_FPKM <- sapply(1:D, function(d){
      N_gene[d]*tau_JSLIDE*10^9 / (N_by_bam[d]*(iso_Len - mu_D[d] ))
    })
    iso_FPKM[ iso_FPKM < 0 ] = 0
    rownames(iso_FPKM) <- names(tau_JSLIDE)
    estED = rep(0, D)
    estED[bamIDs] = estE.D
  }

  Fresult <- list(gene_name = gene_name,
                  tau_result = tau_JSLIDE,
                  #estGamma = round(estGamma, digits = 4),
                  estE.D = estED,
                  FPKM = round(FPKM, digits = 4),
                  iso_FPKM = round(iso_FPKM, digits = 4))

  return(Fresult)
}
