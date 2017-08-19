

### supporting function
# partialSumOperator <- function(d){
#   res = matrix(0,nrow=d, ncol=d);
#   res[lower.tri(res,diag=TRUE)]=1;
#   return(res)
# }



### extract isoform info from gene models
getIsoform <- function(x, nExons, J){
  split_temp <-  strsplit(as.character(x["isoforms"]), split="),")[[1]]
  transcript.01mat <- matrix(rep(0, J*nExons), nrow=J)
  transInd <- list()
  for (i in 1:J){
    if (i == J){
      split_temp[i] <- substr(split_temp[i], 3, nchar(split_temp[i])-2)
      split_temp[i] <- gsub("," , "", split_temp[i])
      trans <- as.numeric(strsplit(split_temp[i], split=" ")[[1]])
    }else{
      split_temp[i] <- substr(split_temp[i], 3, nchar(split_temp[i]))
      split_temp[i] <- gsub("," , "", split_temp[i])
      trans <- as.numeric(strsplit(split_temp[i], split=" ")[[1]])
    }
    transcript.01mat[i,trans] <- 1
  }
#   ind <- !duplicated(transcript.01mat, MARGIN=1)
#   transcript.01mat <- transcript.01mat[ind, , drop=FALSE]
#   J <- nrow(transcript.01mat)
   transcript_name = paste("Transcript",1:J, sep="_"); 
#   
   rownames(transcript.01mat) = transcript_name; 
   colnames(transcript.01mat) = 1:ncol(transcript.01mat)
#  return(list(transcript.01mat,J))
  return(transcript.01mat)
}




##############################################
### main function for reads generation #######
##############################################

generateReads <- function(geneid, gene_models, gene_exons_list, bin_list, D, mu_D, sigma_D){
  
  #print("pass")  
  
  gene = geneid
  x = gene_models[gene, ]
  
  exon_starts <- gene_exons_list[[gene]]$exon_starts
  exon_ends <- gene_exons_list[[gene]]$exon_ends
  exon_Len <- gene_exons_list[[gene]]$exon_Len
  names(exon_Len) = 1:length(exon_Len)
  nExons = as.numeric(x["nExons"])
  
  temp1 = strsplit(as.character(x[5]), split="),")[[1]];
  J = as.numeric(x["transNum"]) # number of isoforms
  transcript.01mat <- getIsoform(x=x, nExons=nExons, J=J)
  #transcript.01mat <- temp[[1]]
  #J <- temp[[2]]
  transcript_name <- rownames(transcript.01mat)
  
  transcript_Len = as.numeric(transcript.01mat %*% matrix(exon_Len, ncol=1)) 
  names(transcript_Len) = transcript_name;  
  
  
  lenOfExonsInBtw.list = list(); 
  for(b.exonId in 1:(nExons-2)){
    for(e.exonId in (b.exonId+2):nExons){
      # indicating whether there are exons in between
      rowIds = sapply(1:J, function(j) 
        sum(transcript.01mat[j,c(b.exonId,e.exonId)]==0)==0);
      if(sum(rowIds)>0){
        lenOfExonsInBtw.list[[paste(b.exonId,e.exonId,sep="_")]] = 
          transcript.01mat[rowIds,(b.exonId+1):(e.exonId-1), drop=FALSE]%*%matrix(exon_Len[(b.exonId+1):(e.exonId-1)], ncol=1); 
      }
    }
  }
  
  
  n.D = sapply(bin_list, function(list){
    length(list[[1]])
  })
  
  #mu = 157
  #sigma=22
  
  readsInfo.byRep.list <- list()
  for (d in 1:D){
    readsInfo.byRep.list[[paste("Rep_",d,sep="")]] = bin_list[[d]]
  }
  
  
  logPrReads.jMinus1.list = list()
  logPrReads.list = list()
  PrReads.list = list()
  for(d in 1:D){
    #print(d)
    #print(paste("d=",d,sep=""))
    obs_list = readsInfo.byRep.list[[d]]
    
    # b1, e1, b2, e2, b.coordinateOnExon, e.coordinateOnExon
    obs <- matrix(NA, ncol=6, nrow=length(obs_list[[1]]))
    colnames(obs) <- c("b1", "e1", "b2", "e2", "b.coordinateOnExon", "e.coordinateOnExon")
    obs[, "b1"] <- sapply(obs_list[[1]], function(l){ min(l[[1]]) })
    obs[, "e1"] <- sapply(obs_list[[1]], function(l){ max(l[[1]]) })
    obs[, "b2"] <- sapply(obs_list[[1]], function(l){ min(l[[2]]) })
    obs[, "e2"] <- sapply(obs_list[[1]], function(l){ max(l[[2]]) })
    obs[, c("b.coordinateOnExon", "e.coordinateOnExon")] <- obs_list[[2]]
    
    ## Step 1: compute the length of the missing part (exons between two ends) --------------------------------------------
    lenOfExonsInBtw.byRead = matrix(0, ncol=J, nrow = nrow(obs)); 
    colnames(lenOfExonsInBtw.byRead) = rownames(transcript.01mat); 
    
    for(b.exonId in 1:nExons){
      for(e.exonId in b.exonId:nExons){
        rowIds = (obs[,"e1"]==as.character(b.exonId))&
          (obs[,"b2"]==as.character(e.exonId));
        
        if( sum(rowIds) > 0 ){
          temp = lenOfExonsInBtw.list[[paste(b.exonId,e.exonId,sep="_")]]; 
          
          if(length(temp)>0){
            lenOfExonsInBtw.byRead[rowIds,rownames(temp)]=
              matrix(1,nrow=sum(rowIds),ncol=1)%*%matrix(temp, nrow=1); 
          }
        }
      }
    }
    
    
    ## step2: Compute the eact length of each fragment! Store in a new matrix -------------------------------------------
    # -Input: "obs" matrix defined as above 
    imputedFragLen.mat = matrix(NA, ncol = J, nrow = nrow(obs)); 
    colnames(imputedFragLen.mat) = transcript_name
    # Type I: b1 & e2 on the same exon ==========================================================
    rowIds.beSameExon = which(obs[,"b1"] == obs[,"e2"]); 
    if (length(rowIds.beSameExon)>0){
      len.vec = as.numeric(obs[rowIds.beSameExon,"e.coordinateOnExon"]) - 
        as.numeric(obs[rowIds.beSameExon,"b.coordinateOnExon"]); 
      exon.vec = obs[rowIds.beSameExon, "b1"]
      ifExonInTranscript = sapply(exon.vec, function(x) {transcript.01mat[,x]==1})
      ifExonInTranscript = matrix(ifExonInTranscript, nrow=J)
      rownames(ifExonInTranscript) <- rownames(transcript.01mat)
      len.mat = t(sapply(1:length(len.vec), function(x) ifExonInTranscript[,x]*len.vec[x]))
      imputedFragLen.mat[rowIds.beSameExon,colnames(len.mat)] = len.mat
    }
    
    
    # Type II: 4-dim on 2-4 exons ==========================================================
    rowIds.be2Exons = which(apply(obs[,1:4, drop=FALSE], 1, function(row){
      identical( sort(unique(row)), as.numeric(row[1]:(row[1]+1))) |
        identical( sort(unique(row)), as.numeric(row[1]:(row[1]+2))) |
        identical( sort(unique(row)), as.numeric(row[1]:(row[1]+3)))
    }))
    
    if (length(rowIds.be2Exons)>0){
      obs.typeII = obs[rowIds.be2Exons, , drop=FALSE]
      #exonPair.typeII.mat = obs.typeII[, c("b1", "e2")]; 
      exonPair.typeII.list <- lapply(1:nrow(obs.typeII), function(row){
        sort(unique(obs.typeII[row,1:4]))
      }) 
    
      ifExonPairInTranscript.typeII = sapply(exonPair.typeII.list, function(x){
        rowSums(transcript.01mat[,x, drop=FALSE]==0)==0
      })
      rownames(ifExonPairInTranscript.typeII) <- rownames(transcript.01mat)
      
      len.vec = as.numeric(obs[rowIds.beSameExon,"e.coordinateOnExon"]) - 
        as.numeric(obs[rowIds.beSameExon,"b.coordinateOnExon"]); 
      
      len.vec = sapply(1:length(rowIds.be2Exons), function(r){
        as.numeric(obs.typeII[r, "e.coordinateOnExon"]) - 
          as.numeric(obs.typeII[r,"b.coordinateOnExon"])+
          sum(exon_Len[ exonPair.typeII.list[[r]][-length(exonPair.typeII.list[[r]])] ])
      })

      
      len.mat = t(sapply(1:length(len.vec), 
                         function(x) ifExonPairInTranscript.typeII[,x]*len.vec[x]))
      imputedFragLen.mat[rowIds.be2Exons,colnames(len.mat)] = len.mat
    }
    
    
    # Type III: (b1,e1) & (b2,e2) different exons ===========================================================
    # take into consideration the missing part in between
    type_ind <- setdiff(1:nrow(obs), union(rowIds.beSameExon, rowIds.be2Exons))
    obs.typeIII = obs[type_ind, , drop=FALSE]
  
    if (length(type_ind)>0){
      #temp.list <- obs_list[[1]][-c(rowIds.beSameExon, rowIds.be2Exons)]
      temp.list <- obs_list[[1]][type_ind]
      exonPair.typeIII.list <- lapply(temp.list, function(l){
        unique(unlist(l))
      })     
      ifExonPairInTranscript.typeIII = sapply(exonPair.typeIII.list, function(x){
        rowSums(transcript.01mat[,x, drop=FALSE]==0)==0
      })
      rownames(ifExonPairInTranscript.typeIII) <- rownames(transcript.01mat)
      
      exonPair.left.list <- lapply(temp.list, function(l){
        setdiff(l[[1]], max(l[[2]]))
      })
      exonPair.right.list <- lapply(1:length(temp.list), function(r){
        setdiff(temp.list[[r]][[2]], exonPair.left.list[[r]])
      })
      
      len.Part1.vec = sapply(1:nrow( obs.typeIII), function(r){
        sum(exon_Len[exonPair.left.list[[r]]]) - obs.typeIII[r, "b.coordinateOnExon"]
      })
      len.Part2.vec = sapply(1:nrow( obs.typeIII), function(r){
        vec <- exonPair.right.list[[r]]
        if (length(vec)==1){
          return( obs.typeIII[r, "e.coordinateOnExon"] )
        }else{
          return (sum(exon_Len[ vec[-length(vec)] ])+obs.typeIII[r, "e.coordinateOnExon"] )
        }
      })
      len.Part12.vec = len.Part1.vec + len.Part2.vec; 
      len.Part12.mat = t(sapply(1:length(len.Part12.vec), 
                                function(x) ifExonPairInTranscript.typeIII[,x]*len.Part12.vec[x]))
      
      #len.Part3.mat = lenOfExonsInBtw.byRead[-c(rowIds.beSameExon, rowIds.be2Exons),];
      len.Part3.mat = lenOfExonsInBtw.byRead[type_ind,]
      #imputedFragLen.mat[-c(rowIds.beSameExon, rowIds.be2Exons),colnames(len.Part12.mat)] = len.Part12.mat + len.Part3.mat;
      imputedFragLen.mat[type_ind, colnames(len.Part12.mat)] = len.Part12.mat + len.Part3.mat;
    }
     
    
    # Impossible reads ====================================================================
    imputedFragLen.mat[imputedFragLen.mat==0]=Inf;
    
    ###########################
    effLen = transcript_Len - mu_D[d]
    effLen[ effLen <= 0 ] = 1e-5
    
    logTranscriptEffectiveLen.vec = log(effLen); 
    logTranscriptEffectiveLen.jMinus1.vec = log(effLen)- 
      logTranscriptEffectiveLen.vec[1]; 
    logTranscriptEffectiveLen.jMinus1.mat = matrix( 1,nrow=n.D[d],ncol=1 ) %*% 
      matrix( logTranscriptEffectiveLen.jMinus1.vec, nrow=1 );
    
    logPrFragLen.mat = dnorm( imputedFragLen.mat, mean=mu_D[d], sd = sigma_D[d], log=TRUE );
    
    # check !!!
    logPrFragLen.mat[logPrFragLen.mat==-Inf] = -300
       
    logPrFragLen.jMinus1.mat = logPrFragLen.mat - logPrFragLen.mat[,"Transcript_1"]%*%matrix(1,ncol=J,nrow=1)
    logPrReads.jMinus1.mat = logPrFragLen.jMinus1.mat - logTranscriptEffectiveLen.jMinus1.mat; 
    
    logPrReads.jMinus1.list[[paste("Rep",d,sep="_")]] = logPrReads.jMinus1.mat
    
       
    ###########################
    logTranscriptEffectiveLen.mat = matrix( 1,nrow=n.D[d],ncol=1 ) %*% 
      matrix( logTranscriptEffectiveLen.vec, nrow=1 );
    
    logPrFragLen.mat = dnorm( imputedFragLen.mat, mean=mu_D[d], sd = sigma_D[d], log=TRUE );
    logPrFragLen.mat[logPrFragLen.mat==-Inf] = -300
    
    logPrReads.mat = logPrFragLen.mat - logTranscriptEffectiveLen.mat; 
    
    logPrReads.list[[paste("Rep",d,sep="_")]] = logPrReads.mat;
    
    #PrReads.list[[paste("Rep",d,sep="_")]] = 10^logPrReads.mat;
  }  
  
  reads_list <- list(
    J=J,
    n.D = n.D,
    trueLogPrReads.list = logPrReads.jMinus1.list,
    logPrReads.list = logPrReads.list,
    readsInfo.byRep.list = readsInfo.byRep.list,
    # PrReads.list = PrReads.list,
    transcript.01mat = transcript.01mat)
  
  return(reads_list)
}

