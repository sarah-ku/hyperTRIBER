getHits <- function(res,stranded=F,fdr=0.1,fold=2,ncore=20,addMeta=TRUE,data_list,edits_of_interest=rbind(c("A","G"),c("T","C")),design_vector=c(rep("control",5),rep("treat",5)),include_ref=T,refGR=locsGR,symm=F)
{

  #filter by edits of interest
  targs <- gsub("E(.)","\\1",res$featureID)
  res <- res[grep(paste(unique(edits_of_interest[,2]),collapse="|"),res$featureID),]
  #res.sig <- res[intersect(which(res$padj<fdr), which(res$log2fold_treat_control>fold)),]
  res.sig <- res[which(res$padj<=fdr),]
  my.hits <- unique(res.sig$groupID)

  print(paste("number of hits = ",length(my.hits)))

  if(stranded)
  {
    my.chrs <- gsub(("(.+)\\_[0-9]+,[+|-]"),"\\1",my.hits)
    my.pos <- as.numeric(as.vector(gsub((".+\\_([0-9]+),[+|-]"),"\\1",my.hits)))
    my.strand <-gsub((".+\\_[0-9]+,([+|-])"),"\\1",my.hits)
    posGR <- GRanges(Rle(my.chrs),strand=Rle(my.strand),IRanges(my.pos,my.pos))
    names(posGR) <- my.hits
  }else{
    my.chrs <- gsub(("(.+)\\_[0-9]+"),"\\1",my.hits)
    my.pos <- as.numeric(as.vector(gsub((".+\\_([0-9]+)"),"\\1",my.hits)))
    posGR <- GRanges(Rle(my.chrs),IRanges(my.pos,my.pos))
    names(posGR) <- my.hits

    if(!addMeta)
    {
      res.sig <- res.sig[order(res.sig$padj),]
      edit_base <- tapply(res.sig$featureID,res.sig$groupID,function(x) x[1])
      posGR$edit_base <- 0
      posGR[names(edit_base)]$edit_base <- gsub("E","",edit_base)
    }

  }


  if(addMeta){

    "adding meta..."

    registerDoParallel(cores=ncore)

    #meta <- list()
    meta <- foreach(i=1:length(posGR)) %dopar% {
      tmp <- do.call(rbind,lapply(data_list,function(x) x[names(posGR)[i],]))
      ref <- names(which(colSums(tmp)==max(colSums(tmp))))[1] #todo better way than choosig the first on event of tie
      hit <- res.sig[which(res.sig$groupID==names(posGR)[i]),]

      #what happens if more than one hit for same location
      if(nrow(hit)>1)
      {
        pos.edits <- matrix(edits_of_interest[edits_of_interest[,1]==ref,],ncol=2)
        pos.targs <- gsub("E(.)","\\1",hit$featureID)
        to.use <- which(pos.targs %in% pos.edits[,2])
        if(length(to.use)==1)
        {
          hit <- hit[to.use,]
        }else{
          print(paste("two edits at base",names(posGR)[i],", selecting most significant",sep=""))
          hit <- hit[to.use,]
          hit <- hit[order(hit$pvalue)[1],]
        }
      }

      targ <-  gsub("E(.)","\\1",hit$featureID)
      #refGR[]
      padj <- hit$padj
      pval <- hit$pvalue
      par_control <- hit$control
      par_treat <- hit$treat
      fc <- hit$log2fold_treat_control

      if(symm)
      {
        if(fc>0)
        {
          prop <- sum(tmp[which(design_vector=="treat"),targ]) / (sum(tmp[which(design_vector=="treat"),ref])+(sum(tmp[which(design_vector=="treat"),targ])))
          prop_control <- sum(tmp[which(design_vector=="control"),targ]) / (sum(tmp[which(design_vector=="treat"),ref])+(sum(tmp[which(design_vector=="treat"),targ])))

        }else{
          prop <- sum(tmp[which(design_vector=="control"),targ]) / (sum(tmp[which(design_vector=="control"),ref])+(sum(tmp[which(design_vector=="control"),targ])))
          prop_control <- sum(tmp[which(design_vector=="treat"),targ]) / (sum(tmp[which(design_vector=="control"),ref])+(sum(tmp[which(design_vector=="control"),targ])))
        }
      }else{
          prop <- sum(tmp[which(design_vector=="treat"),targ]) / (sum(tmp[which(design_vector=="treat"),ref])+(sum(tmp[which(design_vector=="treat"),targ])))
          prop_control <- sum(tmp[which(design_vector=="control"),targ]) / (sum(tmp[which(design_vector=="control"),ref])+(sum(tmp[which(design_vector=="control"),targ])))
      }




      edit_count <- unlist(lapply(data_list,function(x) x[names(posGR)[i],targ]))
      tags_control <- sum(edit_count[names(which(design_vector=="control"))])
      tags_treat <- sum(edit_count[names(which(design_vector=="treat"))])

      data.frame("name"=names(posGR)[i],"ref"=ref,"targ"=targ,"prop"=prop,"prop_ctrl"=prop_control,"padj"=padj,"pvalue"=pval,"control_par"=par_control,"treat_par"=par_treat,"fold_change"=fc,"tags_treat"=tags_treat,"tags_control"=tags_control)
    }
  }
  meta <- do.call(rbind,meta)
  if(include_ref)
  {
    names(refGR) <- refGR$names
    meta$base <- as.vector(refGR[paste(my.chrs,my.pos,sep="_")]$ref)
  }

  colnames(meta)[colnames(meta)=="strand"] <- "gtf_strand"
  mcols(posGR) <- meta
  posGR$ref <- as.vector(posGR$ref)
  posGR$targ <- as.vector(posGR$targ)

  posGR <- posGR[!(posGR$ref==posGR$targ)]

  return(posGR)
}



addStrandForHyperTRIBE <- function(posGR)
{
  strand <- rep("*",length(posGR))
  strand[posGR$ref=="T" & posGR$targ=="C"] <- "-"
  strand[posGR$ref=="A" & posGR$targ=="G"] <- "+"
  strand(posGR) <- Rle(strand)
  return(posGR)
}



#important since UTRs and CDS are always exons
handleAnoType <- function(x)
{
  if( ("3UTR" %in% x & "exon" %in% x) )
  {
    return("3UTR")
  }
  if( ("5UTR" %in% x & "exon" %in% x) )
  {
    return("5UTR")
  }
  if( ("CDS" %in% x & "exon" %in% x) )
  {
    return("CDS")
  }
  if(length(unique(x))==1)
  {
    return(x[1])
  }
  else{
    print("stop")
  }
}
#add in gene annotations/information. Also expression information.

###
#
#


addGenes <- function(gtfGR,posGR,ncore=30,quant=quant.vec,assignStrand=F,geneids)
{
  registerDoParallel(cores=ncore)
  #Exons = gene - introns
  #CDS = gene - introns - UTRs
  #CDS = Exons - UTRs
  feature_names <- names(table(gtfGR$type))

  #which(lapply(strsplit(posGR$genes[,1],","),function(x) length(x))==2)
  #which(lapply(strsplit(posGR$genes[,1],","),function(x) length(x))==0)

  genes <- foreach(i=1:length(posGR)) %dopar%{
    #i <- sample(1:length(posGR),1)
    ols <- gtfGR[subjectHits(findOverlaps(posGR[i],gtfGR,ignore.strand=FALSE))]
    ofr_flag <- F

    if(length(ols)==0)
    {
      ofr_flag <- T
      ols <- gtfGR[subjectHits(findOverlaps(posGR[i]+1000,gtfGR,ignore.strand=FALSE))]
    }

    genes <- unique(ols$gene_id)
    strands <- tapply(as.vector(strand(ols)),ols$gene_id,function(x) x[1])

    if(length(genes)>1)
    {
      gs <- rep(0,length(genes))
      for(q in 1:length(genes))
      {
        gs[q] <- sum(quant[grep(genes[q],row.names(tpm.mat))])
      }
      genes <- genes[rev(order(gs))]
      strands <- (strands[rev(order(gs))])

      gep <- rep(0,length(genes))
      names(gep) <- genes
      for(ge in 1:length(genes)){
        gep[ge] <- length(findOverlaps(ols[ols$gene_id==genes[ge]],posGR))
      }
      if(length(which(gep==max(gep)))==1)
      {
        genes <- names(which(gep==max(gep)))
        strands <- strands[genes]
      }else{
        genes <- genes[1]
        strands <- strands[1]
      }
      ols <- ols[(ols$gene_id %in% genes),]
    }

    #strands <- unique(strands)
    feature_type <- as.vector(ols$type)
    names(feature_type) <- ols$transcript_id
    feature_type <- tapply(feature_type,names(feature_type),handleAnoType)
    transcripts <- names(feature_type)


    if(length(strands)>0)
    {
      new_strand <- strands[1]
    }else{
      new_strand <- "*"
    }



    tquant <- transcripts[transcripts %in% names(quant)]
    #length(tquant)

    if(length(tquant)==1)
    {
      qmat <- quant[tquant]
      abund <- sum(qmat)
      names(abund) <- tquant
      transcript_expr <- rep(NA,length(transcripts))
      names(transcript_expr) <- transcripts
      transcript_expr[names(abund)] <- abund
      transcript_expr <- c(rev(sort(transcript_expr[!is.na(transcript_expr)])),transcript_expr[is.na(transcript_expr)])
      transcripts <- names(transcript_expr)
      feature_type <- feature_type[transcripts]
    }else{
      if(length(tquant)>1){
        qmat <- quant[tquant]
        abund <- ((qmat))
        transcript_expr <- rep(NA,length(transcripts))
        names(transcript_expr) <- transcripts
        transcript_expr[names(abund)] <- abund
        transcript_expr <- rev(sort(transcript_expr))
        transcripts <- names(transcript_expr)
        feature_type <- feature_type[transcripts]
      }else{
        if(length(transcripts)>0)
        {
          transcript_expr <- rep(NA,length(transcripts))
          names(transcript_expr) <- transcripts
          feature_type <- feature_type[transcripts]
        }else{
          transcripts <- ""
          transcript_expr <- ""
          feature_type <- ""
        }
      }
    }


    if(length(ols)>0 & !(new_strand=="*"))
    {
      pifv <- c()
      for(ftl in 1:length(feature_type))
      {
        fR <- ols[ols$transcript_id==names(feature_type)[ftl] & ols$type==feature_type[ftl]]
        if(new_strand=="+")
        {
          point_in_feature <- (start(posGR)[i]-start(fR))/(end(fR)-start(fR))
          names(point_in_feature) <- fR$transcript_id
        }else{
          point_in_feature <- 1-((start(posGR)[i]-start(fR))/(end(fR)-start(fR)))
          names(point_in_feature) <- fR$transcript_id
        }
        pifv <- c(pifv,point_in_feature)
        #print(point_in_feature)
      }
      names(pifv) <- names(feature_type)
      #print(pifv)
      #feature_vec <- point_in_feature

      #feature_span <- max(end(fR)) - min(start(fR))
      #names(point_in_feature) <- feature_type
      #  feature_points <- tapply(point_in_feature[names(transcript_expr[1])],feature_type[1],mean)
      #  feature_vector <- rep(NA,length(feature_names))
      #  names(feature_vector) <- feature_names
      #  feature_vector[names(feature_points)] <- feature_points


    }else{
      #feature_span <- NA
      #feature_vector <- rep(NA,length(feature_names))
      #names(feature_vector) <- feature_names
      #feature_vec <- NA
      pifv <- rep(NA,length(feature_type))
      names(pifv) <- names(feature_type)
    }

    gene.names <- geneids[as.vector(genes)]
    df<- c("genes"=paste(genes,collapse=","),"names"=paste(gene.names,collapse=","),"strand"=paste(strands,collapse=","),"transcripts"=paste(transcripts,collapse=","),"transcript.tpm"=paste(transcript_expr,collapse=","),"transcript.types"=paste(feature_type,collapse=","),"out_of_range"=ofr_flag)
    fdf <- c("fprop"=paste(round(pifv,4),collapse=","))
    #print(unlist(c(df[3],fdf)))
    list(df,fdf,new_strand)

  }
  genes.info <- do.call(rbind,lapply(genes,function(x) x[[1]]))
  feature.locations <- do.call(rbind,lapply(genes,function(x) x[[2]]))


  posGR$gene <- genes.info[,"genes"]
  posGR$name <- genes.info[,"names"]
  posGR$gtf_strand <- genes.info[,"strand"]
  posGR$transcripts <- genes.info[,"transcripts"]
  posGR$transcript.tpm <- genes.info[,"transcript.tpm"]
  posGR$transcript.types <- genes.info[,"transcript.types"]
  posGR$out_of_range <- genes.info[,"out_of_range"]


  #posGR$feature_type <- feature.locations[,"ftype"]
  posGR$feature_prop <- feature.locations[,"fprop"]



  if(assignStrand)
  {
    strands.vec <- do.call(rbind,lapply(genes,function(x) x[[3]]))
    strand(posGR) <- strands.vec
  }

  return(posGR)
}

plotExample <- function(data_list,name="",random=T,give_ref=T,refGR=locsGR)
{
  if(random){
    ab <- sample(row.names(data_list[[1]]),1)
  }else{
    ab <- name #ref is T #atcg
  }

  tmp <- do.call(rbind,lapply(data_list,function(x) x[ab,]))
  print(tmp)
  cc <- colMeans(t(apply(tmp[which(design_vector=="control"),],1,function(x) x/sum(x))))
  ct <- colMeans(t(apply(tmp[which(design_vector=="treat"),],1,function(x) x/sum(x))))
  print(round(rbind(cc/sum(cc),ct/sum(ct)),4))
  #cmat <- rbind(round(colMeans(tmp[which(design_vector=="control"),])),round(colMeans(tmp[which(design_vector=="treat"),])))
  #cmat
  #chisq.test(cmat)
  if(give_ref)
  {
    names(refGR) <- refGR$names
    print(refGR[ab])
  }
}






positionsToGenes <- function(posGR,include.ids=T,gene.ids=NULL,use.meta=T,include.desc=F,gene.descs=NULL)
{
  my.genes <- posGR$gene
  empty <- which(posGR$gene=="")
  posGR$gene[empty] <- names(posGR[empty])
  tpm_per_gene <- tapply(posGR$transcript.tpm,posGR$gene,function(x) paste(unique(unlist(strsplit(paste(x,collapse=","),","))),collapse=","))
  transcripts_per_gene <- tapply(posGR$transcripts,posGR$gene,function(x) paste(unique(unlist(strsplit(paste(x,collapse=","),","))),collapse=","))
  types_per_gene <- tapply(posGR$transcript.types,posGR$gene,function(x) paste(unique(unlist(strsplit(paste(x,collapse=","),","))),collapse=","))
  strand_per_gene <- tapply(posGR$gtf_strand,posGR$gene,function(x) paste(unique(unlist(strsplit(paste(x,collapse=","),","))),collapse=","))

  nhits <- unlist(tapply(names(posGR),posGR$gene,function(x) length(unique(x))))
  names(nhits) <- unlist(lapply(strsplit(names(nhits),","),function(x) x[1]))
  #nhits <- unlist(lapply(hit_padj,length))

  # to_handle <- which(unlist(lapply(strsplit(names(strand_per_gene),","),function(x) length(x)))>1)
  # handled <- rep("",length(to_handle))
  # for(i in 1:length(to_handle))
  # {
  #   hgenes <- strsplit(names(strand_per_gene)[to_handle[i]],",")[[1]]
  #   hlen <- rep(0,length(hgenes))
  #   for(k in 1:length(hgenes))
  #   {
  #     hlen[k] <- length(grep(hgenes[k],my.genes))
  #   }
  #   names(hlen) <- hgenes
  #   if(length(which(hlen==max(hlen)))==1)
  #   {
  #     handled[i] <- hgenes[which(hlen==max(hlen))]
  #   }else{
  #     #max_tpm
  #     nm1 <- strsplit(names(tpm_per_gene)[to_handle[i]],",")[[1]][1]
  #     handled[i] <- gsub("(.+)\\.[0-9]+","\\1",nm1)
  #   }
  # }
  #
  # new.names <- names(strand_per_gene)
  # names(new.names) <- names(strand_per_gene)
  # new.names[to_handle] <- handled
  #

  if(use.meta)
  {
    hit_ids <-  tapply(names(posGR),posGR$gene,function(x) x)
    hit_padj <-  tapply(posGR$padj,posGR$gene,function(x) x)
    hit_prop <-  tapply(posGR$prop,posGR$gene,function(x) x)
    hit_order <- tapply(posGR$padj,posGR$gene,function(x) order(x))
    hit_fcs <- tapply(posGR$fold_change,posGR$gene,function(x) x)
    hit_nature <- tapply(paste(posGR$ref,posGR$targ,sep=":"),posGR$gene,function(x) x)
    genes_order <- names(sort(tapply(posGR$padj,posGR$gene,function(x) min(x))))

    #nhits <- unlist(lapply(hit_padj,length))

    for(i in 1:length(hit_ids))
    {
      new.order <- hit_order[[i]]
      hit_padj[[i]] <- hit_padj[[i]][new.order]
      hit_prop[[i]] <- hit_prop[[i]][new.order]
      hit_ids[[i]] <- hit_ids[[i]][new.order]
      hit_nature[[i]] <- hit_nature[[i]][new.order]
      hit_fcs[[i]] <- hit_fcs[[i]][new.order]
    }

    padj_df <- unlist(lapply(hit_padj[genes_order],function(x) paste(signif(x,3),collapse=",")))
    prop_df <- unlist(lapply(hit_prop[genes_order],function(x) paste(signif(x,3),collapse=",")))
    fc_df <- unlist(lapply(hit_fcs[genes_order],function(x) paste(signif(x,3),collapse=",")))

    ids_df <- unlist(lapply(hit_ids[genes_order],function(x) paste(x,collapse=",")))
    nature_df <- unlist(lapply(hit_nature[genes_order],function(x) paste(x,collapse=",")))

    if(include.ids)
    {
      #gene.ids_order <- as.vector(gene.ids)[as.vector(genes_order)]
      my.df <- data.frame("gene"=genes_order,"name"=NA,"strand"=strand_per_gene[genes_order],
                          "transcripts"=transcripts_per_gene[genes_order],
                          "transcripts.tpm"=tpm_per_gene[genes_order],
                          "hits"=nhits[genes_order],
                          "ids"=ids_df,
                          "padj"=padj_df,
                          "props"=prop_df,
                          "fc"=fc_df,
                          "edit.type"=nature_df)
      my.df$name <- gene.ids[as.vector(my.df$gene)]
      if(include.desc)
      {
        my.df$desc <- as.vector(gene.descs[as.vector(my.df$gene)])
      }
    }else{
      my.df <- data.frame("gene"=genes_order,"strand"=strand_per_gene[genes_order],
                          "transcripts"=transcripts_per_gene[genes_order],
                          "transcripts.tpm"=tpm_per_gene[genes_order],
                          "hits"=nhits[genes_order],
                          "ids"=ids_df,
                          "padj"=padj_df,
                          "props"=prop_df,
                          "fc"=fc_df,
                          "edit.type"=nature_df)
    }
  }else{


    if(include.ids)
    {
      genes <- names(strand_per_gene)
      name <- gene.ids[genes]
      my.df <- data.frame("gene"=genes,"name"=name,"strand"=strand_per_gene,
                          "hits"=nhits,
                          "transcripts"=transcripts_per_gene,
                          "transcripts.tpm"=tpm_per_gene)
      #"hits"=nhits)
      if(include.desc)
      {
        my.df$desc <- gene.descs[as.vector(my.df$gene)]
      }
    }else{
      genes <- names(strand_per_gene)

      my.df <- data.frame("gene"=genes,"strand"=strand_per_gene,
                          "hits"=nhits,
                          "transcripts"=transcripts_per_gene,
                          "transcripts.tpm"=tpm_per_gene)
      #"hits"=nhits[genes_order])
    }
  }


  return(my.df)

}


positionsToDf <- function(posGR,include.ids,gene.ids)
{
  posGR <- posGR[order(posGR$padj)]

  #my.df <- posGR$genes
  my.df  <- data.frame("id"=names(posGR),
                       "gene"=posGR$gene,
                       #"name"=gene.ids[posGR$gene],
                       mcols(posGR)[,c("name","ref","targ","prop","padj","fold_change","tags_treat","gtf_strand","transcripts")])

  return(my.df)

}


positionsToGenes <- function(posGR,include.ids=T,gene.ids,use.meta=T)
{
  my.genes <- as.vector(posGR$gene)
  empty <- which(posGR$gene=="")
  my.genes[empty] <- names(posGR[empty])
  names(my.genes) <- names(posGR)
  posGR$gene <- my.genes
  tpm_per_gene <- tapply(posGR$transcript.tpm,posGR$gene,function(x) paste(unique(unlist(strsplit(paste(x,collapse=","),","))),collapse=","))
  transcripts_per_gene <- tapply(posGR$transcripts,posGR$gene,function(x) paste(unique(unlist(strsplit(paste(x,collapse=","),","))),collapse=","))
  types_per_gene <- tapply(posGR$transcript.types,posGR$gene,function(x) paste(unique(unlist(strsplit(paste(x,collapse=","),","))),collapse=","))
  strand_per_gene <- tapply(posGR$gtf_strand,posGR$gene,function(x) paste(unique(unlist(strsplit(paste(x,collapse=","),","))),collapse=","))

  nhits <- unlist(tapply(names(posGR),posGR$gene,function(x) length(unique(x))))
  names(nhits) <- unlist(lapply(strsplit(names(nhits),","),function(x) x[1]))
  #nhits <- unlist(lapply(hit_padj,length))

  # to_handle <- which(unlist(lapply(strsplit(names(strand_per_gene),","),function(x) length(x)))>1)
  # handled <- rep("",length(to_handle))
  # for(i in 1:length(to_handle))
  # {
  #   hgenes <- strsplit(names(strand_per_gene)[to_handle[i]],",")[[1]]
  #   hlen <- rep(0,length(hgenes))
  #   for(k in 1:length(hgenes))
  #   {
  #     hlen[k] <- length(grep(hgenes[k],my.genes))
  #   }
  #   names(hlen) <- hgenes
  #   if(length(which(hlen==max(hlen)))==1)
  #   {
  #     handled[i] <- hgenes[which(hlen==max(hlen))]
  #   }else{
  #     #max_tpm
  #     nm1 <- strsplit(names(tpm_per_gene)[to_handle[i]],",")[[1]][1]
  #     handled[i] <- gsub("(.+)\\.[0-9]+","\\1",nm1)
  #   }
  # }
  #
  # new.names <- names(strand_per_gene)
  # names(new.names) <- names(strand_per_gene)
  # new.names[to_handle] <- handled
  #

  if(use.meta)
  {
    hit_ids <-  tapply(names(posGR),posGR$gene,function(x) x)
    hit_padj <-  tapply(posGR$padj,posGR$gene,function(x) x)
    hit_prop <-  tapply(posGR$prop,posGR$gene,function(x) x)
    hit_order <- tapply(posGR$padj,posGR$gene,function(x) order(x))
    hit_fcs <- tapply(posGR$fold_change,posGR$gene,function(x) x)
    hit_nature <- tapply(paste(posGR$ref,posGR$targ,sep=":"),posGR$gene,function(x) x)
    genes_order <- names(sort(tapply(posGR$padj,posGR$gene,function(x) min(x))))

    #nhits <- unlist(lapply(hit_padj,length))

    for(i in 1:length(hit_ids))
    {
      new.order <- hit_order[[i]]
      hit_padj[[i]] <- hit_padj[[i]][new.order]
      hit_prop[[i]] <- hit_prop[[i]][new.order]
      hit_ids[[i]] <- hit_ids[[i]][new.order]
      hit_nature[[i]] <- hit_nature[[i]][new.order]
      hit_fcs[[i]] <- hit_fcs[[i]][new.order]
    }

    padj_df <- unlist(lapply(hit_padj[genes_order],function(x) paste(signif(x,3),collapse=",")))
    prop_df <- unlist(lapply(hit_prop[genes_order],function(x) paste(signif(x,3),collapse=",")))
    fc_df <- unlist(lapply(hit_fcs[genes_order],function(x) paste(signif(x,3),collapse=",")))

    ids_df <- unlist(lapply(hit_ids[genes_order],function(x) paste(x,collapse=",")))
    nature_df <- unlist(lapply(hit_nature[genes_order],function(x) paste(x,collapse=",")))

    if(include.ids)
    {
      gene.ids_order <- gene.ids[genes_order]
      my.df <- data.frame("gene"=genes_order,"name"=gene.ids_order,"strand"=strand_per_gene[genes_order],
                          "transcripts"=transcripts_per_gene[genes_order],
                          "transcripts.tpm"=tpm_per_gene[genes_order],
                          "hits"=nhits[genes_order],
                          "ids"=ids_df,
                          "padj"=padj_df,
                          "props"=prop_df,
                          "fc"=fc_df,
                          "edit.type"=nature_df)
    }else{
      my.df <- data.frame("gene"=genes_order,"strand"=strand_per_gene[genes_order],
                          "transcripts"=transcripts_per_gene[genes_order],
                          "transcripts.tpm"=tpm_per_gene[genes_order],
                          "hits"=nhits[genes_order],
                          "ids"=ids_df,
                          "padj"=padj_df,
                          "props"=prop_df,
                          "fc"=fc_df,
                          "edit.type"=nature_df)
    }
  }else{


    if(include.ids)
    {
      genes <- names(strand_per_gene)
      name <- gene.ids[genes]
      my.df <- data.frame("gene"=genes,"name"=name,"strand"=strand_per_gene,
                          "hits"=nhits,
                          "transcripts"=transcripts_per_gene,
                          "transcripts.tpm"=tpm_per_gene)
      #"hits"=nhits)
    }else{
      genes <- names(strand_per_gene)

      my.df <- data.frame("gene"=genes,"strand"=strand_per_gene,
                          "hits"=nhits,
                          "transcripts"=transcripts_per_gene,
                          "transcripts.tpm"=tpm_per_gene)
      #"hits"=nhits[genes_order])
    }
  }


  return(my.df)

}


positionsToDf <- function(posGR,include.ids,gene.ids)
{
  posGR <- posGR[order(posGR$meta$padj)]

  my.df <- posGR$genes
  my.df  <- data.frame("id"=names(posGR),
                       "gene"=my.df[,"genes"],
                       "name"=gnames[my.df[,"genes"]],
                       my.df[,-which(colnames(my.df)=="genes")],
                       posGR$meta[,-which(colnames(posGR$meta)=="name")])

  return(my.df)

}

positionsToDf <- function(posGR,include.ids,gene.ids)
{
  posGR <- posGR[order(posGR$padj)]

  meta.dat <- mcols(posGR)
  meta.dat <- (meta.dat[,-which(colnames(meta.dat) %in% c("name","genes","span","5UTR","exon","start_codon","CDS","stop_codon","3UTR"))])
  #my.df <- posGR$genes
  my.df  <- data.frame("id"=names(posGR),
                       "gene"=posGR$genes,
                       "name"=gnames[posGR$genes],
                       meta.dat)

  return(my.df)

}



getExpressedGenes <- function(cut=0.5,nreps=1,expr=tpm.mat.collapsed,samp_names=names(design_vector[design_vector=="control"]))
{
  expressed_genes <- names(which((apply(expr[,samp_names],1,function(x) length(x[x>cut])>nreps))))
  return(expressed_genes)
}
