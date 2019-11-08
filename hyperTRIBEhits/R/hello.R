#'' library(reshape2)
#' library(foreach)
#' library(doParallel)
#' suppressPackageStartupMessages( library( "DEXSeq" ) )
#'

#' data_list <- data.list
#' design_vector=c(rep("control",5),rep("treat",5))
#' names(design_vector) <- names(data_list)
#' names_vec=row.names(data_list[[1]])
#' min_samp=3
#' min_count=1
#' edits_of_interest=rbind(c("A","G"),c("T","C"))
#'


extractCountData <- function(dat,samp.names)
{
  nsamps <- length(samp.names)
  data.list <- list()
  to_rm <- list()
  for(i in 1:nsamps)
  {
    print(paste(3+c((i*4-3):(i*4))))
    data.samp <- (as.matrix(dat[,3+c((i*4-3):(i*4))]))
    colnames(data.samp) <- c("A","T","C","G")
    row.names(data.samp) <- paste(dat[,1],dat[,2],sep="_")

    suppressWarnings(class(data_list[[i]]) <- "numeric")
    to_rm[[i]] <- which(is.na(data_list[[i]]),arr.ind=T)[,1]

    data.list[[paste(samp.names[[i]])]] <- data.samp
  }

  to.rm <- unique(unlist(to_rm))
  if(length(to.rm)>0)
  {
    data.list <- lapply(data.list,function(x) x[-to.rm,])
  }

  return(data.list)
}

pca_samples <- function(data_list,design_vector)
{
  mypca <- princomp(cont_counts[["T"]][apply(cont_counts[["T"]],1,function(x) length(x[is.na(x)]))==0,])
  plot(mypca$loadings[,1],mypca$loadings[,2])
}

plot_samples <- function(data_list,design_vector)
{
  cont_counts[["A"]] <- do.call(cbind,lapply(data_list[design_vector=="control"],function(x) x$A))
  cont_counts[["T"]]<- do.call(cbind,lapply(data_list[design_vector=="control"],function(x) x$T))
  cont_counts[["C"]] <- do.call(cbind,lapply(data_list[design_vector=="control"],function(x) x$C))
  cont_counts[["G"]] <- do.call(cbind,lapply(data_list[design_vector=="control"],function(x) x$G))

  treat_counts <- list()
  treat_counts[["A"]] <- do.call(cbind,lapply(data_list[design_vector=="treat"],function(x) x$A))
  treat_counts[["T"]] <- do.call(cbind,lapply(data_list[design_vector=="treat"],function(x) x$T))
  treat_counts[["C"]] <- do.call(cbind,lapply(data_list[design_vector=="treat"],function(x) x$C))
  treat_counts[["G"]] <- do.call(cbind,lapply(data_list[design_vector=="treat"],function(x) x$G))


  hist(apply(treat_counts[["G"]],1,sum)/apply(cont_counts[["A"]],1,sum)+apply(treat_counts[["G"]],1,sum))

}

restrict_data <- function(data_list,design_vector,min_samp=3,min_count=0,edits_of_interest=rbind(c("A","G"),c("T","C")))
{

  cont_counts <- list()
  cont_counts[["A"]] <- do.call(cbind,lapply(data_list[design_vector=="control"],function(x) x[,"A"]))
  cont_counts[["T"]]<- do.call(cbind,lapply(data_list[design_vector=="control"],function(x) x[,"T"]))
  cont_counts[["C"]] <- do.call(cbind,lapply(data_list[design_vector=="control"],function(x) x[,"C"]))
  cont_counts[["G"]] <- do.call(cbind,lapply(data_list[design_vector=="control"],function(x) x[,"G"]))

  treat_counts <- list()
  treat_counts[["A"]] <- do.call(cbind,lapply(data_list[design_vector=="treat"],function(x) x[,"A"]))
  treat_counts[["T"]] <- do.call(cbind,lapply(data_list[design_vector=="treat"],function(x) x[,"T"]))
  treat_counts[["C"]] <- do.call(cbind,lapply(data_list[design_vector=="treat"],function(x) x[,"C"]))
  treat_counts[["G"]] <- do.call(cbind,lapply(data_list[design_vector=="treat"],function(x) x[,"G"]))

  depth <- rowSums(do.call(cbind,lapply(data_list[design_vector=="control"],rowSums)))

  tokeep_list <- list()
  for(i in 1:nrow(edits_of_interest))
  {
    my_ref <- edits_of_interest[i,1]
    my_targ <- edits_of_interest[i,2]
    is_max <- rowSums(cont_counts[[my_ref]])/depth > 0.5
    tokeep_treat <- rowSums(treat_counts[[my_targ]]>=min_count)>=min_samp #'need edit to appear in treatment
    tokeep_control <- rowSums(cont_counts[[my_ref]]>=min_count)>=min_samp
    tokeep_list[[i]] <- tokeep_treat & tokeep_control & is_max
    #'tokeep_list[[i]] <- tokeep_treat & is_max

  }

  to_keep <- rowSums(do.call(cbind,tokeep_list))>0

  data_list <- lapply(data_list,function(x) x[which(to_keep),])
  return(data_list)
}


generateCountFiles <- function(data_list,names_vec=row.names(data_list[[1]]),design_vector=c(rep("control",5),rep("treat",5)),out_dir="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/model_tmp_files/")
{
  print("getting ccounts")
  counts <- list()

  for(i in 1:length(data_list))
  {
    ccounts <- data_list[[i]]
    my.names <- row.names(data_list[[i]])
    #'colnames(ccounts) <- c("A","T","C","G")
    #'row.names(ccounts) <- common.names
    edit.rest <- ccounts
    edit.rest.melt <- (melt(t(edit.rest)))
    counts.c1 <- data.frame("value"=edit.rest.melt[,3])
    row.names(counts.c1) <- paste(rep(my.names,each=4),(factor(edit.rest.melt$Var1)),sep=":")
    counts[[i]] <- counts.c1
  }

  names(counts) <- names(data_list)

  #'allcounts.treat <- do.call(cbind,counts[which(design_vector=="treat")])
  #'colnames(allcounts.treat) <- paste("treat",1:ncol(allcounts.treat),sep="_")
  #'allcounts.control <- do.call(cbind,counts[which(design_vector=="control")])
  #'colnames(allcounts.control) <- paste("control",1:ncol(allcounts.control),sep="_")

  #'counts.treat <- lapply(allcounts.treat,function(x) data.frame(x,row.names = names_vec))
  #'lapply(counts.treat,nrow)
  #'counts.control <- lapply(counts.control,function(x) data.frame(x[to.use.new,],row.names = to.use.new))
  #'lapply(counts.control,nrow)

  print("saving ccounts")

  for(i in 1:length(which(design_vector=="treat")))
  {
    k <- which(design_vector=="treat")[i]
    counts.c1 <- counts[[k]]
    write.table(counts.c1,paste0(out_dir,"/countsTreat",i,".txt"),col.names=F,row.names = T,quote=F,sep="\t")
  }

  for(i in 1:length(which(design_vector=="control")))
  {
    k <- which(design_vector=="control")[i]
    counts.c1 <- counts[[k]]
    write.table(counts.c1,paste0(out_dir,"/countsControl",i,".txt"),col.names=F,row.names = T,quote=F,sep="\t")
  }



  print("writing annotation")

  mychr <- gsub("(.*)_[0-9]+:[A|T|C|G]","\\1",row.names(counts.c1))
  head(mychr)
  mypos <- as.vector((gsub(".*_([0-9]+):[A|T|C|G]","\\1",row.names(counts.c1))))
  head(mypos)
  myexon <- as.vector((gsub(".*:([A|T|C|G])","\\1",row.names(counts.c1))))
  head(myexon)

  #'print(table(tapply(myexon,paste(mychr,mypos),length)))

  myname <- as.vector(row.names(counts.c1))
  #'head(myname)
  myv3 <- "exonic_part"
  myv4 <- mypos
  myv5 <- mypos
  myv6 <- rep(".",length(row.names(counts.c1)))
  myv7 <- rep("*",length(row.names(counts.c1)))
  myv8 <- rep(".",length(row.names(counts.c1)))
  myv9 <- paste0("transcripts ",row.names(counts.c1),"; ","exonic_part_number ",myexon,"; ","gene_id ",mychr,"_",mypos)
  #'head(myv9)
  newagg <- data.frame(mychr,myname,myv3,myv4,myv5,myv6,myv7,myv8,myv9)
  #'head(newagg)
  write.table(newagg,paste(out_dir,"/annotation_file_made_up_for_DEXSeq.gff",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

  print(paste("successfully written counts and annotations to",out_dir))

}


make_test <- function(out_dir,design_vector=c(rep("control",5),rep("treat",5)),ncores=1)
{
  inDir <- out_dir
  countFiles = list.files(inDir, pattern="counts*", full.names=TRUE)
  #'countFiles <- countFiles #'remove rep 3 in ECT2 roots
  #'print(basename(countFiles))
  flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
  print(basename(flattenedFile))

  sampleTable = data.frame(
    row.names = c(paste("control",1:length(which(design_vector=="control")),sep="_"),paste("treat",1:length(which(design_vector=="treat")),sep="_")),
    condition =  design_vector
  )

  #'#' ----displaySampleTable----------------------------------------------------
  #'sampleTable

  print("running model")

  #'#' ----makeecs, eval=TRUE----------------------------------------------------

  if(exists("dxd")){
    rm(dxd)
  }

  print("create dxd")

  dxd = DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData=sampleTable,
    design= ~ sample + exon + condition:exon,
    flattenedfile=flattenedFile )

  #'#' ----sizeFactors1----------------------------------------------------------
  print("get size factors")

  dxd = estimateSizeFactors( dxd )

  BPPARAM = MulticoreParam(workers=ncores)
  #'#' ----estDisp1--------------------------------------------------------------

  print("estimate dispersion")
  dxd = estimateDispersions( dxd ,fitType="local",BPPARAM=BPPARAM)

  #'#' ----testForDEU1,cache=TRUE------------------------------------------------
  print("test for DEU")
  dxd = testForDEU( dxd ,BPPARAM=BPPARAM)

  #'#' ----estFC,cache=TRUE------------------------------------------------------
  print("estimate fold changes")
  dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition",BPPARAM=BPPARAM)

  print("returning output")

  res <- DEXSeqResults(my_dxd)

  print("saving to out directory")
  save(res,file=paste0(outdir,"/dxd_results.Rdat"))
  return(res)
}


#'pull results and make GRAnge and annotate

#' getHits <- function()
#' {
#'   #'filter by edits of interest
#'   targs <- gsub("E(.)","\\1",res$featureID)
#'   res <- res[grep(paste(unique(edits_of_interest[,2]),collapse="|"),res$featureID),]
#'   res.sig <- res[intersect(which(res$padj<fdr), which(res$log2fold_treat_control>fold)),]
#'   my.hits <- unique(res.sig$groupID)
#'
#'   print(paste("number of hits = ",length(my.hits)))
#'
#'   my.chrs <- gsub(("(.+)\\_[0-9]+"),"\\1",my.hits)
#'   my.pos <- as.numeric(as.vector(gsub((".+\\_([0-9]+)"),"\\1",my.hits)))
#'   posGR <- GRanges(Rle(my.chrs),IRanges(my.pos,my.pos))
#'   names(posGR) <- my.hits
#' }

getHits <- function(res,fdr=0.01,fold=2,ncore=20,addMeta=TRUE,data_list,edits_of_interest=rbind(c("A","G"),c("T","C")),design_vector=c(rep("control",5),rep("treat",5)))
{

  #'filter by edits of interest
  targs <- gsub("E(.)","\\1",res$featureID)
  res <- res[grep(paste(unique(edits_of_interest[,2]),collapse="|"),res$featureID),]
  res.sig <- res[intersect(which(res$padj<fdr), which(res$log2fold_treat_control>fold)),]
  my.hits <- unique(res.sig$groupID)

  print(paste("number of hits = ",length(my.hits)))

  my.chrs <- gsub(("(.+)\\_[0-9]+"),"\\1",my.hits)
  my.pos <- as.numeric(as.vector(gsub((".+\\_([0-9]+)"),"\\1",my.hits)))
  posGR <- GRanges(Rle(my.chrs),IRanges(my.pos,my.pos))
  names(posGR) <- my.hits

  if(addMeta){

    "adding meta..."

    registerDoParallel(cores=ncore)

    #'meta <- list()
    meta <- foreach(i=1:length(posGR)) %dopar% {
      tmp <- do.call(rbind,lapply(data_list,function(x) x[names(posGR)[i],]))
      ref <- names(which(colSums(tmp)==max(colSums(tmp))))
      hit <- res.sig[which(res.sig$groupID==names(posGR)[i]),]

      #'what happens if more than one hit for same location
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

      prop <- sum(tmp[which(design_vector=="treat"),targ]) / (sum(tmp[which(design_vector=="treat"),ref])+(sum(tmp[which(design_vector=="treat"),targ])))
      pval <- hit$padj
      par_control <- hit$control
      par_treat <- hit$treat
      fc <- hit$log2fold_treat_control

      data.frame("name"=names(posGR)[i],"ref"=ref,"targ"=targ,"prop"=prop,"padj"=pval,"control_par"=par_control,"treat_par"=par_treat,"fold_change"=fc)
    }
  }
  meta <- do.call(rbind,meta)
  posGR$meta <- meta
  return(posGR)
}



addStrandForHyperTRIBE <- function(posGR)
{
  strand <- rep("*",length(posGR))
  strand[posGR$meta$ref=="T" & posGR$meta$targ=="C"] <- "-"
  strand[posGR$meta$ref=="A" & posGR$meta$targ=="G"] <- "+"
  strand(posGR) <- Rle(strand)
  return(posGR)
}



#'add in gene annotations/information. Also expression information.
addGenes <- function(gtfGR,posGR,ncore=30)
{
  #'Exons = gene - introns
  #'CDS = gene - introns - UTRs
  #'CDS = Exons - UTRs
  feature_names <- names(table(gtfGR$type))

  genes <- foreach(i=1:length(posGR)) %dopar%{
    i <- sample(1:length(posGR),1)
    ols <- gtfGR[subjectHits(findOverlaps(posGR[i],gtfGR,ignore.strand=FALSE))]
    ols
    ols <- gtfGR[subjectHits(findOverlaps(posGR[i],gtfGR,ignore.strand=FALSE))]

    if(as.vector(strand(posGR[i]))=="+")
    {
      point_in_feature <- (start(posGR)[i]-start(ols))/(end(ols)-start(ols))
    }else{
      point_in_feature <- 1-((start(posGR)[i]-start(ols))/(end(ols)-start(ols)))
    }
    feature_type <- as.vector(ols$type)
    transcripts <- unique(ols$transcript_id)
    genes <- unique(ols$gene_id)

    feature_span <- max(end(ols)) - min(start(ols))

    #'names(point_in_feature) <- feature_type
    feature_points <- tapply(point_in_feature,feature_type,mean)
    feature_vector <- rep(NA,length(feature_names))
    names(feature_vector) <- feature_names
    feature_vector[names(feature_points)] <- feature_points

    df<- c("genes"=paste(genes,collapse=","),"transcripts"=paste(transcripts,collapse=","))
    df$span <- feature_span
    list(df,round(feature_vector,4))

  }
  genes.info <- do.call(rbind,lapply(genes,function(x) x[[1]]))
  feature.locations <- do.call(rbind,lapply(genes,function(x) x[[2]]))

  posGR$genes <- genes.info
  posGR$feature_locs <- feature.locations
  return(posGR)
}


plotGeneBody <- function(myGR,myGR_list=NULL,list_style=F,my.alpha=.25,my.adjust=1)
{

  if(list_style)
  {
    tp.list <- list()
    for(k in 1:length(myGR_list))
    {
      myGR <- myGR_list[[k]]
      pos_3UTR <- which(!is.na(myGR$feature_locs[,"3UTR"]))
      pos_5UTR <-which(!is.na(myGR$feature_locs[,"5UTR"]))
      pos_CDS <-which(!is.na(myGR$feature_locs[,"CDS"]))

      vals_3UTR <- myGR$feature_locs[pos_3UTR,"3UTR"]
      vals_5UTR <- myGR$feature_locs[pos_5UTR,"5UTR"]
      vals_CDS <- myGR$feature_locs[pos_CDS,"CDS"]

      to.add <- c(vals_5UTR,vals_CDS+1,vals_3UTR+2)
      tp.list[[names(myGR_list)[k]]] <- data.frame("variable"=rep(names(myGR_list)[k],length(to.add)),"value"=to.add)
    }

    df <- do.call(rbind,tp.list)
    colnames(df) <- c("variable","value")

    PC <- ggplot(df, aes(value, fill=variable)) +
      geom_density(alpha=my.alpha,adjust=my.adjust) +
      geom_vline(xintercept=1,lty=2) +
      geom_vline(xintercept=2,lty=2) +
      theme_minimal()

    PC
  }else{
    pos_3UTR <- which(!is.na(myGR$feature_locs[,"3UTR"]))
    pos_5UTR <-which(!is.na(myGR$feature_locs[,"5UTR"]))
    pos_CDS <-which(!is.na(myGR$feature_locs[,"CDS"]))

    #'remove
    #'pos_3UTR <- pos_3UTR[!(pos_3UTR %in% c(pos_CDS,pos_5UTR))]


    vals_3UTR <- myGR$feature_locs[pos_3UTR,"3UTR"]
    vals_5UTR <- myGR$feature_locs[pos_5UTR,"5UTR"]
    vals_CDS <- myGR$feature_locs[pos_CDS,"CDS"]

    tp.dat <- c(vals_5UTR,vals_CDS+1,vals_3UTR+2)

    df <- melt(data.frame(tp.dat))

    PC <- ggplot(df, aes(value, fill=variable)) +
      geom_density(alpha=my.alpha,adjust=my.adjust) +
      geom_vline(xintercept=1,lty=2) +
      geom_vline(xintercept=2,lty=2) +
      theme_minimal()

    PC
  }
}

plotEdits <- function(posGR)
{
  my.df <- table(paste(posGR$meta$ref,posGR$meta$targ,sep=":"))
  barplot(my.df)
}


genelevelStatistics <- function(posGR)
{
  #'hits per gene

  #'hits per expression level

  #'editing proportions per gene
}

peakRelativeDistance <- function(ref=posGR,targ=iCLIP[[3]])
{

}

#' makePCA <- function(data_list,my_names)
#' {
#'   data <- do.call(cbind,lapply(data_list,function(x) x[my_names,]))
#'
#'   prcomp()
#'
#'   load(paste0(server,file_path,"tpm_values_araport.Rdata"))
#'   tpm_values_subset<-tpm_values[which(rowSums(tpm_values)>0),]
#'   samples_info<-read.xlsx(paste0(server,file_path,"samples_design_ect3_quant.xlsx"))
#'
#'   pca_tmp_val<-t(tpm_values_subset)#' pca_tmp_val[,33670:33679]
#'   pca_ab_matrix<-prcomp(pca_tmp_val,center=T,scale.=T)
#'   pca_ab_matrix_summary<-summary(pca_ab_matrix)
#'   pc_df<-full_join(tibble::rownames_to_column(as.data.frame(pca_ab_matrix$x)),
#'                    samples_infopca_plot_quantification<-ggplot(data=pc_df,aes(x=PC1,y=PC2,color=factor(type),label=geom_point(size=2)+labs(color="Type")+geom_text_repel(size=5,fontface='bold')+#'geom_text(aes(label=sample), size=2, hjust=0.9, vjust=1.5)+type_colors+xlab(paste("PC1",round(pca_ab_matrix_summary$importance[2,1],4)*100,"% of variance"ylab(paste("PC2",round(pca_ab_matrix_summary$importance[2,2],4)*100,"% of variance"ggtitle('PC1 and PC2 - log2(TPM values+1)')+facet_wrap(~tissue,scales="free_x")+theme(plot.title=element_text(hjust=0.5),text=element_text(
#' }

