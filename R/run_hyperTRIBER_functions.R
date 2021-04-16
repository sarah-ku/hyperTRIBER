suppressPackageStartupMessages( library( "DEXSeq" ) )
suppressPackageStartupMessages( library( "reshape2" ) )
suppressPackageStartupMessages( library( "doParallel" ) )
suppressPackageStartupMessages( library( "GenomicRanges" ) )
suppressPackageStartupMessages( library( "rtracklayer" ) )

#########################################################
#########################################################
####EXTRACT counts from the mpileup output perl script###
#########################################################
#########################################################

extractCountData <- function(dat,samp.names,stranded=F)
{
  nsamps <- length(samp.names)
  data.list <- list()
  to_rm <- list()
  if(stranded)
  {
    for(i in 1:nsamps)
    {
      print(paste(3+c((i*8-7):(i*8))))
      data.samp <- (as.matrix(dat[,3+c((i*8-7):(i*8))]))
      colnames(data.samp) <- c("A","T","C","G","a","t","c","g")
      print(c(sum(data.samp[,1:4],na.rm=T)/sum(data.samp[,5:8],na.rm=T)))
      #we swap a<->t and c<->g to ensure the correct base name on the negative strand:
      data.samp <- data.frame(rbind(data.samp[,1:4],data.samp[,c("t","a","g","c")]))

      row.names(data.samp) <- paste(rep(dat[,1],times=2),"_",rep(dat[,2],times=2),rep(c(",+",",-"),each=nrow(dat)),sep="")

      data.list[[paste(samp.names[[i]])]] <- data.samp

      #suppressWarnings(class(data_list[[i]]) <- "numeric")
      to_rm[[i]] <- which(is.na(data.list[[i]]),arr.ind=T)[,1]
    }
  }else{
    for(i in 1:nsamps)
    {
      print(paste(3+c((i*4-3):(i*4))))
      data.samp <- (as.matrix(dat[,3+c((i*4-3):(i*4))]))
      colnames(data.samp) <- c("A","T","C","G")
      row.names(data.samp) <- paste(dat[,1],dat[,2],sep="_")

      data.list[[paste(samp.names[[i]])]] <- data.samp

      #suppressWarnings(class(data_list[[i]]) <- "numeric")
      to_rm[[i]] <- which(is.na(data.list[[i]]),arr.ind=T)[,1]
    }
  }

  to.rm <- unique(unlist(to_rm))
  if(length(to.rm)>0)
  {
    data.list <- lapply(data.list,function(x) x[-to.rm,])
  }

  return(data.list)
}

###############################################################################
###############################################################################
#RESTRICT data according to minimum number of replicates in treatment/control##
###############################################################################
###############################################################################

restrict_data <- function(data_list,design_vector,min_samp_control=3,min_samp_treat=3,min_count=0,edits_of_interest=rbind(c("A","G"),c("T","C")))
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
  #head(depth)

  tokeep_list <- list()
  for(i in 1:nrow(edits_of_interest))
  {
    my_ref <- edits_of_interest[i,1]
    my_targ <- edits_of_interest[i,2]
    is_max <- rowSums(cont_counts[[my_ref]])/depth > 0.5
    tokeep_treat <- rowSums(treat_counts[[my_targ]]>=min_count)>=min_samp_treat #need edit to appear in treatment
    tokeep_control <- rowSums(cont_counts[[my_ref]]>=min_count)>=min_samp_control
    #tokeep_list[[i]] <- tokeep_treat & tokeep_control & is_max
    tokeep_list[[i]] <- tokeep_treat & tokeep_control
    #tokeep_list[[i]] <- tokeep_treat & is_max

  }

  to_keep <- rowSums(do.call(cbind,tokeep_list))>0

  data_list <- lapply(data_list,function(x) x[which(to_keep),])
  return(data_list)
}


###############################################################################
###############################################################################
####generate files which can be opened in DEXseq###############################
###############################################################################
###############################################################################

generateCountFiles <- function(data_list,stranded=F,names_vec=row.names(data_list[[1]]),design_vector=c(rep("control",5),rep("treat",5)),out_dir="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/model_tmp_files/")
{
  print("getting ccounts")
  counts <- list()

  for(i in 1:length(data_list))
  {
    ccounts <- data_list[[i]]
    my.names <- row.names(data_list[[i]])
    #colnames(ccounts) <- c("A","T","C","G")
    #row.names(ccounts) <- common.names
    edit.rest <- ccounts
    edit.rest.melt <- (melt(t(edit.rest)))
    counts.c1 <- data.frame("value"=edit.rest.melt[,3])
    row.names(counts.c1) <- paste(rep(my.names,each=4),(factor(edit.rest.melt$Var1)),sep=":")
    counts[[i]] <- counts.c1
  }

  names(counts) <- names(data_list)


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

  if(stranded)
  {
    mychr <- gsub("(.*)_[0-9]+,[+|-]:[A|T|C|G]","\\1",row.names(counts.c1))
    head(mychr)
    mypos <- as.vector((gsub(".*_([0-9]+),[+|-]:[A|T|C|G]","\\1",row.names(counts.c1))))
    head(mypos)
    mystrand <- as.vector((gsub(".*_[0-9]+,([+|-]):[A|T|C|G]","\\1",row.names(counts.c1))))
    head(mystrand)
    myexon <- as.vector((gsub(".*:([A|T|C|G])","\\1",row.names(counts.c1))))
    head(myexon)
  }else{
    mychr <- gsub("(.*)_[0-9]+:[A|T|C|G]","\\1",row.names(counts.c1))
    head(mychr)
    mypos <- as.vector((gsub(".*_([0-9]+):[A|T|C|G]","\\1",row.names(counts.c1))))
    head(mypos)
    myexon <- as.vector((gsub(".*:([A|T|C|G])","\\1",row.names(counts.c1))))
    head(myexon)
  }


  #print(table(tapply(myexon,paste(mychr,mypos),length)))

  myname <- as.vector(row.names(counts.c1))
  #head(myname)
  myv3 <- "exonic_part"
  myv4 <- mypos
  myv5 <- mypos
  myv6 <- rep(".",length(row.names(counts.c1)))
  myv7 <- rep("*",length(row.names(counts.c1)))
  if(stranded)
  {
    myv7 <- mystrand
  }
  myv8 <- rep(".",length(row.names(counts.c1)))
  myv9 <- paste0("transcripts ",row.names(counts.c1),"; ","exonic_part_number ",myexon,"; ","gene_id ",mychr,"_",mypos)
  if(stranded)
  {
    myv9 <- paste0("transcripts ",row.names(counts.c1),"; ","exonic_part_number ",myexon,"; ","gene_id ",mychr,"_",mypos,",",mystrand)
  }
  #head(myv9)
  newagg <- data.frame(mychr,myname,myv3,myv4,myv5,myv6,myv7,myv8,myv9)
  #head(newagg)
  write.table(newagg,paste(out_dir,"/annotation_file_made_up_for_DEXSeq.gff",sep=""),col.names=F,row.names=F,quote=F,sep="\t")

  print(paste("successfully written counts and annotations to",out_dir))

}


###############################################################################
###############################################################################
####MAKE TEST (BASIC)##########################################################
###############################################################################
###############################################################################


make_test <- function(out_dir,design_vector=c(rep("control",5),rep("treat",5)),ncores=1)
{
  inDir <- out_dir
  countFiles = list.files(inDir, pattern="counts*", full.names=TRUE)
  #countFiles <- countFiles #remove rep 3 in ECT2 roots
  #print(basename(countFiles))
  flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
  print(basename(flattenedFile))

  sampleTable = data.frame(
    row.names = c(paste("control",1:length(which(design_vector=="control")),sep="_"),paste("treat",1:length(which(design_vector=="treat")),sep="_")),
    condition =  design_vector
  )

  ## ----displaySampleTable----------------------------------------------------
  #sampleTable

  print("running model")

  ## ----makeecs, eval=TRUE----------------------------------------------------

  if(exists("dxd")){
    rm(dxd)
  }

  print("create dxd")

  dxd = DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData=sampleTable,
    design= ~ sample + exon + condition:exon,
    flattenedfile=flattenedFile )

  ## ----sizeFactors1----------------------------------------------------------
  print("get size factors")

  dxd = estimateSizeFactors( dxd )

  BPPARAM = MulticoreParam(workers=ncores)
  ## ----estDisp1--------------------------------------------------------------

  print("estimate dispersion")
  dxd = estimateDispersions( dxd ,fitType="local",BPPARAM=BPPARAM)

  ## ----testForDEU1,cache=TRUE------------------------------------------------
  print("test for DEU")
  dxd = testForDEU( dxd ,BPPARAM=BPPARAM)

  ## ----estFC,cache=TRUE------------------------------------------------------
  print("estimate fold changes")
  dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition",BPPARAM=BPPARAM)

  print("returning output")

  res <- DEXSeqResults(dxd)

  print("saving to out directory")
  save(res,file=paste0(out_dir,"/dxd_results.Rdat"))
  return(res)
}



###############################################################################
###############################################################################
####MAKE TEST ACCOUNTING FOR ADAR##############################################
###############################################################################
###############################################################################

make_test_ADAR <- function(out_dir,design_vector=c(rep("control",5),rep("treat",5)),ncores=1,adar.vec=rep(1,10))
{
  inDir <- out_dir
  countFiles = list.files(inDir, pattern="counts*", full.names=TRUE)
  #countFiles <- countFiles #remove rep 3 in ECT2 roots
  #print(basename(countFiles))
  flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
  print(basename(flattenedFile))

  sampleTable = data.frame(
    row.names = c(paste("control",1:length(which(design_vector=="control")),sep="_"),paste("treat",1:length(which(design_vector=="treat")),sep="_")),
    condition =  design_vector,
    adar = as.numeric(adar.vec)
  )

  ## ----displaySampleTable----------------------------------------------------
  #sampleTable

  print("running model")

  ## ----makeecs, eval=TRUE----------------------------------------------------

  if(exists("dxd")){
    rm(dxd)
  }

  print("create dxd")

  dxd = DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData=sampleTable,
    #design= ~ sample + exon + condition:exon + adar:exon  + adar:condition:exon,
    design= ~ sample + exon + condition:exon + adar:exon,
    flattenedfile=flattenedFile )

  ## ----sizeFactors1----------------------------------------------------------
  print("get size factors")

  dxd = estimateSizeFactors( dxd )

  BPPARAM = MulticoreParam(workers=ncores)
  ## ----estDisp1--------------------------------------------------------------


   formulaFullModel    =  ~ sample + exon + adar:exon + condition:exon
   formulaReducedModel =  ~ sample + exon + adar:exon

  #formulaFullModel    =  ~ sample + exon + condition:exon + adar:exon + adar:condition:exon
  #formulaReducedModel =  ~ sample + exon + adar:exon

  print("estimate dispersion")
  dxd = estimateDispersions( dxd ,fitType="local",BPPARAM=BPPARAM,formula = formulaFullModel)

  ## ----testForDEU1,cache=TRUE------------------------------------------------
  print("test for DEU")


  #formulaFullModel    =  ~ sample + exon + adar + condition:exon
  #formulaReducedModel =  ~ sample + exon + adar

  dxd = testForDEU( dxd ,BPPARAM=BPPARAM,fullModel = formulaFullModel,reducedModel = formulaReducedModel)

  ## ----estFC,cache=TRUE------------------------------------------------------
  print("estimate fold changes")
  dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition",BPPARAM=BPPARAM)

  print("returning output")

  res <- DEXSeqResults(dxd)

  print("saving to out directory")
  save(res,file=paste0(out_dir,"/dxd_results.Rdat"))
  return(res)
}








###############################################################################
###############################################################################
####Helper function: find potential edits (DEXseq independent)#################
###############################################################################
###############################################################################

getEditsVsReference <- function(data_list,refGR=locsGR,stranded=F,my_edits=my_edits,mytype="treat",design_vector=design_vector,min_ref=3,min_targ=3,max_other=2)
{
  names(refGR) <- paste(as.vector(seqnames(refGR)),start(refGR),sep="_")
  dnames <- row.names(data_list[[1]])
  if(stranded)
  {
    the.chr <- gsub("(.+)_[0-9]+,[-|+]","\\1",dnames)
    the.pos <- as.numeric(gsub(".+_([0-9]+),[-|+]","\\1",dnames))
    the.strand <- gsub(".+_[0-9]+,([-|+])","\\1",dnames)
    refGR <- refGR[paste(the.chr,the.pos,sep="_")]
    pos.dat <- data.frame("chr"=the.chr,"pos"=the.pos,"strand"=the.strand)
    pos.dat$ref <- rep(NA,nrow(pos.dat))
    pos.dat[pos.dat$strand=="+",]$ref <- as.vector(refGR[paste(pos.dat[pos.dat$strand=="+",]$chr,pos.dat[pos.dat$strand=="+",]$pos,sep="_")]$ref)
    pos.dat[pos.dat$strand=="-",]$ref <- as.vector(refGR[paste(pos.dat[pos.dat$strand=="-",]$chr,pos.dat[pos.dat$strand=="-",]$pos,sep="_")]$ref)
    switch.1 <- c("A","T","G","C")
    names(switch.1) <- c("a","t","g","c")
    switch.2 <- c("A","T","G","C")
    names(switch.2) <- c("T","A","C","G")
    pos.dat$ref[pos.dat$ref %in% names(switch.1)] <- switch.1[pos.dat$ref[pos.dat$ref %in% names(switch.1)]]
    pos.dat$ref[pos.dat$strand=="-"] <- switch.2[pos.dat$ref[pos.dat$strand=="-"]]
    row.names(pos.dat) <- dnames
  }else{
    the.chr <- gsub("(.+)_[0-9]+","\\1",dnames)
    the.pos <- as.numeric(gsub(".+_([0-9]+)","\\1",dnames))
    refGR <- refGR[paste(the.chr,the.pos,sep="_")]
    pos.dat <- data.frame("chr"=the.chr,"pos"=the.pos)
    pos.dat$ref <- rep(NA,nrow(pos.dat))
    row.names(pos.dat) <- dnames
    pos.dat$ref <- as.vector(refGR$ref)
    switch.1 <- c("A","T","G","C")
    names(switch.1) <- c("a","t","g","c")
    pos.dat$ref[pos.dat$ref %in% names(switch.1)] <- switch.1[as.vector(pos.dat$ref[pos.dat$ref %in% names(switch.1)])]
  }

  poss_list <- list()
  type <- mytype


  for(i in 1:nrow(my_edits))
  {
    aref <- my_edits[i,1]
    atarg <- my_edits[i,2]
    baseother <- c("A","T","C","G")[!(c("A","T","C","G") %in% c(aref,atarg))]

    aref_count <- (do.call(cbind,lapply(data_list,function(x) x[,aref])))
    atarg_count <- (do.call(cbind,lapply(data_list,function(x) x[,atarg])))
    baseother1_count <- (do.call(cbind,lapply(data_list,function(x) x[,baseother[1]])))
    baseother2_count <- (do.call(cbind,lapply(data_list,function(x) x[,baseother[2]])))

    row.names(aref_count) <- row.names(data_list[[1]])
    row.names(atarg_count) <- row.names(data_list[[1]])
    row.names(baseother1_count) <- row.names(data_list[[1]])
    row.names(baseother2_count) <- row.names(data_list[[1]])

    aref_count_control <- aref_count[,names(design_vector[design_vector==type])]
    atarg_count_control <- atarg_count[,names(design_vector[design_vector==type])]
    baseother1_count_control <- baseother1_count[,names(design_vector[design_vector==type])]
    baseother2_count_control <- baseother2_count[,names(design_vector[design_vector==type])]

    x <- aref_count_control[which(pos.dat$ref==aref),]
    edit_ref <- apply(x,1,function(x) length(x[x>0]))
    edit_ref <- names(edit_ref[edit_ref>=min_ref])

    x <- atarg_count_control[which(pos.dat$ref==aref),]
    edit_targ <- apply(x,1,function(x) length(x[x>0]))
    edit_targ <- names(edit_targ[edit_targ>=min_targ])

    x <- baseother1_count_control[which(pos.dat$ref==aref),]
    edit_other1 <- apply(x,1,function(x) length(x[x>0]))
    edit_other1 <- names(edit_other1[edit_other1<=max_other])

    x <- baseother2_count_control[which(pos.dat$ref==aref),]
    edit_other2 <- apply(x,1,function(x) length(x[x>0]))
    edit_other2 <- names(edit_other2[edit_other2<=max_other])

    poss_sites <- Reduce(intersect,list(edit_ref,edit_targ,edit_other1,edit_other2))
    poss_list[[paste(aref,atarg,sep=":")]] <- (poss_sites)
  }
  return(poss_list)
}


###############################################################################
###############################################################################
####EDITING PROPORTIONS VS REFERENCE###########################################
###############################################################################
###############################################################################

editingProportionsVsReference <- function(data_list,refGR=locsGR,stranded=F)
{
  names(refGR) <- paste(as.vector(seqnames(refGR)),start(refGR),sep="_")
  dnames <- row.names(data_list[[1]])

  if(stranded)
  {
    the.chr <- gsub("(.+)_[0-9]+,[-|+]","\\1",dnames)
    the.pos <- as.numeric(gsub(".+_([0-9]+),[-|+]","\\1",dnames))
    the.strand <- gsub(".+_[0-9]+,([-|+])","\\1",dnames)
    refGR <- refGR[paste(the.chr,the.pos,sep="_")]
    pos.dat <- data.frame("chr"=the.chr,"pos"=the.pos,"strand"=the.strand)
    pos.dat$ref <- rep(NA,nrow(pos.dat))
    pos.dat[pos.dat$strand=="+",]$ref <- as.vector(refGR[paste(pos.dat[pos.dat$strand=="+",]$chr,pos.dat[pos.dat$strand=="+",]$pos,sep="_")]$ref)
    pos.dat[pos.dat$strand=="-",]$ref <- as.vector(refGR[paste(pos.dat[pos.dat$strand=="-",]$chr,pos.dat[pos.dat$strand=="-",]$pos,sep="_")]$ref)
    switch.1 <- c("A","T","G","C")
    names(switch.1) <- c("a","t","g","c")
    switch.2 <- c("A","T","G","C")
    names(switch.2) <- c("T","A","C","G")
    pos.dat$ref[pos.dat$ref %in% names(switch.1)] <- switch.1[pos.dat$ref[pos.dat$ref %in% names(switch.1)]]
    pos.dat$ref[pos.dat$strand=="-"] <- switch.2[pos.dat$ref[pos.dat$strand=="-"]]
    row.names(pos.dat) <- dnames
  }else{
    # the.chr <- gsub("(.+)_[0-9]+","\\1",dnames)
    # the.pos <- as.numeric(gsub(".+_([0-9]+)","\\1",dnames))
    # refGR <- refGR[paste(the.chr,the.pos,sep="_")]
    # pos.dat <- data.frame("chr"=the.chr,"pos"=the.pos)
    # pos.dat$ref <- rep(NA,nrow(pos.dat))
    # row.names(pos.dat) <- dnames
    # pos.dat$ref <- refGR$ref

    the.chr <- as.vector(seqnames(refGR))
    the.pos <- start(refGR)
    pos.dat <- data.frame("chr"=the.chr,"pos"=the.pos)
    pos.dat$ref <- rep(NA,nrow(pos.dat))
    row.names(pos.dat) <- dnames
    pos.dat$ref <- refGR$ref
  }

  A_count <- rowSums(do.call(cbind,lapply(data_list,function(x) x[,"A"])))
  T_count <- rowSums(do.call(cbind,lapply(data_list,function(x) x[,"T"])))
  G_count <- rowSums(do.call(cbind,lapply(data_list,function(x) x[,"G"])))
  C_count <- rowSums(do.call(cbind,lapply(data_list,function(x) x[,"C"])))

  my.df <- data.frame("A"=A_count,"T"=T_count,"G"=G_count,"C"=C_count)
  ref.expr <- rep(0,nrow(my.df))
  ref.expr[pos.dat$ref=="A"] <- my.df[pos.dat$ref=="A","A"]
  ref.expr[pos.dat$ref=="T"] <- my.df[pos.dat$ref=="T","T"]
  ref.expr[pos.dat$ref=="C"] <- my.df[pos.dat$ref=="C","C"]
  ref.expr[pos.dat$ref=="G"] <- my.df[pos.dat$ref=="G","G"]

  tot.expr <- rowSums(my.df)
  edit.A <- my.df[,"A"]/tot.expr
  edit.T <- my.df[,"T"]/tot.expr
  edit.C <- my.df[,"C"]/tot.expr
  edit.G <- my.df[,"G"]/tot.expr

  pos.dat$total <- tot.expr

  pos.dat$edit.A <- edit.A
  pos.dat$edit.T <- edit.T
  pos.dat$edit.C <- edit.C
  pos.dat$edit.G <- edit.G

  return(pos.dat)
}
