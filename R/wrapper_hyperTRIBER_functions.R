###############################################################################
###############################################################################
####WRAPPER SCRIPT FOR WHOLE MODELLING PROCESS#################################
###############################################################################
###############################################################################

makeModel <- function(design_vector,my_data_list,locsGR,gtfGR,model_dir,save_dir,useADAR=F,adar=NULL,ncores=1,min_control=2,min_treat=2,quant.vec,FDR=0.01)
{
  data_list <- my_data_list[names(design_vector)]
  print(length(data_list))
  print(names(data_list))
  print(dim(data_list[[1]]))

  data_list <- restrict_data(data_list=data_list,design_vector=design_vector,min_samp_control=min_control,min_samp_treat=min_treat,min_count=1,edits_of_interest=my_edits)
  print(dim(data_list[[1]]))


  generateCountFiles(data_list = data_list,stranded=F,names_vec = row.names(data_list[[1]]),out_dir = model_dir,design_vector = design_vector)


  if(useADAR)
  {
    dxd.res <- make_test_ADAR(out_dir = model_dir,design_vector = design_vector,ncores=ncores,adar.vec=adar)

  }else{
    dxd.res <- make_test(out_dir = model_dir,design_vector = design_vector,ncores=ncores)

  }

  save(dxd.res,file=paste0(save_dir,"/dxd.res.Rdat"))


  posGR <- getHits(res = dxd.res,stranded=F,fdr = FDR,fold=-1000,addMeta = T,ncore=ncores,
                   include_ref = T,refGR = locsGR,edits_of_interest=my_edits,
                   design_vector=design_vector,data_list=data_list)
  posGR <- addStrandForHyperTRIBE(posGR)


  posGR <- addGenes(gtfGR = gtfGR,posGR=posGR,ncore=ncores,quant=quant.vec,assignStrand=T,geneids=ids)
  posGR$desc <- gene.desc[as.vector(posGR$gene),]$V3

  save(posGR,file=paste0(save_dir,"/posGR.Rdat"))
}
