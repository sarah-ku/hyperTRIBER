plotGeneBody <- function(myGR,myGR_list=NULL,list_style=F,my.alpha=.25,my.adjust=1)
{

  if(list_style)
  {
    tp.list <- list()
    for(k in 1:length(myGR_list))
    {

      myGR <- myGR_list[[k]]
      #myGR <- myGR[!(myGR$out_of_range)]
      tt <- unlist(lapply(strsplit(myGR$transcript.types,","),function(x) x[1]))
      fp <- unlist(lapply(strsplit(myGR$feature_prop,","),function(x) x[1]))
      length(tt)
      length(fp)
      fdf <- data.frame("tt"=tt,"fp"=as.numeric(fp))
      fdf <- fdf[!is.na(fdf$fp),]
      fdf <- fdf[fdf$fp <= 1 & fdf$fp>= 0,]


      vals_3UTR <- fdf[fdf$tt=="3UTR",2]
      vals_5UTR <- fdf[fdf$tt=="5UTR",2]
      vals_CDS <- fdf[fdf$tt=="CDS",2]

      #to.add <- c( 0.5*(vals_5UTR),0.5*(vals_CDS)+0.5,0.5*(vals_3UTR)+1.5)
      to.add <- c((vals_5UTR),(vals_CDS)+1,(vals_3UTR)+2)

      tp.list[[names(myGR_list)[k]]] <- data.frame("variable"=rep(names(myGR_list)[k],length(to.add)),"value"=to.add)
    }

    df <- do.call(rbind,tp.list)
    colnames(df) <- c("variable","value")

    PC <- ggplot(df, aes(value, fill=variable)) +
      geom_density(alpha=my.alpha,adjust=my.adjust) +
      geom_vline(xintercept=1,lty=2) +
      geom_vline(xintercept=2,lty=2) +
      xlim(0,3) +
      theme_minimal()

    PC
  }else{
    pos_3UTR <- which(!is.na(myGR$feature_locs[,"3UTR"]))
    pos_5UTR <-which(!is.na(myGR$feature_locs[,"5UTR"]))
    pos_CDS <-which(!is.na(myGR$feature_locs[,"CDS"]))

    #remove
    #pos_3UTR <- pos_3UTR[!(pos_3UTR %in% c(pos_CDS,pos_5UTR))]


    vals_3UTR <- myGR$feature_locs[pos_3UTR,"3UTR"]
    vals_5UTR <- myGR$feature_locs[pos_5UTR,"5UTR"]
    vals_CDS <- myGR$feature_locs[pos_CDS,"CDS"]

    tp.dat <- c(vals_5UTR,vals_CDS+1,vals_3UTR+2)

    df <- melt(data.frame(tp.dat))

    PC <- ggplot(df, aes(value, fill=variable)) +
      geom_density(alpha=my.alpha,adjust=my.adjust) +
      geom_vline(xintercept=1,lty=2) +
      geom_vline(xintercept=2,lty=2) +
      xlim(0,3) +
      theme_minimal()

    PC
  }
}

plotEdits <- function(posGR)
{
  my.df <- table(paste(posGR$meta$ref,posGR$meta$targ,sep=":"))
  barplot(my.df)
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


makeEditPCA <- function(data_list,editTypes=my_edits,refGR,design_vector,my_title="",refBased=T,posGR=NULL,stranded=F,give.mat=F)
{
  s.names <- names(data_list)


  if(refBased)
  {

    refGR$ref <- toupper(refGR$ref)

    if(stranded)
    {
      conv_vec <- c("A","T","C","G")
      names(conv_vec) <- c("T","A","G","C")

      names(refGR) <- paste0(refGR$names,",+")
      refGR2 <- refGR
      names(refGR2) <- paste0(refGR$names,",-")
      refGR2$ref <- conv_vec[as.vector(refGR2$ref)]
      refGR <- unlist(GRangesList(refGR,refGR2))
      refGR <- refGR[row.names(data_list[[1]])]

    }else{
      names(refGR) <- refGR$names
      refGR <- refGR[row.names(data_list[[1]])]
    }


    tmp_list <- list()
    for(i in 1:nrow(editTypes))
    {
      editType <- editTypes[i,]
      data_red <- lapply(data_list,function(x) x[names(refGR[refGR$ref==editType[1]]),editType])
      data_red_matrix <- do.call(cbind,lapply(data_red,function(x) x[,2]/(x[,1]+x[,2])))
      data_red_matrix[is.na(data_red_matrix)] <- 0
      tmp_list[[i]] <- data_red_matrix
    }

    data_red_matrix <- do.call(rbind,tmp_list)
    data_red_matrix <- data_red_matrix[apply(data_red_matrix,1,max)>0,]
    print(dim(data_red_matrix))

    pca<-prcomp(data_red_matrix,center=T,scale.=T)
    pca_comps <- as.data.frame(pca$rotation)

    #type_colors <- c("red","blue")
  }else{
    tmp_list <- list()
    for(i in 1:nrow(editTypes))
    {
      editType <- editTypes[i,]
      data_red <- lapply(data_list,function(x) x[names(posGR[posGR$ref==editType[1]]),editType])
      data_red_matrix <- do.call(cbind,lapply(data_red,function(x) x[,2]/(x[,1]+x[,2])))
      data_red_matrix[is.na(data_red_matrix)] <- 0
      data_red_matrix <- data_red_matrix[apply(data_red_matrix,1,max)>0,]

      tmp_list[[i]] <- data_red_matrix
    }

    data_red_matrix <- do.call(rbind,tmp_list)

    pca<-prcomp(data_red_matrix,center=T,scale.=T)
    pca_comps <- as.data.frame(pca$rotation)

  }

  type <- design_vector[names(data_list)]
  pca_plot_quantification <- ggplot(data=pca_comps, aes(x=PC1, y=PC2,fill=factor(type), label=s.names)) +
    geom_point(size=6,pch=21)+
    labs(color="Type")+
    geom_text(aes(label=s.names), size=3, hjust=0.9, vjust=1.5,col="black")+
    xlab(paste("PC1",round(summary(pca)$importance[2, 1], 4) * 100, "% of variance")) +
    ylab(paste("PC2", round(summary(pca)$importance[2, 2], 4) * 100, "% of variance")) +
    ggtitle(my_title) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20), axis.text=element_text(size=20)) + theme_classic()

  if(give.mat)
  {
    return(list(pca_plot_quantification,data_red_matrix))
  }else{
    return(pca_plot_quantification)
  }
}



editTypesPlot <- function(posGR,edit_types=my_edits)
{
  evec <- rep(0,nrow(edit_types))
  names(evec) <- paste(edit_types[,1],edit_types[,2],sep=":")
  tp <- table(paste(posGR$meta$ref,posGR$meta$targ,sep=":"))
  evec[names(tp)] <- tp
  evec <- rev(sort(evec))
  my.df <- data.frame("type"=names(evec),"count"=as.vector(evec))
  p1 <- ggplot(my.df,aes(x=factor(type,levels = my.df$type),y=count)) + geom_bar(stat="identity") + theme_bw() +  xlab("") +theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("Edit type counts")
  return(p1)
}
