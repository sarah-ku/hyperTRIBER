









genelevelStatistics <- function(posGR)
{
  #hits per gene
  
  
  #hits per expression level
  
  #editing proportions per gene
}

peakRelativeDistance <- function(ref=posGR,targ=iCLIP[[3]])
{
  
}




FCvsEP <- function(posGR)
{
  plot(posGR$meta$fold_change,posGR$meta$prop,pch=16,cex=0.5,xlab="Fold change",ylab="Editing proportion")
  points(posGR$meta$fold_change[posGR$meta$padj<0.1],posGR$meta$prop[posGR$meta$padj<0.1],col="red",pch=16)
  abline(v=0,lty=2)
}





parControlVSTreat <- function(posGR)
{
  plot(posGR$meta$control_par,posGR$meta$treat_par,cex=0.5)
  abline(0,1,lty=2)
}




