# hyperTRIBER
R package for differential RNA editing analysis

## Summary
hyperTRIBER is an R package used for detecting sites which significant differental editing between conditions using transcriptomics based data. The pipeline was originally developed for hyperTRIBE (targets of RNA binding proteins identified by editing), which allows for the detection of transcripts bound by a given RNA binding protein (RBP) through the detection of edited sites between negative control samples, and cases where the RBP is fused to the hypercatyletic domain of ADAR. However, the approach is highly applicable for general differential editing set-ups.
 
## Installation

The package is installed in R with the following command:

```
library(devtools)
install_github("sarah-ku/hyperTRIBER")
```   

Prior to running the package in R, one needs to run the following command including the custom perl script `hyperTRIBE_mpileup2bases.pl` on their mapped .BAM files for their samples. For example, for a set up with 6 samples (Samp1 to Samp6), would the command look like this:

```wrap
samtools mpileup --max-depth 50000 -Q 30 --skip-indels -f ./reference/reference_genome.fa Samp1.sort.bam Samp2.sort.bam Samp3.sort.bam Samp4.sort.bam Samp5.sort.bam Samp6.sort.bam | perl hyperTRIBE_mpileup2bases.pl> baseCounts_from_mpileup.txt &
```
Note that the reference genome is necessary in order to run samtools mpileup.

## Running the pipeline

Here we prepare the experimental setup for our data
```
#load in our text file that we got as an output from 
data <- read.table("./files/baseCounts_from_mpileup.txt",header=F)
dim(data)

#create a genomics range of the reference genome base for each of the considered positions.
locsGR <- GRanges(Rle(data_drp$V1),IRanges(data_drp$V2,width=1),ref=data_drp$V3,names=paste(data_drp$V1,data_drp$V2,sep="_"))


samp.names <- c("Samp1","Samp2","Samp3","Samp4","Samp5","Samp6")


data_list <- extractCountData(data,samp.names,strand=F)


#now produce one design vector per experiment
design_vector <- c(Samp1 = "control", Samp2 = "control", Samp3 = "control", 
                                Samp4 = "treat", Samp5 = "treat", Samp6 = "treat")

table(design_vector)
```

## R code to Find Edits
This is the R code that filters our data further, reframes the data so it can be run through DEXseq package, and uses the statistical infromation acquired from DEXseq to find significant differential editing sites. Since we are using two different neural populations for our experiment, we expect to find a lot of significant sites.


### Specify our edits of interests.
```
my_edits <- rbind(c("A","G"),
                  c("G","A"),
                  c("T","C"),
                  c("C","T"),
                  c("A","T"),
                  c("T","A"),
                  c("G","T"),
                  c("T","G"),
                  c("G","C"),
                  c("C","G"),
                  c("A","C"),
                  c("C","A"))
```


### Preparations
```
project_id <- "drosophila"
data_list <- data_list[names(design_vector)]

#remember to create folders if they don't exist
model_dir <- "./results/model/"
save_dir <-  "./results/model/saved_output/"
```

### filter according to replication and minimum counts
This function filters out sites based on replication and minimum counts. This is very important as it reduces the run time of functions that we will be running later.
The argument <b>edits_of_interest</b> is used to specify the edits that will be looked at.
The argument <b>min_count</b> only keeps sites that have the minimum amount of the nucleotide of interest. Our minimum amount is set to 2.
The argument <b>min_samp_control</b> looks at the 1st base of each vector contained in <b>my_edits</b> and keeps only sites that contain the set minimum amount of counts in a set number of control samples (2 in our case) while <b>min_samp_treat</b> looks at the 2nd base of each vector and does the same for treatment samples.

For other experimental designs, these can of course be changed but for the sake of the drosophila experiment, we kept the values at 2 for each of the arguments.

```
data_list <- restrict_data(data_list=data_list,design_vector=design_vector,min_samp_control=2,min_samp_treat=2,min_count=2,edits_of_interest=my_edits)
```

 
### Primary comonent analysis
```
my_pca <- makeEditPCA(data_list_restricted,editTypes = editTypes,refGR = locsGR,design_vector = design_vector,my_title = "",refBased = T,stranded=T)
```

### Make model with DEXseq
This creates countfiles that can be used by DEXseq and a fake annotation files used in the modelling process to generate the statistical framework to be used to find significant sites.
```
generateCountFiles(data_list = data_list,stranded=F,names_vec = row.names(data_list[[1]]),out_dir = model_dir,design_vector = design_vector)
dxd.res <- make_test(out_dir = model_dir,design_vector = design_vector,ncores=10)
save(dxd.res,file=paste0(save_dir,"/dxd.res.Rdat"))
```
### Compile hits table from model
Here we use FDR to find significant "hits" from our DEXseq data and compile it into a table.
```
posGR <- getHits(res = dxd.res,stranded=F,fdr = 0.1,fold=-1000,addMeta = T,ncore=10,include_ref = T,refGR = locsGR,edits_of_interest=my_edits,design_vector=design_vector,data_list=data_list)
```

## Annotation

Now that we had our significant hits, we wanted to annotate them to the genome to identify locations and structure of our hits. The GTF file for the annotations was fetched from https://www.ensembl.org/Drosophila_melanogaster/Info/Index
```
#open annotations
gtf <- rtracklayer::import("./reference/Drosophila_melanogaster.BDGP6.86.gtf")
gtfGR <- gtf

ids <- tapply(gtf$gene_id,gtf$gene_name,function(x) x[1])

quant.mat <- tpm.mat[,names(design_vector[design_vector=="control"])]
quant.vec <- rowMeans(quant.mat)

HACK:
at the moment, we are working on fixing this, but the annotation will not work unless these commands are run:
tmp <- as.vector(gtf$type)
tmp[tmp=="three_prime_utr"] <- "3UTR"
tmp[tmp=="five_prime_utr"] <- "5UTR"
tmp[tmp=="transcript"] <- "CDS"
gtf$type <- as.factor(tmp)

#Annotate
posGR <- addGenes(gtfGR = gtf,posGR = posGR,ncore = 10,quant = quant.vec,assignStrand = F,geneids = ids)
save(posGR,file=paste0(save_dir,"/posGR_lenient.Rdat"))
```

## Works cited:
1. Illuminating spatial A-to-I RNA editing signatures within the Drosophila brain
Anne L. Sapiro, Anat Shmueli, Gilbert Lee Henry, Qin Li, Tali Shalit, Orly Yaron, Yoav Paas, Jin Billy Li, Galit Shohat-Ophir
Proceedings of the National Academy of Sciences Feb 2019, 116 (6) 2318-2327; DOI: 10.1073/pnas.1811768116
