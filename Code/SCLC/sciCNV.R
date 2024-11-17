
library(base)
library(robustbase)
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(qlcMatrix)
library(svd)
library(ggplot2)
library(ggridges)
library(viridis)
require(scales)
library(RColorBrewer)
library(beanplot)
 library(gridGraphics)
library(Rtsne)
library(reticulate)
library(umap)

path.code <- "./Code/sciCNV/sciCNV-Analysis/"
source(file.path(path.code, "Mito_umi_gn.R"))
source(file.path(path.code, "RTAM_normalization.R"))
source(file.path(path.code, "sciCNV.R"))
source(file.path(path.code, "Scaling_CNV.R"))
source(file.path(path.code, "CNV_score.R"))
source(file.path(path.code, "sciCNV.R"))
source(file.path(path.code, "Sketch_AveCNV.R"))
source(file.path(path.code, "CNV_htmp_glist.R"))
source(file.path(path.code, "CNV_htmp_gloc.R"))
source(file.path(path.code, "Opt_MeanSD_RTAM1.R"))
source(file.path(path.code, "Opt_MeanSD_RTAM2.R"))
source(file.path(path.code, "heatmap_break_glist.R"))
source(file.path(path.code, "heatmap_break_gloc.R"))



get_cnv_genes <- function(gene_exprs, test_num)
{
	gene_exprs <- gene_exprs
	gene_norm <- RTAM_normalization(mat = gene_exprs, method = "RTAM2", Min_nGn  = 250, Optimizing = FALSE)
	No.test <- test_num
	tst.index  <- seq(1, No.test , 1)                      # No of test cells
	ctrl.index <- seq(No.test+1, ncol(gene_norm), 1)      # No of controcl cells

	## generating infered-CNV data for (test and/or control) cells 
	CNV_data <- sciCNV(norm.mat = gene_norm, No.test = No.test, sharpness  = 1, baseline_adj  = FALSE,  baseline = 0)

	#####################################################################
	##  Scaling the sciCNV curves and setting a noise filter threshold ##
	#####################################################################

	## Scaling CNV-curves to adjust one copy number gain/losses to height +1/-1 if applicable
	CNV.data.scaled <- Scaling_CNV(V7Alt = CNV_data, n.TestCells = No.test, scaling.factor = 1)

	## Defining M_NF as Noise-Free Matrix of test cells (and control cells)
	M_NF <- CNV.data.scaled

	#######  Noise Filteration after scaling
	noise.thr = 0.4   # Noise threshold
	for(w in 1:ncol(M_NF) ){
  	for(j in 1:nrow(M_NF)){
    	if( (M_NF[j,w] > -noise.thr) && (M_NF[j,w] < noise.thr)  ){
      	M_NF[j,w] <- 0
    	}}}
	M_NF <- as.matrix(M_NF)

	#### Taking Square Rate
	for(w in 1:ncol(M_NF)){
	  for(j in 1:nrow(M_NF)){
	    if (M_NF[j,w] > 0){
	      M_NF[j,w]<- sqrt( as.numeric(M_NF[j,w]))
	    } else if (M_NF[j,w] < 0){
	      M_NF[j,w]<- -sqrt( -as.numeric(M_NF[j,w]))
	    }
	  }
	}

	rownames(M_NF) <- rownames(CNV.data.scaled)
	colnames(M_NF) <- c(colnames(CNV.data.scaled)[-length(colnames(CNV.data.scaled))], "AveTest")

	## Assigning chromosome number to each gene sorted based on chromosome number, 
	## starts and ends to sketch the average sciCNV curve of test cells

	Gen.loc <- gene_anno2
	Specific_genes <- which( Gen.loc[, 1] %in% rownames(CNV.data.scaled))
	Assoc.Chr <-  Gen.loc[Specific_genes, 2]

	#### Finalizing the sciCNV-matrix by attaching gene-name and chromosome number lists
	M_NF1 <- cbind.data.frame(as.matrix(Gen.loc[Specific_genes,]), M_NF)
	colnames(M_NF1) <- c(colnames(Gen.loc), colnames(gene_norm), "Ave test")
	return(M_NF1)
}


############
## scicnv ##
############

input <- './Results/SCLC/infercnv/'
input2 <- './Anno/'
output <- './Results/SCLC/scicnv/results/'
dir.create(output,recursive=T)

sample_anno <- read.table(paste0(input,'/sample_anno.txt.gz'))
gene_anno <- read.csv('./Anno/10x_gene_anno.csv')
gene_anno2 <- read.csv('./Code/sciCNV/Dataset/10XGenomics_gen_pos_GRCh38-1.2.0.txt',header=T,sep='\t')
test_num <- as.numeric(table(sample_anno[,2])[2])

gene_exprs <- read.table(paste0(input,'/raw_counts.txt.gz'),header=T,sep='\t')
index <- match(rownames(gene_exprs),gene_anno[,1])
gene_exprs <- gene_exprs[!is.na(index),]
rownames(gene_exprs) <- gene_anno[index[!is.na(index)],2]
index <- match(rownames(gene_exprs),gene_anno2[,1])
gene_exprs <- gene_exprs[!is.na(index),]
cnv_genes <- get_cnv_genes(gene_exprs,test_num)
write.csv(cnv_genes,file=paste0(output,'cnv_genes.csv'),quote=F,row.names=F)

input <- './Results/SCLC/scicnv/results/'
cnv_genes <- read.csv(paste0(input,'cnv_genes.csv'))

pri_cnv <- cnv_genes[,grep('PriT',colnames(cnv_genes))]
re_cnv <- cnv_genes[,grep('Re',colnames(cnv_genes))]
rownames(pri_cnv) <- rownames(re_cnv) <- cnv_genes[,1]
output <- './Results/SCLC/scicnv/Primary/'
output2 <- './Results/SCLC/scicnv/Relapse/'
dir.create(output,recursive=T)
dir.create(output2,recursive=T)
write.csv(pri_cnv,file=paste0(output,'cnv_genes.csv'),quote=F)
write.csv(re_cnv,file=paste0(output2,'cnv_genes.csv'),quote=F)










