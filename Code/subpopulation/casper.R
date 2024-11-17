
#################
## 10x_3cl_310 ##
#################

library(Seurat)
library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Results/Tian/casper/10x_3cl_310/baf/'
output <- './Results/Tian/casper/10x_3cl_310/results/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
gene_exprs <- Read10X('./Data//sc_10x_310/outs/filtered_feature_bc_matrix',gene.column=1)
cell_id <- read.csv('./Results/Tian/cluster/sc_10x_310.csv')
gene_exprs <- gene_exprs[,as.character(cell_id[,1])]
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_exprs), ishg19=F, centromere)

gene_hcc1395 <- as.matrix(gene_exprs)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
gene_exprs <- t(cbind(gene_hcc1395[overlap_gene,],gene_hcc1395bl[overlap_gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
index <- match(annotation$Gene,rownames(gene_exprs))
gene_exprs <- gene_exprs[index[!is.na(index)],]
annotation <- annotation[!is.na(index),]

maf <- read.table(paste0(input,'casper.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))

all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))

l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))



#################
## 10x_3cl_220 ##
#################

library(Seurat)
library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Results/Tian/casper/10x_3cl_220/baf/'
output <- './Results/Tian/casper/10x_3cl_220/results/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
gene_exprs <- Read10X('./Data//sc_10x_220/outs/filtered_gene_bc_matrices/GRCh38',gene.column=1)
cell_id <- read.csv('./Results/Tian/cluster/sc_10x_220.csv')
gene_exprs <- gene_exprs[,as.character(cell_id[,1])]
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_exprs), ishg19=F, centromere)

gene_hcc1395 <- as.matrix(gene_exprs)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
gene_exprs <- t(cbind(gene_hcc1395[overlap_gene,],gene_hcc1395bl[overlap_gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
index <- match(annotation$Gene,rownames(gene_exprs))
gene_exprs <- gene_exprs[index[!is.na(index)],]
annotation <- annotation[!is.na(index),]

maf <- read.table(paste0(input,'casper.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))

all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))

l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))



#################
## 10x_5cl_310 ##
#################

library(Seurat)
library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Results/Tian/casper/10x_5cl_310/baf/'
output <- './Results/Tian/casper/10x_5cl_310/results/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
gene_exprs <- Read10X('./Data//sc_10x_5cl_310/outs/filtered_feature_bc_matrix',gene.column=1)
cell_id <- read.csv('./Results/Tian/cluster/sc_10x_5cl_310.csv')
gene_exprs <- gene_exprs[,as.character(cell_id[,1])]
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_exprs), ishg19=F, centromere)

gene_hcc1395 <- as.matrix(gene_exprs)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
gene_exprs <- t(cbind(gene_hcc1395[overlap_gene,],gene_hcc1395bl[overlap_gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
index <- match(annotation$Gene,rownames(gene_exprs))
gene_exprs <- gene_exprs[index[!is.na(index)],]
annotation <- annotation[!is.na(index),]

maf <- read.table(paste0(input,'casper.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))

all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))

l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))



#################
## 10x_5cl_201 ##
#################

library(Seurat)
library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Results/Tian/casper/10x_5cl_201/baf/'
output <- './Results/Tian/casper/10x_5cl_201/results/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
gene_exprs <- Read10X('./Data//sc_10x_5cl_201/outs/filtered_gene_bc_matrices/GRCh38',gene.column=1)
cell_id <- read.csv('./Results/Tian/cluster/sc_10x_5cl_201.csv')
gene_exprs <- gene_exprs[,as.character(cell_id[,1])]
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_exprs), ishg19=F, centromere)

gene_hcc1395 <- as.matrix(gene_exprs)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
gene_exprs <- t(cbind(gene_hcc1395[overlap_gene,],gene_hcc1395bl[overlap_gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
index <- match(annotation$Gene,rownames(gene_exprs))
gene_exprs <- gene_exprs[index[!is.na(index)],]
annotation <- annotation[!is.na(index),]

maf <- read.table(paste0(input,'casper.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))

all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))

l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


#############
## celseq2 ##
#############

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Results/Tian/casper/celseq2/baf/'
output <- './Results/Tian/casper/celseq2/results/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
gene_exprs <- read.table('/genomics/1_Projects/FDA_QC/Results/Tian/CelSeq2/counts.tsv.gz',header=T,row.names=1)
cell_id <- read.csv('./Results/Tian/cluster/celseq2.csv')
gene_exprs <- gene_exprs[,as.character(cell_id[,1])]
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_exprs), ishg19=F, centromere)

gene_hcc1395 <- gene_exprs
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
gene_exprs <- t(cbind(gene_hcc1395[overlap_gene,],gene_hcc1395bl[overlap_gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
index <- match(annotation$Gene,rownames(gene_exprs))
gene_exprs <- gene_exprs[index[!is.na(index)],]
annotation <- annotation[!is.na(index),]

maf <- read.table(paste0(input,'casper.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))

all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))

l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))



#############
## dropseq ##
#############

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Results/Tian/casper/dropseq/baf/'
output <- './Results/Tian/casper/dropseq/results/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
gene_exprs <- read.table('/genomics/1_Projects/FDA_QC/Results/Tian/Drop/counts.tsv.gz',header=T,row.names=1)
cell_id <- read.csv('./Results/Tian/cluster/dropseq.csv')
gene_exprs <- gene_exprs[,as.character(cell_id[,1])]
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_exprs), ishg19=F, centromere)

gene_hcc1395 <- gene_exprs
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
gene_exprs <- t(cbind(gene_hcc1395[overlap_gene,],gene_hcc1395bl[overlap_gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
index <- match(annotation$Gene,rownames(gene_exprs))
gene_exprs <- gene_exprs[index[!is.na(index)],]
annotation <- annotation[!is.na(index),]

maf <- read.table(paste0(input,'casper.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))

all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))

l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))



