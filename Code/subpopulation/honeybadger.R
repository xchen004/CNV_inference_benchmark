## honeybager analysis ##
## expression analysis is based on breast tissue of GTEx reference data ##
## allele analysis is based on ExAC snp data ##


#################
## 10x_3cl_310 ##
#################

library(Seurat)
library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

t_counts <- Read10X('./Data/sc_10x_310/outs/filtered_feature_bc_matrix',gene.column=1)
colnames(t_counts) <- paste0(colnames(t_counts),'-1')
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]


load('./Results/Tian/hb/10x_3cl_310/allele_base/snp_count.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="10x_3cl_310")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)

t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/Tian/hb/10x_3cl_310/hb.Rdata')




#################
## 10x_3cl_220 ##
#################

library(Seurat)
library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

t_counts <- Read10X('./Data/sc_10x_220/outs/filtered_gene_bc_matrices/GRCh38',gene.column=1)
colnames(t_counts) <- paste0(colnames(t_counts),'-1')
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]


load('./Results/Tian/hb/10x_3cl_220/allele_base/snp_count.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="10x_3cl_220")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)

t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/Tian/hb/10x_3cl_220/hb.Rdata')



#################
## 10x_5cl_310 ##
#################

library(Seurat)
library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

t_counts <- Read10X('./Data/sc_10x_5cl_310/outs/filtered_feature_bc_matrix',gene.column=1)
colnames(t_counts) <- paste0(colnames(t_counts),'-1')
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]


load('./Results/Tian/hb/10x_5cl_310/allele_base/snp_count.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="10x_5cl_310")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)

t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/Tian/hb/10x_5cl_310/hb.Rdata')




#################
## 10x_5cl_201 ##
#################

library(Seurat)
library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

t_counts <- Read10X('./Data/sc_10x_5cl_201/outs/filtered_gene_bc_matrices/GRCh38',gene.column=1)
colnames(t_counts) <- paste0(colnames(t_counts),'-1')
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]


load('./Results/Tian/hb/10x_5cl_201/allele_base/snp_count.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="10x_3cl_220")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)

t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/Tian/hb/10x_5cl_201/hb.Rdata')


#############
## celseq2 ##
#############


library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


t_counts <- read.table('./Data/Tian/CelSeq2/counts.tsv.gz',header=T,row.names=1)
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]


load('./Results/Tian/hb/celseq2/allele_base/snp_count.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="celseq2")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)

t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/Tian/hb/celseq2/hb.Rdata')




#############
## dropseq ##
#############


library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

t_counts <- read.table('./Data/Tian/Drop/counts.tsv.gz',header=T,row.names=1)
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]

load('./Results/Tian/hb/dropseq/allele_base/snp_count.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="dropseq")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)

t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/Tian/hb/dropseq/hb.Rdata')


















