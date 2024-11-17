## version 1 ##
## expression analysis is based on scRNAseq reference data ##
## allele analysis is based on ExAC snp data ##


#########
## 10x ##
#########

library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

load('./Data/10x/gene_counts/HCC1395/gene_counts_featureCounts.rdata')
t_counts <- gene_counts
colnames(t_counts) <- paste0(colnames(t_counts),'-1')
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Data/10x/gene_counts/HCC1395BL/gene_counts_featureCounts.rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]


load('./Results/10x/allele_base/snp_count_hcc1395bl_2.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
results$altCount <- results$altCount[,!is.na(index)]
results$refCount <- results$refCount[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="10x")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)

t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
#t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/10x/all_cell/hb/v1/t_hb.Rdata')




#########
## c1-HT ##
#########

library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

load('./Data/C1_ceber/gene_counts/HCC1395/gene_counts_featureCounts.rdata')
t_counts <- gene_counts
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Data/C1_ceber/gene_counts/HCC1395BL/gene_counts_featureCounts.rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]


load('./Results/c1_fda/allele_base/snp_count.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
results$altCount <- results$altCount[,!is.na(index)]
results$refCount <- results$refCount[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="C1-HT")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)

t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/c1_fda/all_cell/hb/v1/t_hb.Rdata')





#########
## icell8 ##
#########


library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

load('./Data/WaferGen/gene_counts/read_len_SE/150bp/HCC1395/gene_counts_featureCounts.rdata')
t_counts <- gene_counts
t_counts <- t(t_counts)
t_cpm <- log(t_counts/apply(t_counts,1,sum)*1e6+1)
t_cpm <- as.matrix(t(t_cpm))

load('./Data/WaferGen/gene_counts/read_len_SE/150bp/HCC1395BL/gene_counts_featureCounts.rdata')
gene_counts <- t(gene_counts)
ref_cpm <- t(log(gene_counts/apply(gene_counts,1,sum)*1e6+1))
index <- match(rownames(t_cpm),rownames(ref_cpm))
t_ref_cpm <- ref_cpm[index[!is.na(index)],]
t_cpm <- t_cpm[!is.na(index),]

geneid_all <- read.csv('./Anno/10x_gene_anno.csv')
rownames(geneid_all) <- geneid_all[,1]
t_geneid_all <- geneid_all[rownames(t_cpm),]
rownames(t_cpm) <- rownames(t_ref_cpm) <- t_geneid_all[,2]


load('./Results/WaferGen/read_len_SE/150bp/allele_base/snp_count.Rdata')
index <- match(colnames(results$cov),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
results$altCount <- results$altCount[,!is.na(index)]
results$refCount <- results$refCount[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="C1-HT")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=2.6)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)

t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/WaferGen/all_cell/hb/v1/t_hb.Rdata')









