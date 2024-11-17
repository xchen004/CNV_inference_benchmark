## honeybager analysis ##
## expression analysis is based on lung tissue of GTEx reference data ##
## allele analysis is based on ExAC snp data ##


library(HoneyBADGER)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

load('./Results/SCLC/Star_featureCounts/Ensembl_GRCh38/Gtf/gene_counts.rdata')
t_counts <- gene_counts
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


load('./Results/SCLC/hb/allele_base/snp_count.Rdata')
index <- match(paste0(colnames(results$cov),'-1'),colnames(t_cpm))
results$cov <- results$cov[,!is.na(index)]
results$ref <- results$ref[,!is.na(index)]
results$alt <- results$alt[,!is.na(index)]
t_cpm <- t_cpm[,index[!is.na(index)]]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
mart.obj <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host = "apr2019.archive.ensembl.org")

t_hb <- new('HoneyBADGER',name="SCLC")
t_hb$setGexpMats(t_cpm,t_ref_cpm,mart.obj,minMeanBoth=1,minMeanTest=3.5,minMeanRef=3.5)
t_hb$setMvFit(verbose=TRUE)
t_hb$setGexpDev(verbose=TRUE)
t_hb$calcGexpCnvBoundaries(init=TRUE, verbose=TRUE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=FALSE, verbose=FALSE)


t_hb$setAlleleMats(r.init=results$altCount, n.sc.init=results$cov, het.deviance.threshold=0.05)
t_hb$setGeneFactors(txdb)
t_hb$calcAlleleCnvBoundaries(init=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=FALSE, retestBoundSnps=TRUE, verbose=FALSE)
t_hb$retestIdentifiedCnvs(retestBoundGenes=TRUE, retestBoundSnps=TRUE, verbose=FALSE)

save(t_hb,file='./Results/SCLC/hb/hb.Rdata')



load('./Results/SCLC/hb/hb.Rdata')
results_allele <- t_hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
tree_allele <- t_hb$visualizeResults(geneBased=FALSE, alleleBased=TRUE, details=TRUE, margins=c(25,15))
cls_allele <- cutree(tree_allele$hc,3:12)
colnames(cls_allele) <- paste0('prer_cls_',colnames(cls_allele))
cell_cls <- matrix(c(rep('primary',93),rep('relapase',39)),ncol=1)
cls_allele_all <- cbind(cell_cls,cls_allele)

output <- './Results/SCLC/hb/summary/'
dir.create(output,recursive=T)
write.csv(cls_allele_all,file=paste0(output,'allele_cls.csv'),quote=F)









