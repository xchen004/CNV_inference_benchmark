## honeybager analysis ##
## expression analysis is based on breast tissue of GTEx reference data ##
## allele analysis is based on ExAC snp data ##



library(HoneyBADGER)

cell_cls <- read.csv('./Results/Tian/cluster/dropseq.csv',row.names=1)
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,1,0)))
load('./Results/Tian/hb/dropseq/hb.Rdata')

results_exprs <- t_hb$summarizeResults(geneBased=TRUE, alleleBased=FALSE)
tree_exprs <- t_hb$visualizeResults(geneBased=TRUE, alleleBased=FALSE, details=TRUE, margins=c(25,15))
cls_exprs <- cutree(tree_exprs$hc,k=3:12)
colnames(cls_exprs) <- paste0('pred_cls_',colnames(cls_exprs))
cell_cls <- cell_cls[rownames(cls_exprs),]
cls_exprs_all <- cbind(cell_cls,cls_exprs)

results_allele <- t_hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
tree_allele <- t_hb$visualizeResults(geneBased=FALSE, alleleBased=TRUE, details=TRUE, margins=c(25,15))
cls_allele <- cutree(tree_allele$hc,k=3:12)
colnames(cls_allele) <- paste0('pred_cls_',colnames(cls_allele))
cell_cls <- cell_cls[rownames(cls_allele),]
cls_allele_all <- cbind(cell_cls,cls_allele)

output <- './Results/Tian/hb/dropseq/summary/'
dir.create(output,recursive=T)
write.csv(cls_exprs_all,file=paste0(output,'exprs_cls.csv'),quote=F)
write.csv(cls_allele_all,file=paste0(output,'allele_cls.csv'),quote=F)


## cluster accuracy ##

# exprs
H1975:	78/84=92.86%
H2228:	54/56=96.43%
HCC827:	61/61=100%


# allele
H1975:	40/84=47.62%
H2228:	48/56=96.43%
HCC827:	59/61=96.72%


