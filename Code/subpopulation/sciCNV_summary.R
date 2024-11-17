
#################
## 10x_3cl_310 ##
#################

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/scicnv/10x_3cl_310/results/'
input2 <- './Results/Tian/infercnv/10x_3cl_310/'
output <- './Results/Tian/scicnv/10x_3cl_310/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_310.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,0,1)))

## cluster ##

s_score <- read.csv(paste0(input,'cnv_genes.csv'),row.names=1)
sample_anno <- read.table(paste0(input2,'/sample_anno.txt.gz'))
test_num <- table(sample_anno[,2])[2]
s_score2 <- round(s_score[,4:(test_num+3)],3)
s_score2 <- s_score2[apply(s_score2,1,sum)!=0,]
s_score_d <- dist(t(s_score2),method='euclidean')
s_score_fit <- hclust(s_score_d,method='ward.D')
save(s_score_fit,file=paste0(output,'s_score_fit.rdata'))
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
cls_all <- cbind(cell_cls,s_score_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)




#################
## 10x_3cl_220 ##
#################

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/scicnv/10x_3cl_220/results/'
input2 <- './Results/Tian/infercnv/10x_3cl_220/'
output <- './Results/Tian/scicnv/10x_3cl_220/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_220.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))

## cluster ##

s_score <- read.csv(paste0(input,'cnv_genes.csv'),row.names=1)
sample_anno <- read.table(paste0(input2,'/sample_anno.txt.gz'))
test_num <- table(sample_anno[,2])[2]
s_score2 <- round(s_score[,4:(test_num+3)],3)
s_score2 <- s_score2[apply(s_score2,1,sum)!=0,]
s_score_d <- dist(t(s_score2),method='euclidean')
s_score_fit <- hclust(s_score_d,method='ward.D')
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
cls_all <- cbind(cell_cls,s_score_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)





#################
## 10x_5cl_310 ##
#################

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/scicnv/10x_5cl_310/results/'
input2 <- './Results/Tian/infercnv/10x_5cl_310/'
output <- './Results/Tian/scicnv/10x_5cl_310/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_5cl_310.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,3,4,0,1)))

## cluster ##

s_score <- read.csv(paste0(input,'cnv_genes.csv'),row.names=1)
sample_anno <- read.table(paste0(input2,'/sample_anno.txt.gz'))
test_num <- table(sample_anno[,2])[2]
s_score2 <- round(s_score[,4:(test_num+3)],3)
s_score2 <- s_score2[apply(s_score2,1,sum)!=0,]
s_score_d <- dist(t(s_score2),method='euclidean')
s_score_fit <- hclust(s_score_d,method='ward.D')
s_score_cls <- cutree(s_score_fit,5:14,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
cls_all <- cbind(cell_cls,s_score_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)



#################
## 10x_5cl_201 ##
#################

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/scicnv/10x_5cl_201/results/'
input2 <- './Results/Tian/infercnv/10x_5cl_201/'
output <- './Results/Tian/scicnv/10x_5cl_201/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_5cl_201.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,3,4,0,1)))

## cluster ##

s_score <- read.csv(paste0(input,'cnv_genes.csv'),row.names=1)
sample_anno <- read.table(paste0(input2,'/sample_anno.txt.gz'))
test_num <- table(sample_anno[,2])[2]
s_score2 <- round(s_score[,4:(test_num+3)],3)
s_score2 <- s_score2[apply(s_score2,1,sum)!=0,]
s_score_d <- dist(t(s_score2),method='euclidean')
s_score_fit <- hclust(s_score_d,method='ward.D')
s_score_cls <- cutree(s_score_fit,5:14,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
cls_all <- cbind(cell_cls,s_score_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)



#################
## celseq2 ##
#################

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/scicnv/celseq2/results/'
input2 <- './Results/Tian/infercnv/celseq2/'
output <- './Results/Tian/scicnv/celseq2/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/celseq2.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))

## cluster ##

s_score <- read.csv(paste0(input,'cnv_genes.csv'),row.names=1)
sample_anno <- read.table(paste0(input2,'/sample_anno.txt.gz'))
test_num <- table(sample_anno[,2])[2]
s_score2 <- round(s_score[,4:(test_num+3)],3)
s_score2 <- s_score2[apply(s_score2,1,sum)!=0,]
s_score_d <- dist(t(s_score2),method='euclidean')
s_score_fit <- hclust(s_score_d,method='ward.D')
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
cls_all <- cbind(cell_cls,s_score_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)


#################
## dropseq ##
#################

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/scicnv/dropseq/results/'
input2 <- './Results/Tian/infercnv/dropseq/'
output <- './Results/Tian/scicnv/dropseq/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/dropseq.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,1,0)))

## cluster ##

s_score <- read.csv(paste0(input,'cnv_genes.csv'),row.names=1)
sample_anno <- read.table(paste0(input2,'/sample_anno.txt.gz'))
test_num <- table(sample_anno[,2])[2]
s_score2 <- round(s_score[,4:(test_num+3)],3)
s_score2 <- s_score2[apply(s_score2,1,sum)!=0,]
s_score_d <- dist(t(s_score2),method='euclidean')
s_score_fit <- hclust(s_score_d,method='ward.D')
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
cls_all <- cbind(cell_cls,s_score_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)



