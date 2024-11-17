
#################
## 10x_3cl_310 ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/casper/10x_3cl_310/results/'
output <- './Results/Tian/casper/10x_3cl_310/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_310.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,0,1)))

## cluster by large scale cnv ##

load(paste0(input,'/large_cnv.rdata'))
l_cnv <- l_cnv[grep('_1',rownames(l_cnv)),]
l_cnv_d <- dist(l_cnv,method='euclidean')
l_cnv_fit <- hclust(l_cnv_d,method='ward.D')
l_cnv_cls <- cutree(l_cnv_fit,3:12,order_clusters_as_data=F)
colnames(l_cnv_cls) <- paste0('pred_cls_',colnames(l_cnv_cls))
cell_cls <- cell_cls[rownames(l_cnv_cls),]
l_cnv_cls_all <- cbind(cell_cls,l_cnv_cls)
write.csv(l_cnv_cls_all,file=paste0(output,'large_cnv_cls.csv'),quote=F)



## cluster by gene based cnv ##

load(paste0(input,'/casper_segment.rdata'))
load(paste0(input,'/casper_final.rdata'))
gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
rna.matrix <- t(rna.matrix)
gene_d <- dist(rna.matrix,method='euclidean')
gene_fit <- hclust(gene_d,method='ward.D')
gene_cls <- cutree(gene_fit,3:12,order_clusters_as_data=F)
colnames(gene_cls) <- paste0('pred_cls_',colnames(gene_cls))
cell_cls <- cell_cls[rownames(gene_cls),]
gene_cls_all <- cbind(cell_cls,gene_cls)
write.csv(gene_cls_all,file=paste0(output,'gene_cls.csv'),quote=F)


## cluster smoothed gene expression ##

obj <- final.objects[[9]]
s_exprs <- t(obj@control.normalized[[3]])
s_exprs <- s_exprs[grep('_1',rownames(s_exprs)),]
s_exprs_d <- dist(s_exprs,method='euclidean')
s_exprs_fit <- hclust(s_exprs_d,method='ward.D')
s_exprs_cls <- cutree(s_exprs_fit,3:12,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]
cls_all <- cbind(cell_cls,s_exprs_cls)
s_exprs_cls3 <- table(cell_cls[,2],s_exprs_cls[,1])
s_exprs_cls4 <- table(cell_cls[,2],s_exprs_cls[,2])
s_exprs_cls5 <- table(cell_cls[,2],s_exprs_cls[,3])
s_exprs_cls6 <- table(cell_cls[,2],s_exprs_cls[,4])
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)
write.csv(s_exprs_cls3,file=paste0(output,'cls3_table.csv'),quote=F)
write.csv(s_exprs_cls4,file=paste0(output,'cls4_table.csv'),quote=F)
write.csv(s_exprs_cls5,file=paste0(output,'cls5_table.csv'),quote=F)
write.csv(s_exprs_cls6,file=paste0(output,'cls6_table.csv'),quote=F)

col <- brewer.pal(8,"Set1")
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_exprs_fit)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=6, value=col[c(2,2,2,3,1,4)]) %>% plot
colored_bars(colors=bar_col,y_shift=-200,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975'),col=col[1:3],pch=15,bty='n')
dev.off()


## cls3 accuracy ##

H1975:	327/346=94.51
H2228:	325/327=99.39%
HCC827:	612/927=66.02%


## cls4 accuracy ##

H1975:	327/346=94.51
H2228:	325/327=99.39%
HCC827:	612/927=66.02%


## cls5 accuracy ##

H1975:	327/346=94.51%
H2228:	314/327=96.02%
HCC827:	612/927=66.02%


## cls6 accuracy ##

H1975:	326/346=94.22%
H2228:	314/327=96.02%
HCC827:	877/927=94.61%




#################
## 10x_3cl_220 ##
#################


library(CaSpER)
library(biomaRt)
library(GenomicFeatures)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/casper/10x_3cl_220/results/'
output <- './Results/Tian/casper/10x_3cl_220/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_220.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))

## cluster by large scale cnv ##

load(paste0(input,'/large_cnv.rdata'))
l_cnv <- l_cnv[grep('_1',rownames(l_cnv)),]
l_cnv_d <- dist(l_cnv,method='euclidean')
l_cnv_fit <- hclust(l_cnv_d,method='ward.D')
l_cnv_cls <- cutree(l_cnv_fit,3:12,order_clusters_as_data=F)
colnames(l_cnv_cls) <- paste0('pred_cls_',colnames(l_cnv_cls))
cell_cls <- cell_cls[rownames(l_cnv_cls),]
l_cnv_cls_all <- cbind(cell_cls,l_cnv_cls)
write.csv(l_cnv_cls_all,file=paste0(output,'large_cnv_cls.csv'),quote=F)


## cluster by gene based cnv ##

load(paste0(input,'/casper_segment.rdata'))
load(paste0(input,'/casper_final.rdata'))
gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
rna.matrix <- t(rna.matrix)
gene_d <- dist(rna.matrix,method='euclidean')
gene_fit <- hclust(gene_d,method='ward.D')
gene_cls <- cutree(gene_fit,3:12,order_clusters_as_data=F)
colnames(gene_cls) <- paste0('pred_cls_',colnames(gene_cls))
cell_cls <- cell_cls[rownames(gene_cls),]
gene_cls_all <- cbind(cell_cls,gene_cls)
write.csv(gene_cls_all,file=paste0(output,'gene_cls.csv'),quote=F)


## cluster smoothed gene expression ##

obj <- final.objects[[9]]
s_exprs <- t(obj@control.normalized[[3]])
s_exprs <- s_exprs[grep('_1',rownames(s_exprs)),]
s_exprs_d <- dist(s_exprs,method='euclidean')
s_exprs_fit <- hclust(s_exprs_d,method='ward.D')
s_exprs_cls <- cutree(s_exprs_fit,3:12,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]
cls_all <- cbind(cell_cls,s_exprs_cls)
s_exprs_cls3 <- table(cell_cls[,2],s_exprs_cls[,1])
s_exprs_cls4 <- table(cell_cls[,2],s_exprs_cls[,2])
s_exprs_cls5 <- table(cell_cls[,2],s_exprs_cls[,3])
s_exprs_cls6 <- table(cell_cls[,2],s_exprs_cls[,4])
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)
write.csv(s_exprs_cls3,file=paste0(output,'cls3_table.csv'),quote=F)
write.csv(s_exprs_cls4,file=paste0(output,'cls4_table.csv'),quote=F)
write.csv(s_exprs_cls5,file=paste0(output,'cls5_table.csv'),quote=F)
write.csv(s_exprs_cls6,file=paste0(output,'cls6_table.csv'),quote=F)

col <- brewer.pal(8,"Set1")
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_exprs_fit)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=3, value=col[1:3]) %>% plot
colored_bars(colors=bar_col,y_shift=-200,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975'),col=col[1:3],pch=15,bty='n')
dev.off()



#################
## 10x_5cl_310 ##
#################


library(CaSpER)
library(biomaRt)
library(GenomicFeatures)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/casper/10x_5cl_310/results/'
output <- './Results/Tian/casper/10x_5cl_310/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_5cl_310.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,3,4,0,1)))

## cluster by large scale cnv ##

load(paste0(input,'/large_cnv.rdata'))
l_cnv <- l_cnv[grep('_1',rownames(l_cnv)),]
l_cnv_d <- dist(l_cnv,method='euclidean')
l_cnv_fit <- hclust(l_cnv_d,method='ward.D')
l_cnv_cls <- cutree(l_cnv_fit,5:14,order_clusters_as_data=F)
colnames(l_cnv_cls) <- paste0('pred_cls_',colnames(l_cnv_cls))
cell_cls <- cell_cls[rownames(l_cnv_cls),]
l_cnv_cls_all <- cbind(cell_cls,l_cnv_cls)
write.csv(l_cnv_cls_all,file=paste0(output,'large_cnv_cls.csv'),quote=F)


## cluster by gene based cnv ##

load(paste0(input,'/casper_segment.rdata'))
load(paste0(input,'/casper_final.rdata'))
gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
rna.matrix <- t(rna.matrix)
gene_d <- dist(rna.matrix,method='euclidean')
gene_fit <- hclust(gene_d,method='ward.D')
gene_cls <- cutree(gene_fit,5:14,,order_clusters_as_data=F)
colnames(gene_cls) <- paste0('pred_cls_',colnames(gene_cls))
cell_cls <- cell_cls[rownames(gene_cls),]
gene_cls_all <- cbind(cell_cls,gene_cls)
write.csv(gene_cls_all,file=paste0(output,'gene_cls.csv'),quote=F)


## cluster smoothed gene expression ##

obj <- final.objects[[9]]
s_exprs <- t(obj@control.normalized[[3]])
s_exprs <- s_exprs[grep('_1',rownames(s_exprs)),]
s_exprs_d <- dist(s_exprs,method='euclidean')
s_exprs_fit <- hclust(s_exprs_d,method='ward.D')
s_exprs_cls <- cutree(s_exprs_fit,5:14,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]
cls_all <- cbind(cell_cls,s_exprs_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)

col <- brewer.pal(8,"Set1")
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_exprs_fit)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=10, value=col[c(1,5,6,3,2,4,4,7,5,8)]) %>% plot
colored_bars(colors=bar_col,y_shift=-100,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975','A549','H838'),col=col[1:5],pch=15,bty='n')
dev.off()




#################
## 10x_5cl_201 ##
#################


library(CaSpER)
library(biomaRt)
library(GenomicFeatures)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/casper/10x_5cl_201/results/'
output <- './Results/Tian/casper/10x_5cl_201/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_5cl_201.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,3,4,0,1)))

## cluster by large scale cnv ##

load(paste0(input,'/large_cnv.rdata'))
l_cnv <- l_cnv[grep('_1',rownames(l_cnv)),]
l_cnv_d <- dist(l_cnv,method='euclidean')
l_cnv_fit <- hclust(l_cnv_d,method='ward.D')
l_cnv_cls <- cutree(l_cnv_fit,5:14,order_clusters_as_data=F)
colnames(l_cnv_cls) <- paste0('pred_cls_',colnames(l_cnv_cls))
cell_cls <- cell_cls[rownames(l_cnv_cls),]
l_cnv_cls_all <- cbind(cell_cls,l_cnv_cls)
write.csv(l_cnv_cls_all,file=paste0(output,'large_cnv_cls.csv'),quote=F)


## cluster by gene based cnv ##

load(paste0(input,'/casper_segment.rdata'))
load(paste0(input,'/casper_final.rdata'))
gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
rna.matrix <- t(rna.matrix)
gene_d <- dist(rna.matrix,method='euclidean')
gene_fit <- hclust(gene_d,method='ward.D')
gene_cls <- cutree(gene_fit,5:14,,order_clusters_as_data=F)
colnames(gene_cls) <- paste0('pred_cls_',colnames(gene_cls))
cell_cls <- cell_cls[rownames(gene_cls),]
gene_cls_all <- cbind(cell_cls,gene_cls)
write.csv(gene_cls_all,file=paste0(output,'gene_cls.csv'),quote=F)


## cluster smoothed gene expression ##

obj <- final.objects[[9]]
s_exprs <- t(obj@control.normalized[[3]])
s_exprs <- s_exprs[grep('_1',rownames(s_exprs)),]
s_exprs_d <- dist(s_exprs,method='euclidean')
s_exprs_fit <- hclust(s_exprs_d,method='ward.D')
s_exprs_cls <- cutree(s_exprs_fit,5:14,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]
cls_all <- cbind(cell_cls,s_exprs_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)

col <- brewer.pal(8,"Set1")
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_exprs_fit)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=13, value=col[c(4,4,1,4,5,3,2,1,1,5,6,2,7)]) %>% plot
colored_bars(colors=bar_col,y_shift=-100,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975','A549','H838'),col=col[1:5],pch=15,bty='n')
dev.off()



#############
## celseq2 ##
#############

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/casper/celseq2/results/'
output <- './Results/Tian/casper/celseq2/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/celseq2.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))

## cluster by large scale cnv ##

load(paste0(input,'/large_cnv.rdata'))
l_cnv <- l_cnv[grep('_1',rownames(l_cnv)),]
l_cnv_d <- dist(l_cnv,method='euclidean')
l_cnv_fit <- hclust(l_cnv_d,method='ward.D')
l_cnv_cls <- cutree(l_cnv_fit,3:12,order_clusters_as_data=F)
colnames(l_cnv_cls) <- paste0('pred_cls_',colnames(l_cnv_cls))
cell_cls <- cell_cls[rownames(l_cnv_cls),]
l_cnv_cls_all <- cbind(cell_cls,l_cnv_cls)
write.csv(l_cnv_cls_all,file=paste0(output,'large_cnv_cls.csv'),quote=F)


## cluster by gene based cnv ##

load(paste0(input,'/casper_segment.rdata'))
load(paste0(input,'/casper_final.rdata'))
gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
rna.matrix <- t(rna.matrix)
gene_d <- dist(rna.matrix,method='euclidean')
gene_fit <- hclust(gene_d,method='ward.D')
gene_cls <- cutree(gene_fit,3:12,,order_clusters_as_data=F)
colnames(gene_cls) <- paste0('pred_cls_',colnames(gene_cls))
cell_cls <- cell_cls[rownames(gene_cls),]
gene_cls_all <- cbind(cell_cls,gene_cls)
write.csv(gene_cls_all,file=paste0(output,'gene_cls.csv'),quote=F)


## cluster smoothed gene expression ##

obj <- final.objects[[9]]
s_exprs <- t(obj@control.normalized[[3]])
s_exprs <- s_exprs[grep('_1',rownames(s_exprs)),]
s_exprs_d <- dist(s_exprs,method='euclidean')
s_exprs_fit <- hclust(s_exprs_d,method='ward.D')
s_exprs_cls <- cutree(s_exprs_fit,3:12,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]
cls_all <- cbind(cell_cls,s_exprs_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)

col <- brewer.pal(8,"Set1")
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_exprs_fit)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=4, value=col[c(1,4,2,3)]) %>% plot
colored_bars(colors=bar_col,y_shift=-30,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975'),col=col[1:3],pch=15,bty='n')
dev.off()


## cls3 accuracy ##

H1975:	hard to predict
H2228:	74/76=97.37%
HCC827:	hard to predict


## cls4 accuracy ##

H1975:	58/91=63.74%
H2228:	74/76=97.37%
HCC827:	48/75=64%






#############
## dropseq ##
#############

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)
library(RColorBrewer)
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/Tian/casper/dropseq/results/'
output <- './Results/Tian/casper/dropseq/summary/'
dir.create(output,recursive=T)
cell_cls <- read.csv('./Results/Tian/cluster/dropseq.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,1,0)))

## cluster by large scale cnv ##

load(paste0(input,'/large_cnv.rdata'))
l_cnv <- l_cnv[grep('_1',rownames(l_cnv)),]
l_cnv_d <- dist(l_cnv,method='euclidean')
l_cnv_fit <- hclust(l_cnv_d,method='ward.D')
l_cnv_cls <- cutree(l_cnv_fit,3:12,order_clusters_as_data=F)
colnames(l_cnv_cls) <- paste0('pred_cls_',colnames(l_cnv_cls))
cell_cls <- cell_cls[rownames(l_cnv_cls),]
l_cnv_cls_all <- cbind(cell_cls,l_cnv_cls)
write.csv(l_cnv_cls_all,file=paste0(output,'large_cnv_cls.csv'),quote=F)


## cluster by gene based cnv ##

load(paste0(input,'/casper_segment.rdata'))
load(paste0(input,'/casper_final.rdata'))
gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
rna.matrix <- t(rna.matrix)
gene_d <- dist(rna.matrix,method='euclidean')
gene_fit <- hclust(gene_d,method='ward.D')
gene_cls <- cutree(gene_fit,3:12,,order_clusters_as_data=F)
colnames(gene_cls) <- paste0('pred_cls_',colnames(gene_cls))
cell_cls <- cell_cls[rownames(gene_cls),]
gene_cls_all <- cbind(cell_cls,gene_cls)
write.csv(gene_cls_all,file=paste0(output,'gene_cls.csv'),quote=F)



## cluster smoothed gene expression ##

obj <- final.objects[[9]]
s_exprs <- t(obj@control.normalized[[3]])
s_exprs <- s_exprs[grep('_1',rownames(s_exprs)),]
s_exprs_d <- dist(s_exprs,method='euclidean')
s_exprs_fit <- hclust(s_exprs_d,method='ward.D')
s_exprs_cls <- cutree(s_exprs_fit,3:12,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]
cls_all <- cbind(cell_cls,s_exprs_cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)

col <- brewer.pal(8,"Set1")
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_exprs_fit)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=4, value=col[c(3,1,2,4)]) %>% plot
colored_bars(colors=bar_col,y_shift=-30,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975'),col=col[1:3],pch=15,bty='n')
dev.off()



## cls3 accuracy ##

H1975:	21/84=25%
H2228:	54/56=96.43%
HCC827:	61/61=100%


## cls4 accuracy ##

H1975:	81/84=96.43%
H2228:	54/56=96.43%
HCC827:	61/61=100%







