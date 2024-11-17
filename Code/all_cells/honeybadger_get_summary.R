
get_sens_spec_hb_exprs <- function(input, gene_anno, cnv_loss, cnv_gain, version)
{

gene_anno2 <- read.csv('./Anno/gene_anno_inferCNV.csv')
gene_grange <- GRanges(seqnames=Rle(gene_anno2[,2]),ranges=IRanges(start=gene_anno2[,3],end=gene_anno2[,4]),gene=gene_anno2[,5],ensembl_id=gene_anno2[,1])

load(paste0(input,'t_hb.Rdata'))
rgs <- t_hb$cnvs[["gene-based"]][["all"]]
retest <- t_hb$results[["gene-based"]]
amp.gexp.prob <- do.call(rbind, lapply(retest, function(x) x[[1]]))
del.gexp.prob <- do.call(rbind, lapply(retest, function(x) x[[2]]))

prob_cutoff <- 0.5
amp.gexp.prob[amp.gexp.prob > 0.5] <- 1
amp.gexp.prob[amp.gexp.prob <= 0.5] <- 0
del.gexp.prob[del.gexp.prob <= 0.5] <- 0
del.gexp.prob[del.gexp.prob > 0.5] <- -1
amp_del <- amp.gexp.prob + del.gexp.prob

rna.matrix <- c()
for(i in 1:length(rgs))
{
	temp <- findOverlaps(rgs[i],gene_grange)
	temp_amp_del <- t(matrix(amp_del[i,],nrow=ncol(amp_del),ncol=length(temp)))
	rownames(temp_amp_del) <- gene_grange[subjectHits(temp)]$ensembl_id
	rna.matrix <- rbind(rna.matrix,temp_amp_del)
}
colnames(rna.matrix) <- colnames(amp_del)
cell_num <- ncol(rna.matrix)
index <- match(rownames(rna.matrix),gene_anno[,1])
rna.matrix <- rna.matrix[!is.na(index),]
gene_anno <- gene_anno[index[!is.na(index)],]

cutoff <- seq(0.1,0.8,by=0.1)
cell_perct <- apply(rna.matrix,1,sum)/cell_num

for(i in 1:length(cutoff))
{
	cnv_genes_loss <- gene_anno[cell_perct <= -cutoff[i],]
	cnv_genes_gain <- gene_anno[cell_perct >= cutoff[i],]
	index_loss <- match(cnv_loss[,1],cnv_genes_loss[,7])
	index_gain <- match(cnv_gain[,1],cnv_genes_gain[,7])
	status <- rep('detect',length(index_loss))
	status[is.na(index_loss)] <- 'fail'
	cnv_loss <- cbind(cnv_loss,status)
	status <- rep('detect',length(index_gain))
	status[is.na(index_gain)] <- 'fail'
	cnv_gain <- cbind(cnv_gain,status)
}
colnames(cnv_loss)[-1] <- colnames(cnv_gain)[-1] <- c('wgs',paste0('sc_',cutoff))
cnv_all <- rbind(cnv_loss,cnv_gain)

sens_spec_loss <- sens_spec_gain <- sens_spec_all <- c()
for(i in 1:length(cutoff))
{
	temp_loss <- c(sensitivity(cnv_loss[,(i+2)],cnv_loss[,2]),specificity(cnv_loss[,(i+2)],cnv_loss[,2]))
	sens_spec_loss <-rbind(sens_spec_loss,temp_loss)
	temp_gain <- c(sensitivity(cnv_gain[,(i+2)],cnv_gain[,2]),specificity(cnv_gain[,(i+2)],cnv_gain[,2]))
	sens_spec_gain <-rbind(sens_spec_gain,temp_gain)
	temp_all <- c(sensitivity(cnv_all[,(i+2)],cnv_all[,2]),specificity(cnv_all[,(i+2)],cnv_all[,2]))
	sens_spec_all <-rbind(sens_spec_all,temp_all)
}
colnames(sens_spec_loss) <- colnames(sens_spec_gain) <- colnames(sens_spec_all) <- c('sensitivity','specificity')
rownames(sens_spec_loss) <- rownames(sens_spec_gain) <- rownames(sens_spec_all) <- paste0('sc_',cutoff)

output <- paste0(input,'sens_spec_summary/')
dir.create(output,recursive=T)
write.csv(round(sens_spec_loss,3),file=paste0(output,'sens_spec_loss',version,'.csv'),quote=F)
write.csv(round(sens_spec_gain,3),file=paste0(output,'sens_spec_gain',version,'.csv'),quote=F)
write.csv(round(sens_spec_all,3),file=paste0(output,'sens_spec_all',version,'.csv'),quote=F)
write.csv(cnv_loss,file=paste0(output,'cnv_loss',version,'.csv'),quote=F,row.names=F)
write.csv(cnv_gain,file=paste0(output,'cnv_gain',version,'.csv'),quote=F,row.names=F)
write.csv(cnv_all,file=paste0(output,'cnv_all',version,'.csv'),quote=F,row.names=F)

sens_spec <- vector('list',3)
sens_spec[[1]] <- sens_spec_all
sens_spec[[2]] <- sens_spec_loss
sens_spec[[3]] <- sens_spec_gain
names(sens_spec) <- c('All','Loss','Gain')
return(sens_spec)
}




get_sens_spec_hb_allele <- function(input, gene_anno, cnv_loss, cnv_gain, version)
{

gene_anno2 <- read.csv('./Anno/gene_anno_inferCNV.csv')
gene_grange <- GRanges(seqnames=Rle(gene_anno2[,2]),ranges=IRanges(start=gene_anno2[,3],end=gene_anno2[,4]),gene=gene_anno2[,5],ensembl_id=gene_anno2[,1])

load(paste0(input,'t_hb.Rdata'))
rgs <- t_hb$cnvs[["allele-based"]][["all"]]
retest <- t_hb$results[["allele-based"]]
rgs <- rgs[width(ranges(rgs)) > 1]
del.allele.prob <- do.call(rbind, lapply(retest, function(x) x))

prob_cutoff <- 0.5
del.allele.prob[del.allele.prob <= 0.5] <- 0
del.allele.prob[del.allele.prob > 0.5] <- -1
amp_del <- del.allele.prob

rna.matrix <- c()
for(i in 1:nrow(amp_del))
{
	temp <- findOverlaps(rgs[i],gene_grange)
	temp_amp_del <- t(matrix(amp_del[i,],nrow=ncol(amp_del),ncol=length(temp)))
	rownames(temp_amp_del) <- gene_grange[subjectHits(temp)]$ensembl_id
	rna.matrix <- rbind(rna.matrix,temp_amp_del)
}
colnames(rna.matrix) <- colnames(amp_del)
cell_num <- ncol(rna.matrix)
index <- match(rownames(rna.matrix),gene_anno[,1])
rna.matrix <- rna.matrix[!is.na(index),]
gene_anno <- gene_anno[index[!is.na(index)],]

cutoff <- seq(0.1,0.8,by=0.1)
cell_perct <- apply(rna.matrix,1,sum)/cell_num

for(i in 1:length(cutoff))
{
	cnv_genes_loss <- gene_anno[cell_perct <= -cutoff[i],]
	cnv_genes_gain <- gene_anno[cell_perct >= cutoff[i],]
	index_loss <- match(cnv_loss[,1],cnv_genes_loss[,7])
	index_gain <- match(cnv_gain[,1],cnv_genes_gain[,7])
	status <- rep('detect',length(index_loss))
	status[is.na(index_loss)] <- 'fail'
	cnv_loss <- cbind(cnv_loss,status)
	status <- rep('detect',length(index_gain))
	status[is.na(index_gain)] <- 'fail'
	cnv_gain <- cbind(cnv_gain,status)
}
colnames(cnv_loss)[-1] <- colnames(cnv_gain)[-1] <- c('wgs',paste0('sc_',cutoff))
cnv_all <- rbind(cnv_loss,cnv_gain)

sens_spec_loss <- sens_spec_gain <- sens_spec_all <- c()
for(i in 1:length(cutoff))
{
	temp_loss <- c(sensitivity(cnv_loss[,(i+2)],cnv_loss[,2]),specificity(cnv_loss[,(i+2)],cnv_loss[,2]))
	sens_spec_loss <-rbind(sens_spec_loss,temp_loss)
	temp_gain <- c(sensitivity(cnv_gain[,(i+2)],cnv_gain[,2]),specificity(cnv_gain[,(i+2)],cnv_gain[,2]))
	sens_spec_gain <-rbind(sens_spec_gain,temp_gain)
	temp_all <- c(sensitivity(cnv_all[,(i+2)],cnv_all[,2]),specificity(cnv_all[,(i+2)],cnv_all[,2]))
	sens_spec_all <-rbind(sens_spec_all,temp_all)
}
colnames(sens_spec_loss) <- colnames(sens_spec_gain) <- colnames(sens_spec_all) <- c('sensitivity','specificity')
rownames(sens_spec_loss) <- rownames(sens_spec_gain) <- rownames(sens_spec_all) <- paste0('sc_',cutoff)

output <- paste0(input,'sens_spec_summary/')
dir.create(output,recursive=T)
write.csv(round(sens_spec_loss,3),file=paste0(output,'sens_spec_loss',version,'.csv'),quote=F)
write.csv(round(sens_spec_gain,3),file=paste0(output,'sens_spec_gain',version,'.csv'),quote=F)
write.csv(round(sens_spec_all,3),file=paste0(output,'sens_spec_all',version,'.csv'),quote=F)
write.csv(cnv_loss,file=paste0(output,'cnv_loss',version,'.csv'),quote=F,row.names=F)
write.csv(cnv_gain,file=paste0(output,'cnv_gain',version,'.csv'),quote=F,row.names=F)
write.csv(cnv_all,file=paste0(output,'cnv_all',version,'.csv'),quote=F,row.names=F)

sens_spec <- vector('list',3)
sens_spec[[1]] <- sens_spec_all
sens_spec[[2]] <- sens_spec_loss
sens_spec[[3]] <- sens_spec_gain
names(sens_spec) <- c('All','Loss','Gain')
return(sens_spec)
}



#########################################
## get summary ##
#########################################

## expression ## 

library(biomaRt)
library(GenomicFeatures)
library(ggplot2)
library(caret)
library(HoneyBADGER)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_loss <- read.csv(paste0(anno_input,'cnv_cyto_brca_loss_v1.csv'))
cnv_gain <- read.csv(paste0(anno_input,'cnv_cyto_brca_gain_v1.csv'))

protocol <- c('10x','C1-HT','C1','ICELL8')
input <- paste0('./Results/',c('10x','c1_fda','c1_llu','WaferGen'),'/all_cell/hb/v1/')

protocol <- c('ICELL8')
input <- paste0('./Results/',c('WaferGen'),'/all_cell/hb/v3/')

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_hb_exprs(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
	temp_all <- cbind.data.frame(round(temp[[1]],3),rownames(temp[[1]]),protocol[i])
	temp_loss <- cbind.data.frame(round(temp[[2]],3),rownames(temp[[2]]),protocol[i])
	temp_gain <- cbind.data.frame(round(temp[[3]],3),rownames(temp[[3]]),protocol[i])
	sens_spec_all <- rbind(sens_spec_all,temp_all)
	sens_spec_loss <- rbind(sens_spec_loss,temp_loss)
	sens_spec_gain <- rbind(sens_spec_gain,temp_gain)
}

colnames(sens_spec_all) <- colnames(sens_spec_loss) <- colnames(sens_spec_gain) <- c('Sens','Spec','CNV_cutoff','Protocols')
output <- './Results/sens_spec_eva/cyto/all_cell/ref_v3/v1/'
dir.create(output,recursive=T)
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_hb_exprs.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_hb_exprs.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_hb_exprs.csv'),quote=F,row.names=F)



## allele ## 

library(biomaRt)
library(GenomicFeatures)
library(ggplot2)
library(caret)
library(HoneyBADGER)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_loss <- read.csv(paste0(anno_input,'cnv_cyto_brca_loss_v1.csv'))
cnv_gain <- read.csv(paste0(anno_input,'cnv_cyto_brca_gain_v1.csv'))

protocol <- c('10x','C1','ICELL8')
input <- paste0('./Results/',c('10x','c1_llu','WaferGen'),'/all_cell/hb/v1/')

protocol <- c('C1','ICELL8')
input <- paste0('./Results/',c('c1_llu','WaferGen'),'/all_cell/hb/v3/')

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_hb_allele(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
	temp_all <- cbind.data.frame(round(temp[[1]],3),rownames(temp[[1]]),protocol[i])
	temp_loss <- cbind.data.frame(round(temp[[2]],3),rownames(temp[[2]]),protocol[i])
	temp_gain <- cbind.data.frame(round(temp[[3]],3),rownames(temp[[3]]),protocol[i])
	sens_spec_all <- rbind(sens_spec_all,temp_all)
	sens_spec_loss <- rbind(sens_spec_loss,temp_loss)
	sens_spec_gain <- rbind(sens_spec_gain,temp_gain)
}

colnames(sens_spec_all) <- colnames(sens_spec_loss) <- colnames(sens_spec_gain) <- c('Sens','Spec','CNV_cutoff','Protocols')
output <- './Results/sens_spec_eva/cyto/all_cell/ref_v3/v1/'
dir.create(output,recursive=T)
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_hb_allele.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_hb_allele.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_hb_allele.csv'),quote=F,row.names=F)






######################################################################
## get summary for subsampling read depth, read length, and subcell ##
######################################################################


get_sens_spec_hb_allele <- function(input, gene_anno, cnv_loss, cnv_gain, version)
{

gene_anno2 <- read.csv('./Anno/gene_anno_inferCNV.csv')
gene_grange <- GRanges(seqnames=Rle(gene_anno2[,2]),ranges=IRanges(start=gene_anno2[,3],end=gene_anno2[,4]),gene=gene_anno2[,5],ensembl_id=gene_anno2[,1])

load(paste0(input,'t_hb_',version,'_1.Rdata'))
rgs <- t_hb$cnvs[["allele-based"]][["all"]]
retest <- t_hb$results[["allele-based"]]
rgs <- rgs[width(ranges(rgs)) > 1]
del.allele.prob <- do.call(rbind, lapply(retest, function(x) x))

prob_cutoff <- 0.5
del.allele.prob[del.allele.prob <= 0.5] <- 0
del.allele.prob[del.allele.prob > 0.5] <- -1
amp_del <- del.allele.prob

rna.matrix <- c()
for(i in 1:nrow(amp_del))
{
	temp <- findOverlaps(rgs[i],gene_grange)
	temp_amp_del <- t(matrix(amp_del[i,],nrow=ncol(amp_del),ncol=length(temp)))
	rownames(temp_amp_del) <- gene_grange[subjectHits(temp)]$ensembl_id
	rna.matrix <- rbind(rna.matrix,temp_amp_del)
}
colnames(rna.matrix) <- colnames(amp_del)
cell_num <- ncol(rna.matrix)
index <- match(rownames(rna.matrix),gene_anno[,1])
rna.matrix <- rna.matrix[!is.na(index),]
gene_anno <- gene_anno[index[!is.na(index)],]

cutoff <- seq(0.1,0.8,by=0.1)
cell_perct <- apply(rna.matrix,1,sum)/cell_num

for(i in 1:length(cutoff))
{
	cnv_genes_loss <- gene_anno[cell_perct <= -cutoff[i],]
	cnv_genes_gain <- gene_anno[cell_perct >= cutoff[i],]
	index_loss <- match(cnv_loss[,1],cnv_genes_loss[,7])
	index_gain <- match(cnv_gain[,1],cnv_genes_gain[,7])
	status <- rep('detect',length(index_loss))
	status[is.na(index_loss)] <- 'fail'
	cnv_loss <- cbind(cnv_loss,status)
	status <- rep('detect',length(index_gain))
	status[is.na(index_gain)] <- 'fail'
	cnv_gain <- cbind(cnv_gain,status)
}
colnames(cnv_loss)[-1] <- colnames(cnv_gain)[-1] <- c('wgs',paste0('sc_',cutoff))
cnv_all <- rbind(cnv_loss,cnv_gain)

sens_spec_loss <- sens_spec_gain <- sens_spec_all <- c()
for(i in 1:length(cutoff))
{
	temp_loss <- c(sensitivity(cnv_loss[,(i+2)],cnv_loss[,2]),specificity(cnv_loss[,(i+2)],cnv_loss[,2]))
	sens_spec_loss <-rbind(sens_spec_loss,temp_loss)
	temp_gain <- c(sensitivity(cnv_gain[,(i+2)],cnv_gain[,2]),specificity(cnv_gain[,(i+2)],cnv_gain[,2]))
	sens_spec_gain <-rbind(sens_spec_gain,temp_gain)
	temp_all <- c(sensitivity(cnv_all[,(i+2)],cnv_all[,2]),specificity(cnv_all[,(i+2)],cnv_all[,2]))
	sens_spec_all <-rbind(sens_spec_all,temp_all)
}
colnames(sens_spec_loss) <- colnames(sens_spec_gain) <- colnames(sens_spec_all) <- c('sensitivity','specificity')
rownames(sens_spec_loss) <- rownames(sens_spec_gain) <- rownames(sens_spec_all) <- paste0('sc_',cutoff)

sens_spec <- vector('list',3)
sens_spec[[1]] <- sens_spec_all
sens_spec[[2]] <- sens_spec_loss
sens_spec[[3]] <- sens_spec_gain
names(sens_spec) <- c('All','Loss','Gain')
return(sens_spec)
}



## allele ## 

library(biomaRt)
library(GenomicFeatures)
library(ggplot2)
library(caret)
library(HoneyBADGER)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_loss <- read.csv(paste0(anno_input,'cnv_cyto_brca_loss_v1.csv'))
cnv_gain <- read.csv(paste0(anno_input,'cnv_cyto_brca_gain_v1.csv'))

protocol <- c('C1-PE','C1-SE','ICELL8')
type <- paste0(rep(paste0(seq(50,150, by=25),'bp'),rep(5,5)),'/',rep(c('50K','100K','250K','500K','1M'),5))

for(j in 11:length(type))
{
if(j != c(7,11))
{
input <- paste0('./Results/',c(paste0('c1_llu/read_len_subsampling/',type[j],'/combine/'),paste0('c1_llu/read_len_subsampling_SE/',type[j],'/combine/'),paste0('WaferGen/read_len_subsampling_SE_subCell/cell_1/',type[j],'/combine/')))
}

if(j == 7)
{
input <- paste0('./Results/c1_llu/read_len_subsampling/',type[j],'/combine/')
}

if(j == 11)
{
input <- paste0('./Results/',c(paste0('c1_llu/read_len_subsampling/',type[j],'/combine/'),paste0('WaferGen/read_len_subsampling_SE_subCell/cell_1/',type[j],'/combine/')))
}

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_hb_allele(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
	temp_all <- cbind.data.frame(round(temp[[1]],3),rownames(temp[[1]]),protocol[i])
	temp_loss <- cbind.data.frame(round(temp[[2]],3),rownames(temp[[2]]),protocol[i])
	temp_gain <- cbind.data.frame(round(temp[[3]],3),rownames(temp[[3]]),protocol[i])
	sens_spec_all <- rbind(sens_spec_all,temp_all)
	sens_spec_loss <- rbind(sens_spec_loss,temp_loss)
	sens_spec_gain <- rbind(sens_spec_gain,temp_gain)
}

colnames(sens_spec_all) <- colnames(sens_spec_loss) <- colnames(sens_spec_gain) <- c('Sens','Spec','CNV_cutoff','Protocols')
output <- paste0('./Results/sens_spec_eva/cyto/read_len_subsampling_SE_subCell/cell_1/',type[j],'/ref_v1/v1/')
dir.create(output,recursive=T)
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_hb_allele.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_hb_allele.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_hb_allele.csv'),quote=F,row.names=F)
}





