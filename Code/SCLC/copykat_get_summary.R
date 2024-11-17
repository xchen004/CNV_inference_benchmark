
get_sens_spec_copykat <- function(input, gene_anno, cnv_loss, cnv_gain, version)
{

copykat_files <- paste0(input,'rna.matrix.rdata')
load(copykat_files)
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

library(biomaRt)
library(GenomicFeatures)
library(ggplot2)
library(caret)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_primary_loss <- read.csv(paste0(anno_input,'cnv_cyto_sclc_primary_loss_v1.csv'))
cnv_primary_gain <- read.csv(paste0(anno_input,'cnv_cyto_sclc_primary_gain_v1.csv'))
cnv_relapse_loss <- read.csv(paste0(anno_input,'cnv_cyto_sclc_relapse_loss_v1.csv'))
cnv_relapse_gain <- read.csv(paste0(anno_input,'cnv_cyto_sclc_relapse_gain_v1.csv'))

status <- c('Primary','Relapse')
input <- paste0('./Results/SCLC/copykat/',status,'/')

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	if(i ==1 )
	{temp <- get_sens_spec_copykat(input[i], gene_anno, cnv_primary_loss, cnv_primary_gain, 'v1')}
	if(i == 2)
	{temp <- get_sens_spec_copykat(input[i], gene_anno, cnv_relapse_loss, cnv_relapse_gain, 'v1')}
	temp_all <- cbind.data.frame(round(temp[[1]],3),rownames(temp[[1]]),status[i])
	temp_loss <- cbind.data.frame(round(temp[[2]],3),rownames(temp[[2]]),status[i])
	temp_gain <- cbind.data.frame(round(temp[[3]],3),rownames(temp[[3]]),status[i])
	sens_spec_all <- rbind(sens_spec_all,temp_all)
	sens_spec_loss <- rbind(sens_spec_loss,temp_loss)
	sens_spec_gain <- rbind(sens_spec_gain,temp_gain)
}


colnames(sens_spec_all) <- colnames(sens_spec_loss) <- colnames(sens_spec_gain) <- c('Sens','Spec','CNV_cutoff','Status')
output <- './Results/SCLC/copykat/summary/'
dir.create(output,recursive=T)
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_copykat.csv'),quote=F,row.names=F)



## clustering ##


library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)

input <- './Results/SCLC/copykat/'
output <- './Results/SCLC/copykat/summary/'
dir.create(output,recursive=T)

## cluster ##

s_score <- read.table(paste0(input,'results_copykat_CNA_results.txt'),header=T)
s_score2 <- round(s_score[,grep('_1',colnames(s_score))],3)
s_score2 <- s_score2[apply(s_score2,1,sum)!=0,]
cell_cls <- matrix(c(rep('primary',93),rep('relapase',39)),ncol=1)
rownames(cell_cls) <- colnames(s_score2)
colnames(cell_cls) <- 'Status'
s_score_d <- dist(t(s_score2),method='euclidean')
s_score_fit <- hclust(s_score_d,method='ward.D')
save(s_score_fit,file=paste0(output,'s_score_fit.rdata'))
s_score_cls <- cutree(s_score_fit,2:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
cls_all <- cbind(cell_cls,s_score_cls)
colnames(cls_all)[1] <- 'Status'
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)


