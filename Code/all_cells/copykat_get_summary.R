
get_sens_spec_copykat <- function(input, gene_anno, cnv_loss, cnv_gain, version)
{

copykat_files <- paste0(input,'rna.matrix.rdata')
load(copykat_files)
cell_num <- ncol(rna.matrix)
index <- match(rownames(rna.matrix),gene_anno[,1])
rna.matrix <- rna.matrix[!is.na(index),]
gene_anno <- gene_anno[index[!is.na(index)],]

cutoff <- seq(0,1,by=0.05)
cell_perct <- apply(rna.matrix,1,sum)/cell_num

for(i in 1:length(cutoff))
{
	cnv_genes_loss <- gene_anno[cell_perct < -cutoff[i],]
	cnv_genes_gain <- gene_anno[cell_perct > cutoff[i],]
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
## get summary for all cell ##
#########################################

library(biomaRt)
library(GenomicFeatures)
library(ggplot2)
library(caret)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_loss <- read.csv(paste0(anno_input,'cnv_cyto_brca_loss_v1.csv'))
cnv_gain <- read.csv(paste0(anno_input,'cnv_cyto_brca_gain_v1.csv'))

protocol <- c('10x','C1-HT','C1','ICELL8')
input <- paste0('./Results/',c('10x','c1_fda','c1_llu','WaferGen'),'/all_cell/copykat/v3/')

protocol <- c('10x','C1-HT','C1','ICELL8','Bulk')
input <- paste0('./Results/',c('10x','c1_fda','c1_llu','WaferGen','bulk'),'/all_cell/copykat/v2/')

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_copykat(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
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
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_copykat.csv'),quote=F,row.names=F)




######################################################################
## get summary for subsampling read depth, read length, and subcell ##
######################################################################

library(biomaRt)
library(GenomicFeatures)
library(copykat)
library(ggplot2)
library(caret)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_loss <- read.csv(paste0(anno_input,'cnv_cyto_brca_loss_v1.csv'))
cnv_gain <- read.csv(paste0(anno_input,'cnv_cyto_brca_gain_v1.csv'))

protocol <- c('10x','C1-HT','C1-SE','C1-PE','ICELL8')
type <- c('50bp/50K','75bp/50K','100bp/50K','50bp/100K','75bp/100K','100bp/100K','50bp/250K','75bp/250K','100bp/250K')
for(j in 1:length(type))
{
input <- paste0('./Results/',c(paste0('10x/read_len_subsampling_SE_subCell/cell_1/',type[j],'/copykat/v1_r1/'),paste0('c1_fda/read_len_subsampling_SE_subCell/cell_1/',type[j],'/copykat/v1_r1/'),paste0('c1_llu/read_len_subsampling_SE/',type[j],'/copykat/v1_r1/'),paste0('c1_llu/read_len_subsampling/',type[j],'/copykat/v1_r1/'),paste0('WaferGen/read_len_subsampling_SE_subCell/cell_1/',type[j],'/copykat/v1_r1/')))

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_copykat(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
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
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_copykat.csv'),quote=F,row.names=F)
}


protocol <- c('C1-HT','C1-SE','C1-PE','ICELL8')
type <- c('50bp/500K','75bp/500K','100bp/500K','125bp/50K','125bp/100K','125bp/250K','125bp/500K')
for(j in 1:length(type))
{
input <- paste0('./Results/',c(paste0('c1_fda/read_len_subsampling_SE_subCell/cell_1/',type[j],'/copykat/v1_r1/'),paste0('c1_llu/read_len_subsampling_SE/',type[j],'/copykat/v1_r1/'),paste0('c1_llu/read_len_subsampling/',type[j],'/copykat/v1_r1/'),paste0('WaferGen/read_len_subsampling_SE_subCell/cell_1/',type[j],'/copykat/v1_r1/')))

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_copykat(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
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
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_copykat.csv'),quote=F,row.names=F)
}



protocol <- c('C1-SE','C1-PE','ICELL8')
type <- c('50bp/1M','75bp/1M','100bp/1M','125bp/1M','150bp/50K','150bp/100K','150bp/250K','150bp/500K','150bp/1M')
for(j in 1:length(type))
{
input <- paste0('./Results/',c(paste0('c1_llu/read_len_subsampling_SE/',type[j],'/copykat/v1_r1/'),paste0('c1_llu/read_len_subsampling/',type[j],'/copykat/v1_r1/'),paste0('WaferGen/read_len_subsampling_SE_subCell/cell_1/',type[j],'/copykat/v1_r1/')))

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_copykat(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
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
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_copykat.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_copykat.csv'),quote=F,row.names=F)
}




## get cell percentage for highly recurrent cnvs ##


get_perct_copykat <- function(input, gene_anno, cnv_all, version)
{

copykat_files <- paste0(input,'rna.matrix.rdata')
load(copykat_files)
cell_num <- ncol(rna.matrix)
index <- match(rownames(rna.matrix),gene_anno[,1])
rna.matrix <- rna.matrix[!is.na(index),]
gene_anno <- gene_anno[index[!is.na(index)],]
cell_perct <- apply(rna.matrix,1,sum)/cell_num
perct <- c()
for(i in 1:nrow(cnv_all))
{
	temp_index <- gene_anno[,7]==cnv_all[i,1]
	perct <- c(perct,mean(cell_perct[temp_index]))
}
perct[is.na(perct)] <- 0
cell_perct <- cbind(cnv_all,round(perct,2))
colnames(cell_perct)[3] <- 'Cell_perct'
output <- paste0(input,'sens_spec_summary/')
dir.create(output,recursive=T)
write.csv(cell_perct,file=paste0(output,'cell_perct_',version,'.csv'),quote=F,row.names=F)
return(cell_perct)

}


library(biomaRt)
library(GenomicFeatures)
library(copykat)
library(ggplot2)
library(caret)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_all <- read.csv(paste0(anno_input,'cnv_cyto_brca.csv'))

protocol <- c('10x','C1-HT','C1','ICELL8')
input <- paste0('./Results/',c('10x','c1_fda','c1_llu','WaferGen'),'/all_cell/copykat/v2/')

cell_perct_all <- c()
cell_perct_summary <- c()
cut_off <- c(0.01,0.05,0.1)
for(i in 1:length(input))
{
	temp <- get_perct_copykat(input[i], gene_anno, cnv_all, 'v1')
	temp_all <- cbind.data.frame(temp,protocol[i])
	cell_perct_all <- rbind(cell_perct_all,temp_all)
	temp_perct <- temp_all[,2]*temp_all[,3]
	temp_perct <- temp_perct[temp_perct >= 0]
	temp_perct2 <- c()
	for(j in 1:length(cut_off))
	{
	temp_perct2 <- c(temp_perct2,sum(temp_perct<cut_off[j])/length(temp_perct))
	}
	cell_perct_summary <- rbind(cell_perct_summary,temp_perct2)
}

colnames(cell_perct_all)[4] <- 'Protocols'
colnames(cell_perct_summary) <- paste0('Perct_',cut_off)
rownames(cell_perct_summary) <- protocol
output <- './Results/sens_spec_eva/cyto/all_cell/ref_v2/v1/'
dir.create(output,recursive=T)
write.csv(cell_perct_all,file=paste0(output,'cell_perct_all_copykat.csv'),quote=F,row.names=F)
write.csv(round(cell_perct_summary,2),file=paste0(output,'cell_perct_summary_copykat.csv'),quote=F)


