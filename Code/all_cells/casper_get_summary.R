
get_sens_spec_casper <- function(input, gene_anno, cnv_loss, cnv_gain, version)
{

casper_files <- paste0(input,'rna.matrix.rdata')
load(casper_files)
cell_num <- ncol(rna.matrix)
index <- match(rownames(rna.matrix),gene_anno[,1])
rna.matrix <- rna.matrix[!is.na(index),]
gene_anno <- gene_anno[index[!is.na(index)],]

cutoff <- seq(0,1,by=0.05)
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
## get summary for all cell ##
#########################################

library(biomaRt)
library(GenomicFeatures)
library(CaSpER)
library(ggplot2)
library(caret)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_loss <- read.csv(paste0(anno_input,'cnv_cyto_brca_loss_v1.csv'))
cnv_gain <- read.csv(paste0(anno_input,'cnv_cyto_brca_gain_v1.csv'))

#protocol <- c('10x','C1-HT','C1','ICELL8')
#input <- paste0('./Results/',c('10x','c1_fda','c1_llu','WaferGen'),'/all_cell/casper/v1/')

protocol <- c('10x','C1-HT','C1','ICELL8','Bulk')
input <- paste0('./Results/',c('10x','c1_fda','c1_llu','WaferGen','bulk'),'/all_cell/casper/v3/')

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_casper(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
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
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_casper.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_casper.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_casper.csv'),quote=F,row.names=F)




######################################################################
## get summary for subsampling read depth, read length, and subcell ##
######################################################################

library(biomaRt)
library(GenomicFeatures)
library(CaSpER)
library(ggplot2)
library(caret)

## generate rna.matrix ##

get_matrix <- function(segment, obj, output)
{
segment.summary <- segment
final.objects <- obj
output <- output
gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
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
save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))
return(rna.matrix)
}


## 10x, C1-HT,icell8 read depth evalution ##

input <- './Results/10x/read_len_subsampling_SE_subCell/cell_1/'
rl <- dir(input)
rl <- rl[1:3]
rd <- dir(paste0(input,rl[1]))
rd <- rd[1:3]

input <- './Results/c1_fda/read_len_subsampling_SE_subCell/cell_1/'
rl <- dir(input)
rl <- rl[1:4]
rd <- dir(paste0(input,rl[1]))
rd <- rd[c(1,3:5)]

input <- './Results/sens_spec_eva/length_depth/WaferGen/'
rl <- dir(input)[1:2]
rd <- dir(paste0(input,rl[1]))

for(i in 1:length(rl))
{
	for(j in 1:length(rd))
	{
	output <- paste0(input,rl[i],'/',rd[j],'/casper/v1_r1/')
	load(paste0(output,'/casper_2.rdata'))
	load(paste0(output,'/casper_final_2.rdata'))
	temp <- get_matrix(segment.summary,final.objects,output)
	}
}

## C1, ICELL8 ##

input <- './Results/c1_llu/read_len_subsampling_SE/'
input <- './Results/c1_llu/read_len_subsampling/'
input <- './Results/WaferGen/read_len_subsampling_SE_subCell/cell_1/'

rl <- dir(input)
rl <- rl[1:5]
rd <- dir(paste0(input,rl[1]))
rd <- rd[1:5]

for(i in 1:length(rl))
{
	for(j in 1:length(rd))
	{
	output <- paste0(input,rl[i],'/',rd[j],'/casper/v1_r1/')
	load(paste0(output,'/casper.rdata'))
	load(paste0(output,'/casper_final.rdata'))
	temp <- get_matrix(segment.summary,final.objects,output)
	}
}




## get sensitivity and specificity ##


anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_loss <- read.csv(paste0(anno_input,'cnv_cyto_brca_loss_v1.csv'))
cnv_gain <- read.csv(paste0(anno_input,'cnv_cyto_brca_gain_v1.csv'))

protocol <- c('10x','C1-HT','C1-SE','C1-PE','ICELL8')
type <- c('50bp/50K','75bp/50K','100bp/50K','50bp/100K','75bp/100K','100bp/100K','50bp/250K','75bp/250K','100bp/250K')
for(j in 1:length(type))
{
input <- paste0('./Results/',c(paste0('10x/read_len_subsampling_SE_subCell/cell_1/',type[j],'/casper/v1_r1/'),paste0('c1_fda/read_len_subsampling_SE_subCell/cell_1/',type[j],'/casper/v1_r1/'),paste0('c1_llu/read_len_subsampling_SE/',type[j],'/casper/v1_r1/'),paste0('c1_llu/read_len_subsampling/',type[j],'/casper/v1_r1/'),paste0('WaferGen/read_len_subsampling_SE_subCell/cell_1/',type[j],'/casper/v1_r1/')))

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_casper(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
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
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_casper.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_casper.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_casper.csv'),quote=F,row.names=F)
}


protocol <- c('C1-HT','C1-SE','C1-PE','ICELL8')
type <- c('50bp/500K','75bp/500K','100bp/500K','125bp/50K','125bp/100K','125bp/250K','125bp/500K')
for(j in 1:length(type))
{
input <- paste0('./Results/',c(paste0('c1_fda/read_len_subsampling_SE_subCell/cell_1/',type[j],'/casper/v1_r1/'),paste0('c1_llu/read_len_subsampling_SE/',type[j],'/casper/v1_r1/'),paste0('c1_llu/read_len_subsampling/',type[j],'/casper/v1_r1/'),paste0('WaferGen/read_len_subsampling_SE_subCell/cell_1/',type[j],'/casper/v1_r1/')))

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_casper(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
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
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_casper.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_casper.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_casper.csv'),quote=F,row.names=F)
}



protocol <- c('C1-SE','C1-PE','ICELL8')
type <- c('50bp/1M','75bp/1M','100bp/1M','125bp/1M','150bp/50K','150bp/100K','150bp/250K','150bp/500K','150bp/1M')
for(j in 1:length(type))
{
input <- paste0('./Results/',c(paste0('c1_llu/read_len_subsampling_SE/',type[j],'/casper/v1_r1/'),paste0('c1_llu/read_len_subsampling/',type[j],'/casper/v1_r1/'),paste0('WaferGen/read_len_subsampling_SE_subCell/cell_1/',type[j],'/casper/v1_r1/')))

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	temp <- get_sens_spec_casper(input[i], gene_anno, cnv_loss, cnv_gain, 'v1')
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
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_casper.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_casper.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_casper.csv'),quote=F,row.names=F)
}


## get sensitivity and specificity for icell and c1-ht read depth evaluation ##


library(biomaRt)
library(GenomicFeatures)
library(CaSpER)
library(ggplot2)
library(caret)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_loss <- read.csv(paste0(anno_input,'cnv_cyto_brca_loss_v1.csv'))
cnv_gain <- read.csv(paste0(anno_input,'cnv_cyto_brca_gain_v1.csv'))

protocol <- c('ICELL8','C1-HT')
protocol2 <- c('WaferGen','c1_fda')
input <- './Results/sens_spec_eva/length_depth/'


for(k in 1:length(protocol2))
{
	temp_input <- paste0(input,protocol2[k])
	type <- dir(temp_input)[1:4]
	for(j in 1:length(type))
	{
	rd <- dir(paste0(temp_input,'/',type[j]))
	input2 <- paste0(temp_input,'/',type[j],'/',rd,'/casper/v1_r1/')

	sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
	for(i in 1:length(input2))
	{
	temp <- get_sens_spec_casper(input2[i], gene_anno, cnv_loss, cnv_gain, 'v1')
	temp_all <- cbind.data.frame(round(temp[[1]],3),rownames(temp[[1]]),protocol[k],type[j],rd[i])
	temp_loss <- cbind.data.frame(round(temp[[2]],3),rownames(temp[[2]]),protocol[k],type[j],rd[i])
	temp_gain <- cbind.data.frame(round(temp[[3]],3),rownames(temp[[3]]),protocol[k],type[j],rd[i])
	colnames(temp_all) <- colnames(temp_loss) <- colnames(temp_gain) <- c('Sens','Spec','CNV_cutoff','Protocols','Type','Read_depth')
	sens_spec_all <- rbind(sens_spec_all,temp_all)
	sens_spec_loss <- rbind(sens_spec_loss,temp_loss)
	sens_spec_gain <- rbind(sens_spec_gain,temp_gain)
	}


	output <- paste0('./Results/sens_spec_eva/cyto/read_depth/',protocol[k],'/',type[j],'/ref_v1/v1/')
	dir.create(output,recursive=T)
	write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_casper.csv'),quote=F,row.names=F)
	write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_casper.csv'),quote=F,row.names=F)
	write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_casper.csv'),quote=F,row.names=F)

	}

}



## get cell percentage for highly recurrent cnvs ##


get_perct_casper <- function(input, gene_anno, cnv_all, version)
{

casper_files <- paste0(input,'rna.matrix.rdata')
load(casper_files)
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
library(CaSpER)
library(ggplot2)
library(caret)

anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_all <- read.csv(paste0(anno_input,'cnv_cyto_brca.csv'))

protocol <- c('10x','C1-HT','C1','ICELL8')
input <- paste0('./Results/',c('10x','c1_fda','c1_llu','WaferGen'),'/all_cell/casper/v3/')

cell_perct_all <- c()
cell_perct_summary <- c()
cut_off <- c(0.01,0.05,0.1)
for(i in 1:length(input))
{
	temp <- get_perct_casper(input[i], gene_anno, cnv_all, 'v1')
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
output <- './Results/sens_spec_eva/cyto/all_cell/ref_v3/v1/'
dir.create(output,recursive=T)
write.csv(cell_perct_all,file=paste0(output,'cell_perct_all_casper.csv'),quote=F,row.names=F)
write.csv(round(cell_perct_summary,2),file=paste0(output,'cell_perct_summary_casper.csv'),quote=F)













