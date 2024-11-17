
##############
## infercnv ##
##############

## generate rna.matrix ##

input <- './Results/SCLC/infercnv/'
output <- './Results/SCLC/infercnv/summary/'
da_prob <- read.table(paste0(input,'results/BayesNetOutput.HMMi6.rand_trees.hmm_mode-subclusters/infercnv.NormalProbabilities.observations.txt'),header=T, check.names=F)
da_cnv <- matrix(0,nrow=nrow(da_prob),ncol=ncol(da_prob))
colnames(da_cnv) <- colnames(da_prob)
rownames(da_cnv) <- rownames(da_prob)
da_genes <- read.table(paste0(input,'results/12_HMM_preds.pred_cnv_genes.dat'),header=T)
if(length(table(da_genes[,3])) == 2)
{
da_genes[,3] <- da_genes[,3]-2
}
if(length(table(da_genes[,3])) > 2)
{
da_genes[da_genes[,3]<3,3] <- -1
da_genes[da_genes[,3]>3,3] <- 1
}

da_cell <- read.table(paste0(input,'results/12_HMM_preds.cell_groupings'),header=T)
da_cell <- as.matrix(da_cell[grep('observations',da_cell[,1]),])
group <- names(table(da_cell[,1]))

for(i in 1:length(group))
{
	temp_genes <- da_genes[da_genes[,1]==group[i],]
	temp_cell <- da_cell[da_cell[,1]==group[i],2]
	da_cnv[temp_genes[,4],temp_cell] <- temp_genes[,3]
}

da_cnv3 <- da_cnv*da_prob
da_cnv3[da_cnv3>0] <- 1
da_cnv3[da_cnv3<0] <- -1
#rna.matrix <- as.matrix(da_cnv3)
rna.matrix <- da_cnv

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))
rna.matrix.org <- rna.matrix
rna.matrix <- rna.matrix.org[,grep('Pri',colnames(rna.matrix.org))]
temp_output <- paste0(output,'Primary/')
dir.create(temp_output,recursive=T)
save(rna.matrix,file=paste0(temp_output,'rna.matrix.rdata'))

rna.matrix <- rna.matrix.org[,grep('Re',colnames(rna.matrix.org))]
temp_output <- paste0(output,'Relapse/')
dir.create(temp_output,recursive=T)
save(rna.matrix,file=paste0(temp_output,'rna.matrix.rdata'))


## get sensitivity and specificity ##


get_sens_spec_infercnv <- function(input, gene_anno, cnv_loss, cnv_gain, version)
{

load(paste0(input,'rna.matrix.rdata'))
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


anno_input <- './Anno/cyto/'
gene_anno <- read.csv(paste0(anno_input,'gene_anno_cyto.csv'))
gene_anno[,7] <- paste0(gene_anno[,2],gene_anno[,7])
cnv_primary_loss <- read.csv(paste0(anno_input,'cnv_cyto_sclc_primary_loss_v1.csv'))
cnv_primary_gain <- read.csv(paste0(anno_input,'cnv_cyto_sclc_primary_gain_v1.csv'))
cnv_relapse_loss <- read.csv(paste0(anno_input,'cnv_cyto_sclc_relapse_loss_v1.csv'))
cnv_relapse_gain <- read.csv(paste0(anno_input,'cnv_cyto_sclc_relapse_gain_v1.csv'))

status <- c('Primary','Relapse')
input <- paste0('./Results/SCLC/infercnv/summary/',status,'/')

sens_spec_all <- sens_spec_loss <- sens_spec_gain <- c()
for(i in 1:length(input))
{
	if(i ==1 )
	{temp <- get_sens_spec_infercnv(input[i], gene_anno, cnv_primary_loss, cnv_primary_gain, 'v1')}
	if(i == 2)
	{temp <- get_sens_spec_infercnv(input[i], gene_anno, cnv_relapse_loss, cnv_relapse_gain, 'v1')}
	temp_all <- cbind.data.frame(round(temp[[1]],3),rownames(temp[[1]]),status[i])
	temp_loss <- cbind.data.frame(round(temp[[2]],3),rownames(temp[[2]]),status[i])
	temp_gain <- cbind.data.frame(round(temp[[3]],3),rownames(temp[[3]]),status[i])
	sens_spec_all <- rbind(sens_spec_all,temp_all)
	sens_spec_loss <- rbind(sens_spec_loss,temp_loss)
	sens_spec_gain <- rbind(sens_spec_gain,temp_gain)
}


colnames(sens_spec_all) <- colnames(sens_spec_loss) <- colnames(sens_spec_gain) <- c('Sens','Spec','CNV_cutoff','Status')
output <- './Results/SCLC/infercnv/summary/'
dir.create(output,recursive=T)
write.csv(sens_spec_all,file=paste0(output,'sens_spec_all_infercnv.csv'),quote=F,row.names=F)
write.csv(sens_spec_loss,file=paste0(output,'sens_spec_loss_infercnv.csv'),quote=F,row.names=F)
write.csv(sens_spec_gain,file=paste0(output,'sens_spec_gain_infercnv.csv'),quote=F,row.names=F)







## figure generation ##

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)
library(RColorBrewer)

input <- './Results/SCLC/infercnv/results2/'
output <- './Results/SCLC/infercnv/summary/'
dir.create(output,recursive=T)
cell_cls <- matrix(c(rep('primary',93),rep('relapase',39)),ncol=1)

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,2:12)
rownames(cell_cls) <- rownames(cls)
colnames(cell_cls) <- 'Status'
colnames(cls) <- paste0('pred_cls',colnames(cls))
cls_all <- cbind(cell_cls,cls)
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)


## fig ##

col <- hue_pal()(5)
col2 <- brewer.pal(10,'Set3')

cell_cls[,1] <- factor(cell_cls[,1],levels=c('primary','relapse'))


bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,1])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p2 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'inferCNV') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p2_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col,col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p3_col <- ggplot(bar_col, aes(x=id, y=pos, fill=cell_line)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + labs(fill = 'Cell line') +
scale_fill_manual(values=col) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.position = 'bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p4_col <- ggplot(bar_col2, aes(x=id, y=pos, fill=protocol)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + labs(fill = 'Batch') +
scale_fill_manual(values=col2[1:4]) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), , legend.position = 'bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))


p_legend <- get_legend(p3_col)
p_legend2 <- get_legend(p4_col)
p_legend_all <- plot_grid(p_legend,p_legend2,nrow=1)

p2_col <- p2_col + theme(legend.position='none')

gp2 <- ggplotGrob(p2)
gp2_col <- ggplotGrob(p2_col)  

maxWidth <- grid::unit.pmax(gp2$widths[2:5], gp2_col$widths[2:5])
gp2$widths[2:5] <- as.list(maxWidth)
gp2_col$widths[2:5] <- as.list(maxWidth)

p2_all <- plot_grid(gp2, gp2_col, ncol=1,rel_heights=c(6,1))














