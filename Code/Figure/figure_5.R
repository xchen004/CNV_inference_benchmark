
library(ggpubr)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape)
library(cowplot)
library(ggdendro)

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)
library(RColorBrewer)
library(scales)


col_type <- brewer.pal(9,'Set1')

## ARI ##

## version 1: overall ari ##

ari_infercnv <- read.csv('./Results/Tian/cross_protocols/infercnv/summary/adj_rand_index.csv',row.name=1)
ari_infercnv <- melt(ari_infercnv[1:8,1:2])
ari_infercnv <- cbind(rep(3:10,2),ari_infercnv,'inferCNV')
colnames(ari_infercnv) <- c('Cluster','ARI_type','ARI','Methods')
ari_casper <- read.csv('./Results/Tian/cross_protocols/casper/summary/dropseq_baf/adj_rand_index.csv',row.name=1)
ari_casper <- melt(ari_casper[1:8,1:2])
ari_casper <- cbind(rep(3:10,2),ari_casper,'CaSpER')
colnames(ari_casper) <- c('Cluster','ARI_type','ARI','Methods')
ari_scicnv <- read.csv('./Results/Tian/cross_protocols/scicnv/summary/adj_rand_index.csv',row.name=1)
ari_scicnv <- melt(ari_scicnv[1:8,1:2])
ari_scicnv <- cbind(rep(3:10,2),ari_scicnv,'sciCNV')
colnames(ari_scicnv) <- c('Cluster','ARI_type','ARI','Methods')
ari_copykat <- read.csv('./Results/Tian/cross_protocols/copykat/summary/adj_rand_index.csv',row.name=1)
ari_copykat <- melt(ari_copykat[1:8,1:2])
ari_copykat <- cbind(rep(3:10,2),ari_copykat,'CopyKAT')
colnames(ari_copykat) <- c('Cluster','ARI_type','ARI','Methods')
ari_hb_exprs <- read.csv('./Results/Tian/cross_protocols/hb/exprs_adj_rand_index.csv',row.name=1)
ari_hb_exprs <- melt(ari_hb_exprs[1:8,1:2])
ari_hb_exprs <- cbind(rep(3:10,2),ari_hb_exprs,'HoneyBADGER-expression')
colnames(ari_hb_exprs) <- c('Cluster','ARI_type','ARI','Methods')
ari_hb_allele <- read.csv('./Results/Tian/cross_protocols/hb/allele_adj_rand_index.csv',row.name=1)
ari_hb_allele <- melt(ari_hb_allele[1:8,1:2])
ari_hb_allele <- cbind(rep(3:10,2),ari_hb_allele,'HoneyBADGER-allele')
colnames(ari_hb_allele) <- c('Cluster','ARI_type','ARI','Methods')

ari <- rbind(ari_infercnv,ari_casper,ari_scicnv,ari_copykat,ari_hb_exprs,ari_hb_allele)
levels(ari[,2])[levels(ari[,2]) == 'celltype'] <- 'Cell line'
levels(ari[,2])[levels(ari[,2]) == 'protocol'] <- 'Batch'

p1 <- ggplot(ari,aes(x = Cluster, y = ARI, color=ARI_type, shape=Methods, linetype=Methods)) + geom_line() + geom_point(size=2) + 
labs(x = 'Number of clusters', color = 'ARI', shape = 'CNV methods', linetype = 'CNV methods') + 
scale_color_manual(values=col_type) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(3,10), breaks = seq(3,10,by=1)) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times'), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p1_legend <- get_legend(p1)
p1 <- p1+ theme(legend.position='none')


## version 2: ari within celltype and protocol ##

ari_infercnv <- read.csv('./Results/Tian/cross_protocols/infercnv/summary/adj_rand_index.csv',row.name=1)
ari_infercnv <- cbind.data.frame(apply(ari_infercnv[,c(4,5,7)],1,mean),apply(ari_infercnv[,8:11],1,mean))
colnames(ari_infercnv) <- c('protocol','celltype')
ari_infercnv <- melt(ari_infercnv[1:8,1:2])
ari_infercnv <- cbind(rep(3:10,2),ari_infercnv,'inferCNV')
colnames(ari_infercnv) <- c('Cluster','ARI_type','ARI','Methods')
ari_casper <- read.csv('./Results/Tian/cross_protocols/casper/summary/dropseq_baf/adj_rand_index.csv',row.name=1)
ari_casper <- cbind.data.frame(apply(ari_casper[,c(4,5,7)],1,mean),apply(ari_casper[,8:11],1,mean))
colnames(ari_casper) <- c('protocol','celltype')
ari_casper <- melt(ari_casper[1:8,1:2])
ari_casper <- cbind(rep(3:10,2),ari_casper,'CaSpER')
colnames(ari_casper) <- c('Cluster','ARI_type','ARI','Methods')
ari_scicnv <- read.csv('./Results/Tian/cross_protocols/scicnv/summary/adj_rand_index.csv',row.name=1)
ari_scicnv <- cbind.data.frame(apply(ari_scicnv[,c(4,5,7)],1,mean),apply(ari_scicnv[,8:11],1,mean))
colnames(ari_scicnv) <- c('protocol','celltype')
ari_scicnv <- melt(ari_scicnv[1:8,1:2])
ari_scicnv <- cbind(rep(3:10,2),ari_scicnv,'sciCNV')
colnames(ari_scicnv) <- c('Cluster','ARI_type','ARI','Methods')
ari_copykat <- read.csv('./Results/Tian/cross_protocols/copykat/summary/adj_rand_index.csv',row.name=1)
ari_copykat <- cbind.data.frame(apply(ari_copykat[,c(4,5,7)],1,mean),apply(ari_copykat[,8:11],1,mean))
colnames(ari_copykat) <- c('protocol','celltype')
ari_copykat <- melt(ari_copykat[1:8,1:2])
ari_copykat <- cbind(rep(3:10,2),ari_copykat,'CopyKAT')
colnames(ari_copykat) <- c('Cluster','ARI_type','ARI','Methods')
ari_hb_exprs <- read.csv('./Results/Tian/cross_protocols/hb/exprs_adj_rand_index.csv',row.name=1)
ari_hb_exprs <- cbind.data.frame(apply(ari_hb_exprs[,c(4,5,7)],1,mean),apply(ari_hb_exprs[,8:11],1,mean))
colnames(ari_hb_exprs) <- c('protocol','celltype')
ari_hb_exprs <- melt(ari_hb_exprs[1:8,1:2])
ari_hb_exprs <- cbind(rep(3:10,2),ari_hb_exprs,'HoneyBADGER-expression')
colnames(ari_hb_exprs) <- c('Cluster','ARI_type','ARI','Methods')
ari_hb_allele <- read.csv('./Results/Tian/cross_protocols/hb/allele_adj_rand_index.csv',row.name=1)
ari_hb_allele <- cbind.data.frame(apply(ari_hb_allele[,c(4,5,7)],1,mean),apply(ari_hb_allele[,8:11],1,mean))
colnames(ari_hb_allele) <- c('protocol','celltype')
ari_hb_allele <- melt(ari_hb_allele[1:8,1:2])
ari_hb_allele <- cbind(rep(3:10,2),ari_hb_allele,'HoneyBADGER-allele')
colnames(ari_hb_allele) <- c('Cluster','ARI_type','ARI','Methods')

ari <- rbind(ari_infercnv,ari_casper,ari_scicnv,ari_copykat,ari_hb_exprs,ari_hb_allele)
levels(ari[,2])[levels(ari[,2]) == 'celltype'] <- 'Cell line'
levels(ari[,2])[levels(ari[,2]) == 'protocol'] <- 'Batch'
ari[,2] <- factor(ari[,2], levels=c('Cell line','Batch'))

p1_v2 <- ggplot(ari,aes(x = Cluster, y = ARI, color=ARI_type, shape=Methods, linetype=Methods)) + geom_line() + geom_point(size=2) + 
labs(title = 'Before batch correction', x = 'Number of clusters', color = 'ARI', shape = 'CNV methods', linetype = 'CNV methods') +
scale_color_manual(values=col_type) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(3,10), breaks = seq(3,10,by=1)) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times',hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))




## version 3: ari within celltype and protocol after limma ##

ari_infercnv <- read.csv('./Results/Tian/cross_protocols/infercnv/limma/summary/adj_rand_index.csv',row.name=1)
ari_infercnv <- cbind.data.frame(apply(ari_infercnv[,c(4,5,7)],1,mean),apply(ari_infercnv[,8:11],1,mean))
colnames(ari_infercnv) <- c('protocol','celltype')
ari_infercnv <- melt(ari_infercnv[1:8,1:2])
ari_infercnv <- cbind(rep(3:10,2),ari_infercnv,'inferCNV')
colnames(ari_infercnv) <- c('Cluster','ARI_type','ARI','Methods')
ari_casper <- read.csv('./Results/Tian/cross_protocols/casper/limma/summary/dropseq_baf/adj_rand_index.csv',row.name=1)
ari_casper <- cbind.data.frame(apply(ari_casper[,c(4,5,7)],1,mean),apply(ari_casper[,8:11],1,mean))
colnames(ari_casper) <- c('protocol','celltype')
ari_casper <- melt(ari_casper[1:8,1:2])
ari_casper <- cbind(rep(3:10,2),ari_casper,'CaSpER')
colnames(ari_casper) <- c('Cluster','ARI_type','ARI','Methods')
ari_scicnv <- read.csv('./Results/Tian/cross_protocols/scicnv/limma/summary/adj_rand_index.csv',row.name=1)
ari_scicnv <- cbind.data.frame(apply(ari_scicnv[,c(4,5,7)],1,mean),apply(ari_scicnv[,8:11],1,mean))
colnames(ari_scicnv) <- c('protocol','celltype')
ari_scicnv <- melt(ari_scicnv[1:8,1:2])
ari_scicnv <- cbind(rep(3:10,2),ari_scicnv,'sciCNV')
colnames(ari_scicnv) <- c('Cluster','ARI_type','ARI','Methods')
ari_copykat <- read.csv('./Results/Tian/cross_protocols/copykat/limma/summary/adj_rand_index.csv',row.name=1)
ari_copykat <- cbind.data.frame(apply(ari_copykat[,c(4,5,7)],1,mean),apply(ari_copykat[,8:11],1,mean))
colnames(ari_copykat) <- c('protocol','celltype')
ari_copykat <- melt(ari_copykat[1:8,1:2])
ari_copykat <- cbind(rep(3:10,2),ari_copykat,'CopyKAT')
colnames(ari_copykat) <- c('Cluster','ARI_type','ARI','Methods')
ari_hb_exprs <- read.csv('./Results/Tian/cross_protocols/hb/limma/exprs_adj_rand_index.csv',row.name=1)
ari_hb_exprs <- cbind.data.frame(apply(ari_hb_exprs[,c(4,5,7)],1,mean),apply(ari_hb_exprs[,8:11],1,mean))
colnames(ari_hb_exprs) <- c('protocol','celltype')
ari_hb_exprs <- melt(ari_hb_exprs[1:8,1:2])
ari_hb_exprs <- cbind(rep(3:10,2),ari_hb_exprs,'HoneyBADGER-expression')
colnames(ari_hb_exprs) <- c('Cluster','ARI_type','ARI','Methods')

ari <- rbind(ari_infercnv,ari_casper,ari_scicnv,ari_copykat,ari_hb_exprs)
levels(ari[,2])[levels(ari[,2]) == 'celltype'] <- 'Cell line'
levels(ari[,2])[levels(ari[,2]) == 'protocol'] <- 'Batch'
ari[,2] <- factor(ari[,2], levels=c('Cell line','Batch'))

p1_v3 <- ggplot(ari,aes(x = Cluster, y = ARI, color=ARI_type, shape=Methods, linetype=Methods)) + geom_line() + geom_point(size=2) + 
labs(title = 'limma batch correction',x = 'Number of clusters', color = 'ARI', shape = 'CNV methods', linetype = 'CNV methods') +
scale_color_manual(values=col_type) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(3,10), breaks = seq(3,10,by=1)) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times',hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p1_v3 <- p1_v3 + theme(legend.position = 'none')



## version 4: ari within celltype and protocol after combat ##

ari_infercnv <- read.csv('./Results/Tian/cross_protocols/infercnv/combat/summary/adj_rand_index.csv',row.name=1)
ari_infercnv <- cbind.data.frame(apply(ari_infercnv[,c(4,5,7)],1,mean),apply(ari_infercnv[,8:11],1,mean))
colnames(ari_infercnv) <- c('protocol','celltype')
ari_infercnv <- melt(ari_infercnv[1:8,1:2])
ari_infercnv <- cbind(rep(3:10,2),ari_infercnv,'inferCNV')
colnames(ari_infercnv) <- c('Cluster','ARI_type','ARI','Methods')
ari_casper <- read.csv('./Results/Tian/cross_protocols/casper/combat/summary/dropseq_baf/adj_rand_index.csv',row.name=1)
ari_casper <- cbind.data.frame(apply(ari_casper[,c(4,5,7)],1,mean),apply(ari_casper[,8:11],1,mean))
colnames(ari_casper) <- c('protocol','celltype')
ari_casper <- melt(ari_casper[1:8,1:2])
ari_casper <- cbind(rep(3:10,2),ari_casper,'CaSpER')
colnames(ari_casper) <- c('Cluster','ARI_type','ARI','Methods')
ari_scicnv <- read.csv('./Results/Tian/cross_protocols/scicnv/combat/summary/adj_rand_index.csv',row.name=1)
ari_scicnv <- cbind.data.frame(apply(ari_scicnv[,c(4,5,7)],1,mean),apply(ari_scicnv[,8:11],1,mean))
colnames(ari_scicnv) <- c('protocol','celltype')
ari_scicnv <- melt(ari_scicnv[1:8,1:2])
ari_scicnv <- cbind(rep(3:10,2),ari_scicnv,'scicnv')
colnames(ari_scicnv) <- c('Cluster','ARI_type','ARI','Methods')
ari_copykat <- read.csv('./Results/Tian/cross_protocols/copykat/combat/summary/adj_rand_index.csv',row.name=1)
ari_copykat <- cbind.data.frame(apply(ari_copykat[,c(4,5,7)],1,mean),apply(ari_copykat[,8:11],1,mean))
colnames(ari_copykat) <- c('protocol','celltype')
ari_copykat <- melt(ari_copykat[1:8,1:2])
ari_copykat <- cbind(rep(3:10,2),ari_copykat,'CopyKAT')
colnames(ari_copykat) <- c('Cluster','ARI_type','ARI','Methods')
ari_hb_exprs <- read.csv('./Results/Tian/cross_protocols/hb/combat/exprs_adj_rand_index.csv',row.name=1)
ari_hb_exprs <- cbind.data.frame(apply(ari_hb_exprs[,c(4,5,7)],1,mean),apply(ari_hb_exprs[,8:11],1,mean))
colnames(ari_hb_exprs) <- c('protocol','celltype')
ari_hb_exprs <- melt(ari_hb_exprs[1:8,1:2])
ari_hb_exprs <- cbind(rep(3:10,2),ari_hb_exprs,'HoneyBADGER-expression')
colnames(ari_hb_exprs) <- c('Cluster','ARI_type','ARI','Methods')

ari <- rbind(ari_infercnv,ari_casper,ari_scicnv,ari_copykat,ari_hb_exprs)
levels(ari[,2])[levels(ari[,2]) == 'celltype'] <- 'Cell line'
levels(ari[,2])[levels(ari[,2]) == 'protocol'] <- 'Batch'
ari[,2] <- factor(ari[,2], levels=c('Cell line','Batch'))

p1_v4 <- ggplot(ari,aes(x = Cluster, y = ARI, color=ARI_type, shape=Methods, linetype=Methods)) + geom_line() + geom_point(size=2) + 
labs(title = 'ComBat batch correction', x = 'Number of clusters', color = 'ARI', shape = 'CNV methods', linetype = 'CNV methods') +
scale_color_manual(values=col_type) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(3,10), breaks = seq(3,10,by=1)) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times',hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p1_v4 <- p1_v4 + theme(legend.position = 'none')




##############
## inferCNV ##
##############

cell_cls <- read.csv('./Results/Tian/cross_protocols/data/cls.csv',row.names=1)
protocol_id <- rep(1,nrow(cell_cls))
protocol_id[cell_cls[,3] == 'Cel-Seq2'] <- 2
protocol_id[cell_cls[,3] == '10x_3cl'] <- 3
protocol_id[cell_cls[,3] == '10x_5cl'] <- 4
cell_cls <- cbind.data.frame(cell_cls,protocol_id)
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,0,2,3,4)))
cell_cls[,2] <- factor(cell_cls[,2],levels=c('H2228','H1975','HCC827','A549','H838'))
cell_cls[,3] <- factor(cell_cls[,3],levels=c('10x_5cl','10x_3cl','Cel-Seq2','Drop-Seq'))


col <- brewer.pal(9,'Set1')
col2 <- brewer.pal(8,'Set2')


## no batch correction ##

input <- './Results/Tian/cross_protocols/infercnv/results/1_core/'
output <- './Results/Tian/cross_protocols/infercnv/summary/'

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p2 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'inferCNV') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p2_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
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



## limma batch correction ##

input <- './Results/Tian/cross_protocols/infercnv/limma/results/'
output <- './Results/Tian/cross_protocols/infercnv/limma/summary/'

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p2 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'inferCNV after limma') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p2_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p2_col <- p2_col + theme(legend.position='none')

gp2 <- ggplotGrob(p2)
gp2_col <- ggplotGrob(p2_col)  

maxWidth <- grid::unit.pmax(gp2$widths[2:5], gp2_col$widths[2:5])
gp2$widths[2:5] <- as.list(maxWidth)
gp2_col$widths[2:5] <- as.list(maxWidth)

p2_all_limma <- plot_grid(gp2, gp2_col, ncol=1,rel_heights=c(6,1))



## combat batch correction ##

input <- './Results/Tian/cross_protocols/infercnv/combat/results/'
output <- './Results/Tian/cross_protocols/infercnv/combat/summary/'

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p2 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'inferCNV after ComBat') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p2_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p2_col <- p2_col + theme(legend.position='none')

gp2 <- ggplotGrob(p2)
gp2_col <- ggplotGrob(p2_col)  

maxWidth <- grid::unit.pmax(gp2$widths[2:5], gp2_col$widths[2:5])
gp2$widths[2:5] <- as.list(maxWidth)
gp2_col$widths[2:5] <- as.list(maxWidth)

p2_all_combat <- plot_grid(gp2, gp2_col, ncol=1,rel_heights=c(6,1))



############
## casper ##
############

## no batch correction ##

input <- './Results/Tian/cross_protocols/casper/results/dropseq_baf/'
load(paste0(input,'s_exprs_fit.rdata'))
dend <- as.dendrogram(s_exprs_fit)
s_exprs_cls <- cutree(s_exprs_fit,3:12,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p3 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'CaSpER') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p3_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp3 <- ggplotGrob(p3)
gp3_col <- ggplotGrob(p3_col)  

maxWidth <- grid::unit.pmax(gp3$widths[2:5], gp3_col$widths[2:5])
gp3$widths[2:5] <- as.list(maxWidth)
gp3_col$widths[2:5] <- as.list(maxWidth)

p3_all <- plot_grid(gp3, gp3_col, ncol=1,rel_heights=c(6,1))



## limma batch correction ##

input <- './Results/Tian/cross_protocols/casper/limma/results/dropseq_baf/'
load(paste0(input,'s_exprs_fit.rdata'))
dend <- as.dendrogram(s_exprs_fit)
s_exprs_cls <- cutree(s_exprs_fit,3:12,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p3 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'CaSpER after limma') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p3_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp3 <- ggplotGrob(p3)
gp3_col <- ggplotGrob(p3_col)  

maxWidth <- grid::unit.pmax(gp3$widths[2:5], gp3_col$widths[2:5])
gp3$widths[2:5] <- as.list(maxWidth)
gp3_col$widths[2:5] <- as.list(maxWidth)

p3_all_limma <- plot_grid(gp3, gp3_col, ncol=1,rel_heights=c(6,1))



## combat batch correction ##

input <- './Results/Tian/cross_protocols/casper/combat/results/dropseq_baf/'
load(paste0(input,'s_exprs_fit.rdata'))
dend <- as.dendrogram(s_exprs_fit)
s_exprs_cls <- cutree(s_exprs_fit,3:12,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
cell_cls <- cell_cls[rownames(s_exprs_cls),]

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p3 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'CaSpER after ComBat') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p3_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp3 <- ggplotGrob(p3)
gp3_col <- ggplotGrob(p3_col)  

maxWidth <- grid::unit.pmax(gp3$widths[2:5], gp3_col$widths[2:5])
gp3$widths[2:5] <- as.list(maxWidth)
gp3_col$widths[2:5] <- as.list(maxWidth)

p3_all_combat <- plot_grid(gp3, gp3_col, ncol=1,rel_heights=c(6,1))



########
## hb ##
########


library(HoneyBADGER)

## before batch correction ##

load('./Results/Tian/cross_protocols/hb/hb.Rdata')
results_exprs <- t_hb$summarizeResults(geneBased=TRUE, alleleBased=FALSE)
tree_exprs <- t_hb$visualizeResults(geneBased=TRUE, alleleBased=FALSE, details=TRUE, margins=c(25,15))
cls_exprs <- cutree(tree_exprs$hc,k=3:12,order_clusters_as_data=F)
colnames(cls_exprs) <- paste0('pred_cls_',colnames(cls_exprs))
cell_cls <- cell_cls[rownames(cls_exprs),]
dend_exprs <- as.dendrogram(tree_exprs$hc)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p4 <- dend_exprs %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'HoneyBADGER-expression') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p4_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp4 <- ggplotGrob(p4)
gp4_col <- ggplotGrob(p4_col)  
maxWidth <- grid::unit.pmax(gp4$widths[2:5], gp4_col$widths[2:5])
gp4$widths[2:5] <- as.list(maxWidth)
gp4_col$widths[2:5] <- as.list(maxWidth)

p4_all <- plot_grid(gp4, gp4_col, ncol=1,rel_heights=c(6,1))


## limma batch correction ##

load('./Results/Tian/cross_protocols/hb/limma/hb.Rdata')

results_exprs <- t_hb$summarizeResults(geneBased=TRUE, alleleBased=FALSE)
tree_exprs <- t_hb$visualizeResults(geneBased=TRUE, alleleBased=FALSE, details=TRUE, margins=c(25,15))
cls_exprs <- cutree(tree_exprs$hc,k=3:12,order_clusters_as_data=F)
colnames(cls_exprs) <- paste0('pred_cls_',colnames(cls_exprs))
cell_cls <- cell_cls[rownames(cls_exprs),]
dend_exprs <- as.dendrogram(tree_exprs$hc)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p4 <- dend_exprs %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'HoneyBADGER-expression after limma') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p4_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp4 <- ggplotGrob(p4)
gp4_col <- ggplotGrob(p4_col)  
maxWidth <- grid::unit.pmax(gp4$widths[2:5], gp4_col$widths[2:5])
gp4$widths[2:5] <- as.list(maxWidth)
gp4_col$widths[2:5] <- as.list(maxWidth)

p4_all_limma <- plot_grid(gp4, gp4_col, ncol=1,rel_heights=c(6,1))



## combat batch correction ##

load('./Results/Tian/cross_protocols/hb/combat/hb.Rdata')

results_exprs <- t_hb$summarizeResults(geneBased=TRUE, alleleBased=FALSE)
tree_exprs <- t_hb$visualizeResults(geneBased=TRUE, alleleBased=FALSE, details=TRUE, margins=c(25,15))
cls_exprs <- cutree(tree_exprs$hc,k=3:12,order_clusters_as_data=F)
colnames(cls_exprs) <- paste0('pred_cls_',colnames(cls_exprs))
cell_cls <- cell_cls[rownames(cls_exprs),]
dend_exprs <- as.dendrogram(tree_exprs$hc)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p4 <- dend_exprs %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'HoneyBADGER-expression after ComBat') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p4_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp4 <- ggplotGrob(p4)
gp4_col <- ggplotGrob(p4_col)  
maxWidth <- grid::unit.pmax(gp4$widths[2:5], gp4_col$widths[2:5])
gp4$widths[2:5] <- as.list(maxWidth)
gp4_col$widths[2:5] <- as.list(maxWidth)

p4_all_combat <- plot_grid(gp4, gp4_col, ncol=1,rel_heights=c(6,1))



## allele based ##

load('./Results/Tian/cross_protocols/hb/hb.Rdata')
results_allele <- t_hb$summarizeResults(geneBased=FALSE, alleleBased=TRUE)
tree_allele <- t_hb$visualizeResults(geneBased=FALSE, alleleBased=TRUE, details=TRUE, margins=c(25,15))
cls_allele <- cutree(tree_allele$hc,k=3:12,order_clusters_as_data=F)
colnames(cls_allele) <- paste0('pred_cls_',colnames(cls_allele))
cell_cls <- cell_cls[rownames(cls_allele),]
dend_allele <- as.dendrogram(tree_allele$hc)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p5 <- dend_allele %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'HoneyBADGER-allele') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p5_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp5 <- ggplotGrob(p5)
gp5_col <- ggplotGrob(p5_col)  

maxWidth <- grid::unit.pmax(gp5$widths[2:5], gp5_col$widths[2:5])
gp5$widths[2:5] <- as.list(maxWidth)
gp5_col$widths[2:5] <- as.list(maxWidth)

p5_all <- plot_grid(gp5, gp5_col, ncol=1,rel_heights=c(6,1))



############
## sciCNV ##
############

## before batch correction ##

input <- './Results/Tian/cross_protocols/scicnv/results/'
input2 <- './Results/Tian/cross_protocols/infercnv/'
output <- './Results/Tian/cross_protocols/scicnv/summary/'

load(paste0(input,'s_score_fit.rdata'))
cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
cls3_table <- table(cell_cls[,2],cls[,1])
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_score_fit)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p6 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'sciCNV') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p6_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp6 <- ggplotGrob(p6)
gp6_col <- ggplotGrob(p6_col)  

maxWidth <- grid::unit.pmax(gp6$widths[2:5], gp6_col$widths[2:5])
gp6$widths[2:5] <- as.list(maxWidth)
gp6_col$widths[2:5] <- as.list(maxWidth)

p6_all <- plot_grid(gp6, gp6_col, ncol=1,rel_heights=c(6,1))



## limma batch correction ##

input <- './Results/Tian/cross_protocols/scicnv/limma/results/'
input2 <- './Results/Tian/cross_protocols/infercnv/limma/'
output <- './Results/Tian/cross_protocols/scicnv/limma/summary/'

load(paste0(input,'s_score_fit.rdata'))
cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
cls3_table <- table(cell_cls[,2],cls[,1])
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_score_fit)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p6 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'sciCNV after limma') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p6_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp6 <- ggplotGrob(p6)
gp6_col <- ggplotGrob(p6_col)  

maxWidth <- grid::unit.pmax(gp6$widths[2:5], gp6_col$widths[2:5])
gp6$widths[2:5] <- as.list(maxWidth)
gp6_col$widths[2:5] <- as.list(maxWidth)

p6_all_limma <- plot_grid(gp6, gp6_col, ncol=1,rel_heights=c(6,1))


## combat batch correction ##

input <- './Results/Tian/cross_protocols/scicnv/combat/results/'
input2 <- './Results/Tian/cross_protocols/infercnv/combat/'
output <- './Results/Tian/cross_protocols/scicnv/combat/summary/'

load(paste0(input,'s_score_fit.rdata'))
cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
cls3_table <- table(cell_cls[,2],cls[,1])
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_score_fit)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p6 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'sciCNV after ComBat') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p6_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp6 <- ggplotGrob(p6)
gp6_col <- ggplotGrob(p6_col)  

maxWidth <- grid::unit.pmax(gp6$widths[2:5], gp6_col$widths[2:5])
gp6$widths[2:5] <- as.list(maxWidth)
gp6_col$widths[2:5] <- as.list(maxWidth)

p6_all_combat <- plot_grid(gp6, gp6_col, ncol=1,rel_heights=c(6,1))





############
## copykat ##
############

## before batch correction ##

input <- './Results/Tian/cross_protocols/copykat/results/'
output <- './Results/Tian/cross_protocols/copykat/summary/'

load(paste0(input,'s_score_fit.rdata'))
cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
cls3_table <- table(cell_cls[,2],cls[,1])
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_score_fit)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p7 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'CopyKAT') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p7_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp7 <- ggplotGrob(p7)
gp7_col <- ggplotGrob(p7_col)  

maxWidth <- grid::unit.pmax(gp7$widths[2:5], gp7_col$widths[2:5])
gp7$widths[2:5] <- as.list(maxWidth)
gp7_col$widths[2:5] <- as.list(maxWidth)

p7_all <- plot_grid(gp7, gp7_col, ncol=1,rel_heights=c(6,1))



## limma batch correction ##

input <- './Results/Tian/cross_protocols/copykat/limma/results/'
output <- './Results/Tian/cross_protocols/copykat/limma/summary/'

load(paste0(input,'s_score_fit.rdata'))
cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
cls3_table <- table(cell_cls[,2],cls[,1])
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_score_fit)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p7 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'CopyKAT after limma') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p7_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp7 <- ggplotGrob(p7)
gp7_col <- ggplotGrob(p7_col)  

maxWidth <- grid::unit.pmax(gp7$widths[2:5], gp7_col$widths[2:5])
gp7$widths[2:5] <- as.list(maxWidth)
gp7_col$widths[2:5] <- as.list(maxWidth)

p7_all_limma <- plot_grid(gp7, gp7_col, ncol=1,rel_heights=c(6,1))


## combat batch correction ##


input <- './Results/Tian/cross_protocols/copykat/combat/results/'
output <- './Results/Tian/cross_protocols/copykat/combat/summary/'

load(paste0(input,'s_score_fit.rdata'))
cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
cls3_table <- table(cell_cls[,2],cls[,1])
bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(s_score_fit)

bar_col <- data.frame(id=1:nrow(cell_cls),pos=0.4,cell_cls,cell_batch=cell_cls[,2])
bar_col2 <- data.frame(id=1:nrow(cell_cls),pos=0.2,cell_cls,cell_batch=cell_cls[,3])
bar_col3 <- rbind(bar_col,bar_col2)

p7 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=10, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'CopyKAT after ComBat') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p7_col <- ggplot(bar_col3, aes(x=id, y=pos, fill=cell_batch)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=c(col[1:5],col2[1:4])) + theme(axis.title = element_blank(), axis.ticks = element_blank(), 
axis.text = element_blank(), legend.position='none', panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

gp7 <- ggplotGrob(p7)
gp7_col <- ggplotGrob(p7_col)  

maxWidth <- grid::unit.pmax(gp7$widths[2:5], gp7_col$widths[2:5])
gp7$widths[2:5] <- as.list(maxWidth)
gp7_col$widths[2:5] <- as.list(maxWidth)

p7_all_combat <- plot_grid(gp7, gp7_col, ncol=1,rel_heights=c(6,1))



p1_all <- plot_grid(p1_v2,p1_v3,p1_v4,p1_legend, labels= c('a','b','c'), nrow=1, rel_widths=c(3,3,3,2))
p8 <- plot_grid(p2_all, p3_all, p6_all, p7_all, labels=c('d','e','f','g'), nrow=1)
p9 <- plot_grid(p2_all_limma, p3_all_limma, p6_all_limma, p7_all_limma, labels=c('h','i','j','k'), nrow=1)
p10 <- plot_grid(p2_all_combat, p3_all_combat, p6_all_combat, p7_all_combat, labels=c('l','m','n','o') , nrow=1)
p11 <- plot_grid(p4_all, p5_all, labels=c('p','q'), ncol=4, nrow=1)


jpeg('./Manuscript/figure_table/figure5_v4.jpeg',width=12,height=14,res=300,units='in')
grid.arrange(p1_all, p8, p9, p10, p11, p_legend_all, nrow=6, heights=c(5,5,5,5,5,1))
dev.off()

pdf('./Manuscript/figure_table/figure5_v4.pdf',width=12,height=14)
grid.arrange(p1_all, p8, p9, p10, p11, p_legend_all, nrow=6, heights=c(5,5,5,5,5,1))
dev.off()









