
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

col_type <- brewer.pal(9,'Set1')

## barplot for true clusters ##

da_process <- function(da)
{

Score <- da[,3]
Score[!is.na(da[,1])] <- da[!is.na(da[,1]),1]
da <- cbind(da[,9:10],Score)
levels(da[,2])[levels(da[,2])=='10X_3cl_220'] <- '10x_3cl_CellRangerV2'
levels(da[,2])[levels(da[,2])=='10X_3cl_310'] <- '10x_3cl_CellRangerV3'
levels(da[,2])[levels(da[,2])=='10X_5cl_201'] <- '10x_5cl_CellRangerV2'
levels(da[,2])[levels(da[,2])=='10X_5cl_310'] <- '10x_5cl_CellRangerV3'
da[,2] <- factor(da[,2],levels=c('Drop-Seq','Cel-Seq2','10x_3cl_CellRangerV2','10x_3cl_CellRangerV3','10x_5cl_CellRangerV2','10x_5cl_CellRangerV3'))
da[,1] <- factor(da[,1],levels=c('inferCNV','CaSpER','sciCNV','CopyKAT','HoneyBADGER-expression','HoneyBADGER-allele'))
return(da)
}



ari <- read.csv('./Results/Tian/ARI.csv')
da_ari <- da_process(ari)
temp <- cbind(expand.grid(Methods=levels(da_ari$Methods), Datasets=levels(da_ari$Datasets)),Score=NA)
da_ari <- rbind(da_ari,temp[c(17,23,29,30,35,36),])
da_ari <- cbind(da_ari,'ARI')
colnames(da_ari)[4] <- 'Metrics'

fm <- read.csv('./Results/Tian/FM.csv')
da_fm <- da_process(fm)
temp <- cbind(expand.grid(Methods=levels(da_fm$Methods), Datasets=levels(da_fm$Datasets)),Score=NA)
da_fm <- rbind(da_fm,temp[c(17,23,29,30,35,36),])
da_fm <- cbind(da_fm,'FM')
colnames(da_fm)[4] <- 'Metrics'

nmi <- read.csv('./Results/Tian/NMI.csv')
da_nmi <- da_process(nmi)
temp <- cbind(expand.grid(Methods=levels(da_nmi$Methods), Datasets=levels(da_nmi$Datasets)),Score=NA)
da_nmi <- rbind(da_nmi,temp[c(17,23,29,30,35,36),])
da_nmi <- cbind(da_nmi,'NMI')
colnames(da_nmi)[4] <- 'Metrics'

vm <- read.csv('./Results/Tian/VM.csv')
da_vm <- da_process(vm)
temp <- cbind(expand.grid(Methods=levels(da_vm$Methods), Datasets=levels(da_vm$Datasets)),Score=NA)
da_vm <- rbind(da_vm,temp[c(17,23,29,30,35,36),])
da_vm <- cbind(da_vm,'VM')
colnames(da_vm)[4] <- 'Metrics'

da <- rbind(da_ari,da_fm,da_nmi,da_vm)

p1 <- ggplot(da, aes(x=Methods, y=Score, fill=Datasets)) + geom_bar(stat='identity',position=position_dodge()) +
scale_fill_manual(values=col_type) + facet_grid(Metrics ~.) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times'), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'),  
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))



col <- brewer.pal(9,'Set1')
col2 <- brewer.pal(8,'Set2')

##############
## inferCNV ##
##############

## cellrangerv2 ##

input <- './Results/Tian/infercnv/10x_3cl_220/results/'
output <- './Results/Tian/infercnv/10x_3cl_220/summary/'
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_220.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,2] <- paste0(cell_cls[,2],'_CellRangerV2')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))
cell_cls[,2] <- factor(cell_cls[,2],levels=c('H2228_CellRangerV2','H1975_CellRangerV2','HCC827_CellRangerV2'))

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
bar_col <- col[cell_cls[,1]]
bar_col2 <- data.frame(id=1:nrow(cell_cls),cell_cls,bar_col)

p2 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=3, value=col2[1:3]) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title='inferCNV: 10x_3cl_CellRangerV2') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p2_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=cell_line)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + labs(fill = 'Cell lines') + scale_fill_manual(values=col) + 
theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p_legend <- get_legend(p2_col)
p2_col <- p2_col + theme(legend.position='none')

gp2 <- ggplotGrob(p2)
gp2_col <- ggplotGrob(p2_col)  

maxWidth <- grid::unit.pmax(gp2$widths[2:5], gp2_col$widths[2:5])
gp2$widths[2:5] <- as.list(maxWidth)
gp2_col$widths[2:5] <- as.list(maxWidth)

p2_all <- plot_grid(gp2, gp2_col, ncol=1,rel_heights=c(6,1))




## cellrangerv3 ##

input <- './Results/Tian/infercnv/10x_3cl_310/results/'
output <- './Results/Tian/infercnv/10x_3cl_310/summary/'
cell_cls2 <- read.csv('./Results/Tian/cluster/sc_10x_310.csv',row.names=1)
rownames(cell_cls2) <- paste0(rownames(cell_cls2),'_1')
index <- match(rownames(cell_cls2),rownames(cell_cls))
cell_cls2[,2] <- paste0(cell_cls2[,2],'_CellRangerV2')
cell_cls2[is.na(index),2] <- sub('V2','V3',(cell_cls2[is.na(index),2]))

cell_cls2[,1] <- as.integer(factor(cell_cls2[,1],levels=c(2,0,1)))
cell_cls2[is.na(index),1] <- cell_cls2[is.na(index),1]+3
cell_cls2[,2] <- factor(cell_cls2[,2],levels=c('H2228_CellRangerV2','H1975_CellRangerV2','HCC827_CellRangerV2', 'H2228_CellRangerV3','H1975_CellRangerV3','HCC827_CellRangerV3'))
col <- brewer.pal(9,'Set1')

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls2 <- cell_cls2[rownames(cls),]
bar_col <- col[cell_cls2[,1]]
bar_col2 <- data.frame(id=1:nrow(cell_cls2),cell_cls2,bar_col)

p3 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=4, value=col2[1:4]) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title='inferCNV: 10x_3cl_CellRangerV3') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))


p3_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=cell_line)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + labs(fill = 'Cell lines') + scale_fill_manual(values=col) +
theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p_legend <- get_legend(p3_col)
p3_col <- p3_col + theme(legend.position='none')

gp3 <- ggplotGrob(p3)
gp3_col <- ggplotGrob(p3_col)  

maxWidth <- grid::unit.pmax(gp3$widths[2:5], gp3_col$widths[2:5])
gp3$widths[2:5] <- as.list(maxWidth)
gp3_col$widths[2:5] <- as.list(maxWidth)

p3_all <- plot_grid(gp3, gp3_col, ncol=1,rel_heights=c(6,1))



############
## scicnv ##
############

## cellrangerv2 ##

input <- './Results/Tian/scicnv/10x_3cl_220/summary/'
input2 <- './Results/Tian/infercnv/10x_3cl_220/'
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_220.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))

## cluster ##

load(paste0(input,'s_score_fit.rdata'))
dend <- as.dendrogram(s_score_fit)
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
bar_col <- col[cell_cls[,1]]
bar_col2 <- data.frame(id=1:nrow(cell_cls),cell_cls,bar_col)

p4 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=3, value=col2[1:3]) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title='sciCNV: 10x_3cl_CellRangerV2') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))


p4_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=cell_line)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + labs(fill = 'Cell lines') + scale_fill_manual(values=col) +
theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p_legend <- get_legend(p4_col)
p4_col <- p4_col + theme(legend.position='none')

gp4 <- ggplotGrob(p4)
gp4_col <- ggplotGrob(p4_col)  

maxWidth <- grid::unit.pmax(gp4$widths[2:5], gp4_col$widths[2:5])
gp4$widths[2:5] <- as.list(maxWidth)
gp4_col$widths[2:5] <- as.list(maxWidth)

p4_all <- plot_grid(gp4, gp4_col, ncol=1,rel_heights=c(6,1))




## cellrangerv3 ##

input <- './Results/Tian/scicnv/10x_3cl_310/summary/'
cell_cls2 <- read.csv('./Results/Tian/cluster/sc_10x_310.csv',row.names=1)
rownames(cell_cls2) <- paste0(rownames(cell_cls2),'_1')
index <- match(rownames(cell_cls2),rownames(cell_cls))
cell_cls2[,2] <- paste0(cell_cls2[,2],'_CellRangerV2')
cell_cls2[is.na(index),2] <- sub('V2','V3',(cell_cls2[is.na(index),2]))

cell_cls2[,1] <- as.integer(factor(cell_cls2[,1],levels=c(2,0,1)))
cell_cls2[is.na(index),1] <- cell_cls2[is.na(index),1]+3
cell_cls2[,2] <- factor(cell_cls2[,2],levels=c('H2228_CellRangerV2','H1975_CellRangerV2','HCC827_CellRangerV2', 'H2228_CellRangerV3','H1975_CellRangerV3','HCC827_CellRangerV3'))
col <- brewer.pal(9,'Set1')

load(paste0(input,'s_score_fit.rdata'))
dend <- as.dendrogram(s_score_fit)
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls2 <- cell_cls2[rownames(s_score_cls),]
bar_col <- col[cell_cls2[,1]]
bar_col2 <- data.frame(id=1:nrow(cell_cls2),cell_cls2,bar_col)

p5 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=4, value=col2[1:4]) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title='sciCNV: 10x_3cl_CellRangerV3') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))


p5_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=cell_line)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + labs(fill = 'Cell lines') + scale_fill_manual(values=col) +
theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p_legend <- get_legend(p5_col)
p5_col <- p5_col + theme(legend.position='none')

gp5 <- ggplotGrob(p5)
gp5_col <- ggplotGrob(p5_col)  

maxWidth <- grid::unit.pmax(gp5$widths[2:5], gp5_col$widths[2:5])
gp5$widths[2:5] <- as.list(maxWidth)
gp5_col$widths[2:5] <- as.list(maxWidth)

p5_all <- plot_grid(gp5, gp5_col, ncol=1,rel_heights=c(6,1))




#############
## copykat ##
#############


## cellrangerv2 ##

input <- './Results/Tian/copykat/10x_3cl_220/summary/'
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_220.csv',row.names=1)
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))

## cluster ##

load(paste0(input,'s_score_fit.rdata'))
dend <- as.dendrogram(s_score_fit)
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls <- cell_cls[rownames(s_score_cls),]
bar_col <- col[cell_cls[,1]]
bar_col2 <- data.frame(id=1:nrow(cell_cls),cell_cls,bar_col)

p6 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=3, value=col2[1:3]) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title='CopyKAT: 10x_3cl_CellRangerV2') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))


p6_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=cell_line)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + labs(fill = 'Cell lines') + scale_fill_manual(values=col) +
theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), 
panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p_legend <- get_legend(p6_col)
p6_col <- p6_col + theme(legend.position='none')

gp6 <- ggplotGrob(p6)
gp6_col <- ggplotGrob(p6_col)  

maxWidth <- grid::unit.pmax(gp6$widths[2:5], gp6_col$widths[2:5])
gp6$widths[2:5] <- as.list(maxWidth)
gp6_col$widths[2:5] <- as.list(maxWidth)

p6_all <- plot_grid(gp6, gp6_col, ncol=1,rel_heights=c(6,1))



## cellrangerv3 ##

input <- './Results/Tian/copykat/10x_3cl_310/summary/'
cell_cls2 <- read.csv('./Results/Tian/cluster/sc_10x_310.csv',row.names=1)
rownames(cell_cls2) <- paste0(rownames(cell_cls2),'_1')
index <- match(rownames(cell_cls2),rownames(cell_cls))
cell_cls2[,2] <- paste0(cell_cls2[,2],'_CellRangerV2')
cell_cls2[is.na(index),2] <- sub('V2','V3',(cell_cls2[is.na(index),2]))

cell_cls2[,1] <- as.integer(factor(cell_cls2[,1],levels=c(2,0,1)))
cell_cls2[is.na(index),1] <- cell_cls2[is.na(index),1]+3
cell_cls2[,2] <- factor(cell_cls2[,2],levels=c('H2228_CellRangerV2','H1975_CellRangerV2','HCC827_CellRangerV2', 'H2228_CellRangerV3','H1975_CellRangerV3','HCC827_CellRangerV3'))
col <- brewer.pal(9,'Set1')

load(paste0(input,'s_score_fit.rdata'))
dend <- as.dendrogram(s_score_fit)
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
cell_cls2 <- cell_cls2[rownames(s_score_cls),]
bar_col <- col[cell_cls2[,1]]
bar_col2 <- data.frame(id=1:nrow(cell_cls2),cell_cls2,bar_col)

p7 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=4, value=col2[1:4]) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title='CopyKAT: 10x_3cl_CellRangerV3') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))


p7_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=cell_line)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + labs(fill = 'Cell lines') + scale_fill_manual(values=col) +
theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position='bottom',
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p_legend <- get_legend(p7_col)
p7_col <- p7_col + theme(legend.position='none')

gp7 <- ggplotGrob(p7)
gp7_col <- ggplotGrob(p7_col)  

maxWidth <- grid::unit.pmax(gp7$widths[2:5], gp7_col$widths[2:5])
gp7$widths[2:5] <- as.list(maxWidth)
gp7_col$widths[2:5] <- as.list(maxWidth)

p7_all <- plot_grid(gp7, gp7_col, ncol=1,rel_heights=c(6,1))



## ari score associated with cluster 3-10 ##

ari <- read.csv('./Results/Tian/ARI.csv')
ari <- ari[-grep('Honey',ari[,9]),]
Cluster <- 3:10
ari2 <- c()
for(i in 1:nrow(ari))
{
	temp <- cbind.data.frame(Cluster,as.numeric(ari[i,1:8]),ari[i,9],ari[i,10])
	ari2 <- rbind(ari2,temp)
}

colnames(ari2)[2:4] <- c('ARI','Methods','Datasets')
levels(ari2[,4])[levels(ari2[,4])=='10X_3cl_220'] <- '10x_3cl_CellRangerV2'
levels(ari2[,4])[levels(ari2[,4])=='10X_3cl_310'] <- '10x_3cl_CellRangerV3'
levels(ari2[,4])[levels(ari2[,4])=='10X_5cl_201'] <- '10x_5cl_CellRangerV2'
levels(ari2[,4])[levels(ari2[,4])=='10X_5cl_310'] <- '10x_5cl_CellRangerV3'
ari2[,4] <- factor(ari2[,4],levels=c('Drop-Seq','Cel-Seq2','10x_3cl_CellRangerV2','10x_3cl_CellRangerV3','10x_5cl_CellRangerV2','10x_5cl_CellRangerV3'))

ari_infer <- ari2[ari2[,3]=='inferCNV',]
ari_casper <- ari2[ari2[,3]=='CaSpER',]
ari_sci <- ari2[ari2[,3]=='sciCNV',]
ari_copy <- ari2[ari2[,3]=='CopyKAT',]

p8 <- ggplot(ari_infer, aes(x=Cluster,y=ARI, col=Datasets))+ geom_line() + labs(title='inferCNV') + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(3,10), breaks = seq(3,10,by=1)) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times',hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), legend.position='bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))
p8_legend <- get_legend(p8)

p9 <- ggplot(ari_casper, aes(x=Cluster,y=ARI, col=Datasets))+ geom_line() + labs(title='CaSpER') + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(3,10), breaks = seq(3,10,by=1)) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times',hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p10 <- ggplot(ari_sci, aes(x=Cluster,y=ARI, col=Datasets))+ geom_line() + labs(title='sciCNV') + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(3,10), breaks = seq(3,10,by=1)) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times',hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p11 <- ggplot(ari_copy, aes(x=Cluster,y=ARI, col=Datasets))+ geom_line() + labs(title='CopyKAT') + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(3,10), breaks = seq(3,10,by=1)) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times',hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p8 <- p8 + theme(legend.position='none')
p9 <- p9 + theme(legend.position='none')
p10 <- p10 + theme(legend.position='none')
p11 <- p11 + theme(legend.position='none')



## cell contamination check ##

input <- './Results/Tian/decontx/'

load(paste0(input,'res_3cl_220_celltype.rdata'))
conp_3cl_220 <- cbind.data.frame(results[[2]][[3]],'10x_3cl_CellRangerV2')
load(paste0(input,'res_3cl_310_celltype.rdata'))
conp_3cl_310 <- cbind.data.frame(results[[2]][[3]],'10x_3cl_CellRangerV3')
index <- match(rownames(conp_3cl_310),rownames(conp_3cl_220))
conp_3cl_310_1 <- cbind.data.frame(conp_3cl_310[!is.na(index),1],'Cells_CellRangerV2')
conp_3cl_310_2 <- cbind.data.frame(conp_3cl_310[is.na(index),1],'Cells_CellRangerV3')
load(paste0(input,'res_5cl_201_celltype.rdata'))
conp_5cl_201 <- cbind.data.frame(results[[2]][[3]],'10x_5cl_CellRangerV2')
load(paste0(input,'res_5cl_310_celltype.rdata'))
conp_5cl_310 <- cbind.data.frame(results[[2]][[3]],'10x_5cl_CellRangerV3')

colnames(conp_3cl_220) <- colnames(conp_3cl_310) <- colnames(conp_3cl_310_1) <- colnames(conp_3cl_310_2) <- colnames(conp_5cl_201) <- colnames(conp_5cl_310) <- c('Est_conp','Datasets')
da <- rbind(conp_3cl_220,conp_3cl_310,conp_5cl_201,conp_5cl_310)


p16 <- ggplot(da, aes(x=Datasets,y=Est_conp, col=Datasets))+ geom_violin() + labs(title='Estimated contamination of 10x data', x='',y='Estimated contamination (%)') + scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
scale_color_manual(values=col_type) + theme(axis.text = element_text(size = 9,family='Times'), 
plot.title = element_text(size = 9,family='Times',hjust=0.5), axis.title.x = element_text(size = 9,family='Times'), 
axis.text.x = element_text(angle=10, hjust=0.5,vjust=0.5), axis.title.y = element_text(size = 9,family='Times'), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none',
panel.background = element_rect(fill='white',colour='black'))


p16 <- ggplot(da, aes(x=Datasets,y=Est_conp, col=Datasets))+ geom_violin() + labs(title='Estimated contamination of 10x data', x='',y='Estimated contamination (%)') + scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
scale_color_manual(values=col_type) + theme(axis.text = element_text(size = 9,family='Times'), 
plot.title = element_text(size = 9,family='Times',hjust=0.5), axis.title.x = element_text(size = 9,family='Times'), 
axis.title.y = element_text(size = 9,family='Times'), 
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position='none',
panel.background = element_rect(fill='white',colour='black'))









jpeg('./Manuscript/figure_table/figure4_v2.jpeg',width=10,height=6,res=300,units='in')
p1
dev.off()

pdf('./Manuscript/figure_table/figure4_v2.pdf',width=10,height=6)
p1
dev.off()



p12 <- plot_grid(p2_all,p3_all,p4_all,p5_all,p6_all,p7_all, labels=c('a','b','c','d','e','f'),ncol=2)
p13 <- plot_grid(p12,p_legend,nrow=3,p16,rel_heights=c(6,1,2),labels=c('','','g'))

jpeg('./Manuscript/figure_table/sup_figure/sup_fig6_v3.jpeg',width=8,height=10,res=300,units='in')
p13
dev.off()

pdf('./Manuscript/figure_table/sup_figure/sup_fig6_v3.pdf',width=8,height=10)
p13
dev.off()


p14 <- plot_grid(p8,p9,p10,p11, labels=c('a','b','c','d'),ncol=2)
p15 <- plot_grid(p14,p8_legend,nrow=2,rel_heights=c(6,1))


