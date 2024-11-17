
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


input <- './Results/sens_spec_eva/cyto/sclc/'
method <- c('inferCNV','CaSpER','sciCNV','CopyKAT')
method2 <- tolower(method)

da_all <- c()
da_gain <- c()
da_loss <- c()

for(i in 1:length(method))
{
	temp <- read.csv(paste0(input,'sens_spec_all_',method2[i],'.csv'))
	temp <- cbind.data.frame(temp,method[i])
	da_all <- rbind(da_all,temp[c(1,9),])
	temp <- read.csv(paste0(input,'sens_spec_gain_',method2[i],'.csv'))
	temp <- cbind.data.frame(temp,method[i])
	da_gain <- rbind(da_gain,temp[c(1,9),])
	temp <- read.csv(paste0(input,'sens_spec_loss_',method2[i],'.csv'))
	temp <- cbind.data.frame(temp,method[i])
	da_loss <- rbind(da_loss,temp[c(1,9),])
}
colnames(da_all)[c(1,2,5)] <- colnames(da_loss)[c(1,2,5)] <- colnames(da_gain)[c(1,2,5)] <- c('Sensitivity','Specificity','Methods')


p1 <- ggplot(da_all,aes(x = Specificity, y = Sensitivity, color=Status, shape=Methods)) + geom_point(size=4) +
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) +
scale_colour_manual(values=col_type) + labs(title='Cytobands for all recurrent CNVs') +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))
p1_legend <- get_legend(p1)


p2 <- ggplot(da_gain,aes(x = Specificity, y = Sensitivity, color=Status, shape=Methods)) + geom_point(size=4) +
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) +
scale_colour_manual(values=col_type) + labs(title='Cytobands for CNV gains') +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


p3 <- ggplot(da_loss,aes(x = Specificity, y = Sensitivity, color=Status, shape=Methods)) + geom_point(size=4) +
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) +
scale_colour_manual(values=col_type) + labs(title='Cytobands for CNV losses') +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))



################
## clustering ##
################

col <- brewer.pal(9,'Set1')
col2 <- brewer.pal(8,'Set2')

## infercnv ##

input <- './Results/SCLC/infercnv/results2/'

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
cell_cls <- data.frame(matrix(c(rep('Primary',93),rep('Relapse',39)),ncol=1))
colnames(cell_cls) <- 'Status'
rownames(cell_cls) <- rownames(cls)
cell_cls_org <- cell_cls
bar_col <- col[as.integer(cell_cls[,1])]
bar_col2 <- data.frame(id=1:nrow(cell_cls),cell_cls,bar_col)


p4 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=5, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title='inferCNV') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p4_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=Status)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=col) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), 
legend.justification = c("left", "center"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))

p4_legend <- get_legend(p4_col)
p4_col <- p4_col + theme(legend.position='none')

gp4 <- ggplotGrob(p4)
gp4_col <- ggplotGrob(p4_col)  

maxWidth <- grid::unit.pmax(gp4$widths[2:5], gp4_col$widths[2:5])
gp4$widths[2:5] <- as.list(maxWidth)
gp4_col$widths[2:5] <- as.list(maxWidth)

p4_all <- plot_grid(gp4, gp4_col, ncol=1,rel_heights=c(6,1))



## casper ##

input <- './Results/SCLC/casper/summary/'
load(paste0(input,'s_exprs_fit.rdata'))
dend <- as.dendrogram(s_exprs_fit)
s_exprs_cls <- cutree(s_exprs_fit,3:12,order_clusters_as_data=F)
colnames(s_exprs_cls) <- paste0('pred_cls',colnames(s_exprs_cls))
s_exprs_cls2 <- read.csv(paste0(input,'cls_all.csv'),row.names=1)
cell_cls <- s_exprs_cls2[rownames(s_exprs_cls),]
bar_col <- col[as.integer(cell_cls[,1])]
bar_col2 <- data.frame(id=1:nrow(cell_cls),cell_cls,bar_col)

p5 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=5, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'CaSpER') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p5_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=Status)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=col) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))
p5_col <- p5_col + theme(legend.position='none')

gp5 <- ggplotGrob(p5)
gp5_col <- ggplotGrob(p5_col)  

maxWidth <- grid::unit.pmax(gp5$widths[2:5], gp5_col$widths[2:5])
gp5$widths[2:5] <- as.list(maxWidth)
gp5_col$widths[2:5] <- as.list(maxWidth)

p5_all <- plot_grid(gp5, gp5_col, ncol=1,rel_heights=c(6,1))



## sciCNV ##

input <- './Results/SCLC/scicnv/summary/'
load(paste0(input,'s_score_fit.rdata'))
dend <- as.dendrogram(s_score_fit)
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
s_score_cls2 <- read.csv(paste0(input,'cls_all.csv'),row.names=1)
cell_cls <- s_score_cls2[rownames(s_score_cls),]
bar_col <- col[as.integer(cell_cls[,1])]
bar_col2 <- data.frame(id=1:nrow(cell_cls),cell_cls,bar_col)

p6 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=5, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'sciCNV') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p6_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=Status)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=col) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))
p6_col <- p6_col + theme(legend.position='none')

gp6 <- ggplotGrob(p6)
gp6_col <- ggplotGrob(p6_col)  

maxWidth <- grid::unit.pmax(gp6$widths[2:5], gp6_col$widths[2:5])
gp6$widths[2:5] <- as.list(maxWidth)
gp6_col$widths[2:5] <- as.list(maxWidth)

p6_all <- plot_grid(gp6, gp6_col, ncol=1,rel_heights=c(6,1))




## copykat ##

input <- './Results/SCLC/copykat/summary/'
load(paste0(input,'s_score_fit.rdata'))
dend <- as.dendrogram(s_score_fit)
s_score_cls <- cutree(s_score_fit,3:12,order_clusters_as_data=F)
colnames(s_score_cls) <- paste0('pred_cls',colnames(s_score_cls))
s_score_cls2 <- read.csv(paste0(input,'cls_all.csv'),row.names=1)
cell_cls <- s_score_cls2[rownames(s_score_cls),]
bar_col <- col[as.integer(cell_cls[,1])]
bar_col2 <- data.frame(id=1:nrow(cell_cls),cell_cls,bar_col)

p7 <- dend %>% as.dendrogram() %>% set('labels','') %>% set('branches_k_color',k=5, value=col2) %>% set('branches_lwd', 0.5) %>% as.ggdend() %>% ggplot() + labs(title = 'CopyKAT') + theme(plot.title = element_text(size=9,family='Times',hjust=0.5))

p7_col <- ggplot(bar_col2, aes(x=id, y=0.2, fill=Status)) + geom_tile() + scale_y_continuous(expand=c(0,0)) + 
scale_fill_manual(values=col) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='white'))
p7_col <- p7_col + theme(legend.position='none')

gp7 <- ggplotGrob(p7)
gp7_col <- ggplotGrob(p7_col)  

maxWidth <- grid::unit.pmax(gp7$widths[2:5], gp7_col$widths[2:5])
gp7$widths[2:5] <- as.list(maxWidth)
gp7_col$widths[2:5] <- as.list(maxWidth)

p7_all <- plot_grid(gp7, gp7_col, ncol=1,rel_heights=c(6,1))


p1 <- p1 + theme(legend.position='none')
p2 <- p2 + theme(legend.position='none')
p3 <- p3 + theme(legend.position='none')
p8 <- plot_grid(p1,p2,p3,p1_legend,nrow=1,labels=c('a','b','c'),rel_widths=c(4,4,4,1.5))

p9 <- plot_grid(p4_all,p5_all,p6_all,p7_all,nrow=2,labels=c('d','e','f','g'))
p10 <- plot_grid(p9,p4_legend,nrow=1,rel_widths=c(12,1.5))


jpeg('./Manuscript/figure_table/figure6_v3.jpeg',width=8,height=8,res=300,units='in')
grid.arrange(p8, p10, nrow=2, heights=c(1,2))
dev.off()

pdf('./Manuscript/figure_table/figure6_v3.pdf',width=8,height=8)
grid.arrange(p8, p10, nrow=2, heights=c(1,2))
dev.off()

p11 <- plot_grid(p1,p1_legend,nrow=1,labels=c('a'),rel_widths=c(4,1))
p12 <- plot_grid(p4_all,p5_all,p6_all,p7_all,nrow=2,labels=c('b','c','d','e'))
p13 <- plot_grid(p12,p4_legend,nrow=1,rel_widths=c(4,1))


jpeg('./Manuscript/figure_table/figure6_v4.jpeg',width=8,height=8,res=300,units='in')
grid.arrange(p11, p13, nrow=2, heights=c(1.5,2))
dev.off()

pdf('./Manuscript/figure_table/figure6_v4.pdf',width=8,height=8)
grid.arrange(p11, p13, nrow=2, heights=c(1.5,2))
dev.off()



