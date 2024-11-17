

##############################################
## protocol and method of all cells ##
##############################################

library(ggpubr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(fossil)

col_type <- brewer.pal(9,'Set1')

input <- './Results/sens_spec_eva/cyto/all_cell/ref_v1/v1/'
input2 <- './Results/sens_spec_eva/cyto/all_cell/ref_v2/v1/'


## all cnvs ##

infercnv_summary <- read.csv(paste0(input,'sens_spec_all_infercnv.csv'))
infercnv_summary <- cbind.data.frame(infercnv_summary[infercnv_summary[,3]=='sc_0.1',-3],'inferCNV')
casper_summary <- read.csv(paste0(input,'sens_spec_all_casper.csv'))
casper_summary <- cbind.data.frame(casper_summary[casper_summary[,3]=='sc_0.1',-3],'CaSpER')
scicnv_summary <- read.csv(paste0(input,'sens_spec_all_scicnv.csv'))
scicnv_summary <- cbind.data.frame(scicnv_summary[scicnv_summary[,3]=='sc_0.1',-3],'sciCNV')
copykat_summary <- read.csv(paste0(input,'sens_spec_all_copykat.csv'))
copykat_summary <- cbind.data.frame(copykat_summary[copykat_summary[,3]=='sc_0.1',-3],'CopyKAT')
hb_exprs_summary <- read.csv(paste0(input,'sens_spec_all_hb_exprs.csv'))
hb_exprs_summary <- cbind.data.frame(hb_exprs_summary[hb_exprs_summary[,3]=='sc_0.1',-3],'HoneyBADGER-expression')
hb_allele_summary <- read.csv(paste0(input,'sens_spec_all_hb_allele.csv'))
hb_allele_summary <- cbind.data.frame(hb_allele_summary[hb_allele_summary[,3]=='sc_0.1',-3],'HoneyBADGER-allele')

casper_summary2 <- read.csv(paste0(input2,'sens_spec_all_casper.csv'))
casper_summary2 <- cbind.data.frame(casper_summary2[casper_summary2[,3]=='sc_0.1',-3],'CaSpER')
copykat_summary2 <- read.csv(paste0(input2,'sens_spec_all_copykat.csv'))
copykat_summary2 <- cbind.data.frame(copykat_summary2[copykat_summary2[,3]=='sc_0.1',-3],'CopyKAT')
casper_summary <- rbind(casper_summary,casper_summary2[5,])
copykat_summary <- rbind(copykat_summary,copykat_summary2[5,])

colnames(infercnv_summary)[-3] <- colnames(casper_summary)[-3] <- colnames(scicnv_summary)[-3] <- colnames(copykat_summary)[-3] <- colnames(hb_exprs_summary)[-3] <- colnames(hb_allele_summary)[-3] <- c('Sensitivity','Specificity','Methods')
da <- rbind(infercnv_summary,casper_summary,scicnv_summary,copykat_summary,hb_exprs_summary,hb_allele_summary)
levels(da[,3])[5] <- 'Bulk RNA-seq'
da[,3] <- factor(da[,3],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq'))

p1 <- ggplot(da,aes(x = Specificity, y = Sensitivity, color=Protocols, shape=Methods)) + geom_point(size=4) +
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) +
scale_colour_manual(values=col_type) + labs(title='Cytobands for all recurrent CNVs') +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


## gain cnvs ##

infercnv_summary <- read.csv(paste0(input,'sens_spec_gain_infercnv.csv'))
infercnv_summary <- cbind.data.frame(infercnv_summary[infercnv_summary[,3]=='sc_0.1',-3],'inferCNV')
casper_summary <- read.csv(paste0(input,'sens_spec_gain_casper.csv'))
casper_summary <- cbind.data.frame(casper_summary[casper_summary[,3]=='sc_0.1',-3],'CaSpER')
scicnv_summary <- read.csv(paste0(input,'sens_spec_gain_scicnv.csv'))
scicnv_summary <- cbind.data.frame(scicnv_summary[scicnv_summary[,3]=='sc_0.1',-3],'sciCNV')
copykat_summary <- read.csv(paste0(input,'sens_spec_gain_copykat.csv'))
copykat_summary <- cbind.data.frame(copykat_summary[copykat_summary[,3]=='sc_0.1',-3],'CopyKAT')
hb_exprs_summary <- read.csv(paste0(input,'sens_spec_gain_hb_exprs.csv'))
hb_exprs_summary <- cbind.data.frame(hb_exprs_summary[hb_exprs_summary[,3]=='sc_0.1',-3],'HoneyBADGER-expression')
hb_allele_summary <- read.csv(paste0(input,'sens_spec_gain_hb_allele.csv'))
hb_allele_summary <- cbind.data.frame(hb_allele_summary[hb_allele_summary[,3]=='sc_0.1',-3],'HoneyBADGER-allele')

casper_summary2 <- read.csv(paste0(input2,'sens_spec_gain_casper.csv'))
casper_summary2 <- cbind.data.frame(casper_summary2[casper_summary2[,3]=='sc_0.1',-3],'CaSpER')
copykat_summary2 <- read.csv(paste0(input2,'sens_spec_gain_copykat.csv'))
copykat_summary2 <- cbind.data.frame(copykat_summary2[copykat_summary2[,3]=='sc_0.1',-3],'CopyKAT')
casper_summary <- rbind(casper_summary,casper_summary2[5,])
copykat_summary <- rbind(copykat_summary,copykat_summary2[5,])

colnames(infercnv_summary)[-3] <- colnames(casper_summary)[-3] <- colnames(scicnv_summary)[-3] <- colnames(copykat_summary)[-3] <- colnames(hb_exprs_summary)[-3] <- colnames(hb_allele_summary)[-3] <- c('Sensitivity','Specificity','Methods')
da <- rbind(infercnv_summary,casper_summary,scicnv_summary,copykat_summary,hb_exprs_summary,hb_allele_summary)
levels(da[,3])[5] <- 'Bulk RNA-seq'
da[,3] <- factor(da[,3],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq'))

p2 <- ggplot(da,aes(x = Specificity, y = Sensitivity, color=Protocols, shape=Methods)) + geom_point(size=4) +
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) +
scale_colour_manual(values=col_type) + labs(title='Cytobands for CNV gains') + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))



## loss cnvs ##

infercnv_summary <- read.csv(paste0(input,'sens_spec_loss_infercnv.csv'))
infercnv_summary <- cbind.data.frame(infercnv_summary[infercnv_summary[,3]=='sc_0.1',-3],'inferCNV')
casper_summary <- read.csv(paste0(input,'sens_spec_loss_casper.csv'))
casper_summary <- cbind.data.frame(casper_summary[casper_summary[,3]=='sc_0.1',-3],'CaSpER')
scicnv_summary <- read.csv(paste0(input,'sens_spec_loss_scicnv.csv'))
scicnv_summary <- cbind.data.frame(scicnv_summary[scicnv_summary[,3]=='sc_0.1',-3],'sciCNV')
copykat_summary <- read.csv(paste0(input,'sens_spec_loss_copykat.csv'))
copykat_summary <- cbind.data.frame(copykat_summary[copykat_summary[,3]=='sc_0.1',-3],'CopyKAT')
hb_exprs_summary <- read.csv(paste0(input,'sens_spec_loss_hb_exprs.csv'))
hb_exprs_summary <- cbind.data.frame(hb_exprs_summary[hb_exprs_summary[,3]=='sc_0.1',-3],'HoneyBADGER-expression')
hb_allele_summary <- read.csv(paste0(input,'sens_spec_loss_hb_allele.csv'))
hb_allele_summary <- cbind.data.frame(hb_allele_summary[hb_allele_summary[,3]=='sc_0.1',-3],'HoneyBADGER-allele')

casper_summary2 <- read.csv(paste0(input2,'sens_spec_gain_casper.csv'))
casper_summary2 <- cbind.data.frame(casper_summary2[casper_summary2[,3]=='sc_0.1',-3],'CaSpER')
copykat_summary2 <- read.csv(paste0(input2,'sens_spec_gain_copykat.csv'))
copykat_summary2 <- cbind.data.frame(copykat_summary2[copykat_summary2[,3]=='sc_0.1',-3],'CopyKAT')
casper_summary <- rbind(casper_summary,casper_summary2[5,])
copykat_summary <- rbind(copykat_summary,copykat_summary2[5,])

colnames(infercnv_summary)[-3] <- colnames(casper_summary)[-3] <- colnames(scicnv_summary)[-3] <- colnames(copykat_summary)[-3] <- colnames(hb_exprs_summary)[-3] <- colnames(hb_allele_summary)[-3] <- c('Sensitivity','Specificity','Methods')
da <- rbind(infercnv_summary,casper_summary,scicnv_summary,copykat_summary,hb_exprs_summary,hb_allele_summary)
levels(da[,3])[5] <- 'Bulk RNA-seq'
da[,3] <- factor(da[,3],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq'))

p3 <- ggplot(da,aes(x = Specificity, y = Sensitivity, color=Protocols, shape=Methods)) + geom_point(size=4) +
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) +
scale_colour_manual(values=col_type) + labs(title='Cytobands for CNV losses') + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


##########################
## protocol consistency ##
##########################

## casper ##

## all cnvs ##

input <- './Results/10x/all_cell/casper/v1/sens_spec_summary/'
input2 <- './Results/c1_fda/all_cell/casper/v1/sens_spec_summary/'
input3 <- './Results/c1_llu/all_cell/casper/v1/sens_spec_summary/'
input4 <- './Results/WaferGen/all_cell/casper/v1/sens_spec_summary/'
input5 <- './Results/bulk/all_cell/casper/v2/sens_spec_summary/'

cnv_anno <- read.csv('./Anno/cyto/cnv_cyto_brca.csv')
cnv_10x <- read.csv(paste0(input,'cnv_allv1.csv'))
cnv_10x[,1] <- as.character(cnv_10x[,1])
cnv_10x[63,1] <- paste0(cnv_10x[63,1],'.2')
cnv_c1_ht <- read.csv(paste0(input2,'cnv_allv1.csv'))
cnv_c1_ht[,1] <- as.character(cnv_c1_ht[,1])
cnv_c1_ht[63,1] <- paste0(cnv_c1_ht[63,1],'.2')
cnv_c1 <- read.csv(paste0(input3,'cnv_allv1.csv'))
cnv_c1[,1] <- as.character(cnv_c1[,1])
cnv_c1[63,1] <- paste0(cnv_c1[63,1],'.2')
cnv_icell <- read.csv(paste0(input4,'cnv_allv1.csv'))
cnv_icell[,1] <- as.character(cnv_icell[,1])
cnv_icell[63,1] <- paste0(cnv_icell[63,1],'.2')
cnv_bulk <- read.csv(paste0(input5,'cnv_allv1.csv'))
cnv_bulk[,1] <- as.character(cnv_bulk[,1])
cnv_bulk[63,1] <- paste0(cnv_bulk[63,1],'.2')

cnv_10x <- cnv_10x[match(cnv_anno[,1],cnv_10x[,1]),]
cnv_c1_ht <- cnv_c1_ht[match(cnv_anno[,1],cnv_c1_ht[,1]),]
cnv_c1 <- cnv_c1[match(cnv_anno[,1],cnv_c1[,1]),]
cnv_icell <- cnv_icell[match(cnv_anno[,1],cnv_icell[,1]),]
cnv_bulk <- cnv_bulk[match(cnv_anno[,1],cnv_bulk[,1]),]

score_wgs <- as.integer(cnv_10x[,2])%%2*cnv_anno[,2]
score_10x <- as.integer(cnv_10x[,5])%%2*cnv_anno[,2]
score_c1_ht <- as.integer(cnv_c1_ht[,5])%%2*cnv_anno[,2]
score_c1 <- as.integer(cnv_c1[,5])%%2*cnv_anno[,2]
score_icell <- as.integer(cnv_icell[,5])%%2*cnv_anno[,2]
score_bulk <- as.integer(cnv_bulk[,5])%%2*cnv_anno[,2]

ari_wgs_10x <- adj.rand.index(score_wgs+1,score_10x+1)
ari_wgs_c1_ht <- adj.rand.index(score_wgs+1,score_c1_ht+1)
ari_wgs_c1 <- adj.rand.index(score_wgs+1,score_c1+1)
ari_wgs_icell <- adj.rand.index(score_wgs+1,score_icell+1)
ari_wgs_bulk <- adj.rand.index(score_wgs+1,score_bulk+1)
ari_10x_c1_ht <- adj.rand.index(score_10x+1,score_c1_ht+1)
ari_10x_c1 <- adj.rand.index(score_10x+1,score_c1+1)
ari_10x_icell <- adj.rand.index(score_10x+1,score_icell+1)
ari_10x_bulk <- adj.rand.index(score_10x+1,score_bulk+1)
ari_c1_ht_c1 <- adj.rand.index(score_c1_ht+1,score_c1+1)
ari_c1_ht_icell <- adj.rand.index(score_c1_ht+1,score_icell+1)
ari_c1_ht_bulk <- adj.rand.index(score_c1_ht+1,score_bulk+1)
ari_c1_icell <- adj.rand.index(score_c1+1,score_icell+1)
ari_c1_bulk <- adj.rand.index(score_c1+1,score_bulk+1)
ari_icell_bulk <- adj.rand.index(score_icell+1,score_bulk+1)

ARI <- c(ari_wgs_10x,ari_wgs_c1_ht,ari_wgs_c1,ari_wgs_icell,ari_wgs_bulk,ari_10x_c1_ht,ari_10x_c1,ari_10x_icell,ari_10x_bulk, ari_c1_ht_c1,ari_c1_ht_icell,ari_c1_ht_bulk,ari_c1_icell,ari_c1_bulk,ari_icell_bulk)
ARI <- rep(ARI,2)
Protocols <- c(rep('WGS',5),rep('10x',4),rep('C1-HT',3),rep('C1',2),'ICELL8','10x','C1-HT','C1','ICELL8','Bulk RNA-seq',
'C1-HT','C1','ICELL8','Bulk RNA-seq','C1','ICELL8','Bulk RNA-seq','ICELL8',rep('Bulk RNA-seq',2))
da_ari <- cbind.data.frame(ARI,Protocols)
da_ari[,2] <- factor(da_ari[,2],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq','WGS'))

da_wgs <- cbind.data.frame(cnv_anno[,1],score_wgs,'WGS')
da_10x <- cbind.data.frame(cnv_anno[,1],score_10x,'10x')
da_c1_ht <- cbind.data.frame(cnv_anno[,1],score_c1_ht,'C1-HT')
da_c1 <- cbind.data.frame(cnv_anno[,1],score_c1,'C1')
da_icell <- cbind.data.frame(cnv_anno[,1],score_icell,'ICELL8')
da_bulk <- cbind.data.frame(cnv_anno[,1],score_bulk,'Bulk RNA-seq')

colnames(da_wgs) <- colnames(da_10x) <- colnames(da_c1_ht) <- colnames(da_c1) <- colnames(da_icell) <- colnames(da_bulk) <- c('Cytoband','Status','Protocols')
da <- rbind(da_wgs,da_10x,da_c1_ht,da_c1,da_icell,da_bulk)
da[,1] <- factor(da[,1],levels=as.character(cnv_anno[,1]))
levels(da[,1])[27] <- levels(da[,1])[26]
da[,3] <- factor(da[,3],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq','WGS'))

p4 <- ggplot(da,aes(x=Cytoband,y=Status,fill=Protocols)) + labs(x = 'Cytoband',y='') + geom_bar(stat="identity") + 
scale_fill_manual(values=col_type) + theme(axis.text = element_text(size = 9,family='Times'), 
plot.title = element_text(size = 9,family='Times', hjust=0.5), axis.text.x = element_text(angle=90, hjust=1),
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "center"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p5 <-  ggplot(da_ari,aes(x=Protocols,y=ARI,fill=Protocols)) + labs(x = '',y='') + geom_boxplot() + 
scale_fill_manual(values=col_type) + scale_y_continuous(limits=c(0.6,0.9), breaks = seq(0.6,0.9,by=0.1)) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5),
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.position = 'bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


## gain cnvs ##

input <- './Results/10x/all_cell/casper/v1/sens_spec_summary/'
input2 <- './Results/c1_fda/all_cell/casper/v1/sens_spec_summary/'
input3 <- './Results/c1_llu/all_cell/casper/v1/sens_spec_summary/'
input4 <- './Results/WaferGen/all_cell/casper/v1/sens_spec_summary/'
input5 <- './Results/bulk/all_cell/casper/v2/sens_spec_summary/'

cnv_anno <- read.csv('./Anno/cyto/cnv_cyto_brca.csv')
cnv_10x <- read.csv(paste0(input,'cnv_gainv1.csv'))
cnv_10x[,1] <- as.character(cnv_10x[,1])
cnv_c1_ht <- read.csv(paste0(input2,'cnv_gainv1.csv'))
cnv_c1_ht[,1] <- as.character(cnv_c1_ht[,1])
cnv_c1 <- read.csv(paste0(input3,'cnv_gainv1.csv'))
cnv_c1[,1] <- as.character(cnv_c1[,1])
cnv_icell <- read.csv(paste0(input4,'cnv_gainv1.csv'))
cnv_icell[,1] <- as.character(cnv_icell[,1])
cnv_bulk <- read.csv(paste0(input5,'cnv_gainv1.csv'))
cnv_bulk[,1] <- as.character(cnv_bulk[,1])

index <- match(cnv_anno[,1],cnv_10x[,1])
cnv_anno <- cnv_anno[!is.na(index),]
cnv_10x <- cnv_10x[match(cnv_anno[,1],cnv_10x[,1]),]
cnv_c1_ht <- cnv_c1_ht[match(cnv_anno[,1],cnv_c1_ht[,1]),]
cnv_c1 <- cnv_c1[match(cnv_anno[,1],cnv_c1[,1]),]
cnv_icell <- cnv_icell[match(cnv_anno[,1],cnv_icell[,1]),]
cnv_bulk <- cnv_bulk[match(cnv_anno[,1],cnv_bulk[,1]),]

score_wgs <- as.integer(cnv_10x[,2])
score_10x <- as.integer(cnv_10x[,5])
score_c1_ht <- as.integer(cnv_c1_ht[,5])
score_c1 <- as.integer(cnv_c1[,5])
score_icell <- as.integer(cnv_icell[,5])
score_bulk <- as.integer(cnv_bulk[,5])

ari_wgs_10x <- adj.rand.index(score_wgs,score_10x)
ari_wgs_c1_ht <- adj.rand.index(score_wgs,score_c1_ht)
ari_wgs_c1 <- adj.rand.index(score_wgs,score_c1)
ari_wgs_icell <- adj.rand.index(score_wgs,score_icell)
ari_wgs_bulk <- adj.rand.index(score_wgs,score_bulk)
ari_10x_c1_ht <- adj.rand.index(score_10x,score_c1_ht)
ari_10x_c1 <- adj.rand.index(score_10x,score_c1)
ari_10x_icell <- adj.rand.index(score_10x,score_icell)
ari_10x_bulk <- adj.rand.index(score_10x,score_bulk)
ari_c1_ht_c1 <- adj.rand.index(score_c1_ht,score_c1)
ari_c1_ht_icell <- adj.rand.index(score_c1_ht,score_icell)
ari_c1_ht_bulk <- adj.rand.index(score_c1_ht,score_bulk)
ari_c1_icell <- adj.rand.index(score_c1,score_icell)
ari_c1_bulk <- adj.rand.index(score_c1,score_bulk)
ari_icell_bulk <- adj.rand.index(score_icell,score_bulk)

ARI <- c(ari_wgs_10x,ari_wgs_c1_ht,ari_wgs_c1,ari_wgs_icell,ari_wgs_bulk,ari_10x_c1_ht,ari_10x_c1,ari_10x_icell,ari_10x_bulk, ari_c1_ht_c1,ari_c1_ht_icell,ari_c1_ht_bulk,ari_c1_icell,ari_c1_bulk,ari_icell_bulk)
ARI <- rep(ARI,2)
Protocols <- c(rep('WGS',5),rep('10x',4),rep('C1-HT',3),rep('C1',2),'ICELL8','10x','C1-HT','C1','ICELL8','Bulk RNA-seq',
'C1-HT','C1','ICELL8','Bulk RNA-seq','C1','ICELL8','Bulk RNA-seq','ICELL8',rep('Bulk RNA-seq',2))
da_ari_gain <- cbind.data.frame(ARI,Protocols)
da_ari_gain[,2] <- factor(da_ari_gain[,2],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq','WGS'))

p5_gain <-  ggplot(da_ari_gain,aes(x=Protocols,y=ARI,fill=Protocols)) + labs(x = '',y='') + geom_boxplot() + 
scale_fill_manual(values=col_type) + scale_y_continuous(limits=c(0,0.7), breaks = seq(0,0.7,by=0.1)) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5),
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.position = 'bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))



## loss cnvs ##

input <- './Results/10x/all_cell/casper/v1/sens_spec_summary/'
input2 <- './Results/c1_fda/all_cell/casper/v1/sens_spec_summary/'
input3 <- './Results/c1_llu/all_cell/casper/v1/sens_spec_summary/'
input4 <- './Results/WaferGen/all_cell/casper/v1/sens_spec_summary/'
input5 <- './Results/bulk/all_cell/casper/v2/sens_spec_summary/'

cnv_anno <- read.csv('./Anno/cyto/cnv_cyto_brca.csv')
cnv_10x <- read.csv(paste0(input,'cnv_lossv1.csv'))
cnv_10x[,1] <- as.character(cnv_10x[,1])
cnv_10x[30,1] <- paste0(cnv_10x[30,1],'.2')
cnv_c1_ht <- read.csv(paste0(input2,'cnv_lossv1.csv'))
cnv_c1_ht[,1] <- as.character(cnv_c1_ht[,1])
cnv_c1_ht[30,1] <- paste0(cnv_c1_ht[30,1],'.2')
cnv_c1 <- read.csv(paste0(input3,'cnv_lossv1.csv'))
cnv_c1[,1] <- as.character(cnv_c1[,1])
cnv_c1[30,1] <- paste0(cnv_c1[30,1],'.2')
cnv_icell <- read.csv(paste0(input4,'cnv_lossv1.csv'))
cnv_icell[,1] <- as.character(cnv_icell[,1])
cnv_icell[30,1] <- paste0(cnv_icell[30,1],'.2')
cnv_bulk <- read.csv(paste0(input5,'cnv_lossv1.csv'))
cnv_bulk[,1] <- as.character(cnv_bulk[,1])
cnv_bulk[30,1] <- paste0(cnv_bulk[30,1],'.2')

index <- match(cnv_anno[,1],cnv_10x[,1])
cnv_anno <- cnv_anno[!is.na(index),]
cnv_10x <- cnv_10x[match(cnv_anno[,1],cnv_10x[,1]),]
cnv_c1_ht <- cnv_c1_ht[match(cnv_anno[,1],cnv_c1_ht[,1]),]
cnv_c1 <- cnv_c1[match(cnv_anno[,1],cnv_c1[,1]),]
cnv_icell <- cnv_icell[match(cnv_anno[,1],cnv_icell[,1]),]
cnv_bulk <- cnv_bulk[match(cnv_anno[,1],cnv_bulk[,1]),]

score_wgs <- as.integer(cnv_10x[,2])
score_10x <- as.integer(cnv_10x[,5])
score_c1_ht <- as.integer(cnv_c1_ht[,5])
score_c1 <- as.integer(cnv_c1[,5])
score_icell <- as.integer(cnv_icell[,5])
score_bulk <- as.integer(cnv_bulk[,5])

ari_wgs_10x <- adj.rand.index(score_wgs,score_10x)
ari_wgs_c1_ht <- adj.rand.index(score_wgs,score_c1_ht)
ari_wgs_c1 <- adj.rand.index(score_wgs,score_c1)
ari_wgs_icell <- adj.rand.index(score_wgs,score_icell)
ari_wgs_bulk <- adj.rand.index(score_wgs,score_bulk)
ari_10x_c1_ht <- adj.rand.index(score_10x,score_c1_ht)
ari_10x_c1 <- adj.rand.index(score_10x,score_c1)
ari_10x_icell <- adj.rand.index(score_10x,score_icell)
ari_10x_bulk <- adj.rand.index(score_10x,score_bulk)
ari_c1_ht_c1 <- adj.rand.index(score_c1_ht,score_c1)
ari_c1_ht_icell <- adj.rand.index(score_c1_ht,score_icell)
ari_c1_ht_bulk <- adj.rand.index(score_c1_ht,score_bulk)
ari_c1_icell <- adj.rand.index(score_c1,score_icell)
ari_c1_bulk <- adj.rand.index(score_c1,score_bulk)
ari_icell_bulk <- adj.rand.index(score_icell,score_bulk)

ARI <- c(ari_wgs_10x,ari_wgs_c1_ht,ari_wgs_c1,ari_wgs_icell,ari_wgs_bulk,ari_10x_c1_ht,ari_10x_c1,ari_10x_icell,ari_10x_bulk, ari_c1_ht_c1,ari_c1_ht_icell,ari_c1_ht_bulk,ari_c1_icell,ari_c1_bulk,ari_icell_bulk)
ARI <- rep(ARI,2)
Protocols <- c(rep('WGS',5),rep('10x',4),rep('C1-HT',3),rep('C1',2),'ICELL8','10x','C1-HT','C1','ICELL8','Bulk RNA-seq',
'C1-HT','C1','ICELL8','Bulk RNA-seq','C1','ICELL8','Bulk RNA-seq','ICELL8',rep('Bulk RNA-seq',2))
da_ari_loss <- cbind.data.frame(ARI,Protocols)
da_ari_loss[,2] <- factor(da_ari_loss[,2],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq','WGS'))

p5_loss <-  ggplot(da_ari_loss,aes(x=Protocols,y=ARI,fill=Protocols)) + labs(x = '',y='') + geom_boxplot() + 
scale_fill_manual(values=col_type) + scale_y_continuous(limits=c(0,0.7), breaks = seq(0,0.7,by=0.1)) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5),
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.position = 'bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


## copykat ##

input <- './Results/10x/all_cell/copykat/v1/sens_spec_summary/'
input2 <- './Results/c1_fda/all_cell/copykat/v1/sens_spec_summary/'
input3 <- './Results/c1_llu/all_cell/copykat/v1/sens_spec_summary/'
input4 <- './Results/WaferGen/all_cell/copykat/v1/sens_spec_summary/'
input5 <- './Results/bulk/all_cell/copykat/v2/sens_spec_summary/'

## all cnvs ##

cnv_anno <- read.csv('./Anno/cyto/cnv_cyto_brca.csv')
cnv_10x <- read.csv(paste0(input,'cnv_allv1.csv'))
cnv_10x[,1] <- as.character(cnv_10x[,1])
cnv_10x[63,1] <- paste0(cnv_10x[63,1],'.2')
cnv_c1_ht <- read.csv(paste0(input2,'cnv_allv1.csv'))
cnv_c1_ht[,1] <- as.character(cnv_c1_ht[,1])
cnv_c1_ht[63,1] <- paste0(cnv_c1_ht[63,1],'.2')
cnv_c1 <- read.csv(paste0(input3,'cnv_allv1.csv'))
cnv_c1[,1] <- as.character(cnv_c1[,1])
cnv_c1[63,1] <- paste0(cnv_c1[63,1],'.2')
cnv_icell <- read.csv(paste0(input4,'cnv_allv1.csv'))
cnv_icell[,1] <- as.character(cnv_icell[,1])
cnv_icell[63,1] <- paste0(cnv_icell[63,1],'.2')
cnv_bulk <- read.csv(paste0(input5,'cnv_allv1.csv'))
cnv_bulk[,1] <- as.character(cnv_bulk[,1])
cnv_bulk[63,1] <- paste0(cnv_bulk[63,1],'.2')

cnv_10x <- cnv_10x[match(cnv_anno[,1],cnv_10x[,1]),]
cnv_c1_ht <- cnv_c1_ht[match(cnv_anno[,1],cnv_c1_ht[,1]),]
cnv_c1 <- cnv_c1[match(cnv_anno[,1],cnv_c1[,1]),]
cnv_icell <- cnv_icell[match(cnv_anno[,1],cnv_icell[,1]),]
cnv_bulk <- cnv_bulk[match(cnv_anno[,1],cnv_bulk[,1]),]

score_wgs <- as.integer(cnv_10x[,2])%%2*cnv_anno[,2]
score_10x <- as.integer(cnv_10x[,5])%%2*cnv_anno[,2]
score_c1_ht <- as.integer(cnv_c1_ht[,5])%%2*cnv_anno[,2]
score_c1 <- as.integer(cnv_c1[,5])%%2*cnv_anno[,2]
score_icell <- as.integer(cnv_icell[,5])%%2*cnv_anno[,2]
score_bulk <- as.integer(cnv_bulk[,5])%%2*cnv_anno[,2]

ari_wgs_10x <- adj.rand.index(score_wgs+1,score_10x+1)
ari_wgs_c1_ht <- adj.rand.index(score_wgs+1,score_c1_ht+1)
ari_wgs_c1 <- adj.rand.index(score_wgs+1,score_c1+1)
ari_wgs_icell <- adj.rand.index(score_wgs+1,score_icell+1)
ari_wgs_bulk <- adj.rand.index(score_wgs+1,score_bulk+1)
ari_10x_c1_ht <- adj.rand.index(score_10x+1,score_c1_ht+1)
ari_10x_c1 <- adj.rand.index(score_10x+1,score_c1+1)
ari_10x_icell <- adj.rand.index(score_10x+1,score_icell+1)
ari_10x_bulk <- adj.rand.index(score_10x+1,score_bulk+1)
ari_c1_ht_c1 <- adj.rand.index(score_c1_ht+1,score_c1+1)
ari_c1_ht_icell <- adj.rand.index(score_c1_ht+1,score_icell+1)
ari_c1_ht_bulk <- adj.rand.index(score_c1_ht+1,score_bulk+1)
ari_c1_icell <- adj.rand.index(score_c1+1,score_icell+1)
ari_c1_bulk <- adj.rand.index(score_c1+1,score_bulk+1)
ari_icell_bulk <- adj.rand.index(score_icell+1,score_bulk+1)

ARI <- c(ari_wgs_10x,ari_wgs_c1_ht,ari_wgs_c1,ari_wgs_icell,ari_wgs_bulk,ari_10x_c1_ht,ari_10x_c1,ari_10x_icell,ari_10x_bulk, ari_c1_ht_c1,ari_c1_ht_icell,ari_c1_ht_bulk,ari_c1_icell,ari_c1_bulk,ari_icell_bulk)
ARI <- rep(ARI,2)
Protocols <- c(rep('WGS',5),rep('10x',4),rep('C1-HT',3),rep('C1',2),'ICELL8','10x','C1-HT','C1','ICELL8','Bulk RNA-seq',
'C1-HT','C1','ICELL8','Bulk RNA-seq','C1','ICELL8','Bulk RNA-seq','ICELL8',rep('Bulk RNA-seq',2))
da_ari2 <- cbind.data.frame(ARI,Protocols)
da_ari2[,2] <- factor(da_ari[,2],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq','WGS'))

da_wgs <- cbind.data.frame(cnv_anno[,1],score_wgs,'WGS')
da_10x <- cbind.data.frame(cnv_anno[,1],score_10x,'10x')
da_c1_ht <- cbind.data.frame(cnv_anno[,1],score_c1_ht,'C1-HT')
da_c1 <- cbind.data.frame(cnv_anno[,1],score_c1,'C1')
da_icell <- cbind.data.frame(cnv_anno[,1],score_icell,'ICELL8')
da_bulk <- cbind.data.frame(cnv_anno[,1],score_bulk,'Bulk RNA-seq')

colnames(da_wgs) <- colnames(da_10x) <- colnames(da_c1_ht) <- colnames(da_c1) <- colnames(da_icell) <- colnames(da_bulk) <- c('Cytoband','Status','Protocols')
da2 <- rbind(da_wgs,da_10x,da_c1_ht,da_c1,da_icell,da_bulk)
da2[,1] <- factor(da2[,1],levels=as.character(cnv_anno[,1]))
levels(da2[,1])[27] <- levels(da2[,1])[26]
da2[,3] <- factor(da2[,3],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq','WGS'))

p6 <- ggplot(da2,aes(x=Cytoband,y=Status,fill=Protocols)) + labs(x = 'Cytoband',y='') + geom_bar(stat="identity") + 
scale_fill_manual(values=col_type) + theme(axis.text = element_text(size = 9,family='Times'), 
plot.title = element_text(size = 9,family='Times', hjust=0.5), axis.text.x = element_text(angle=90, hjust=1),
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "center"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))

p7 <-  ggplot(da_ari2,aes(x=Protocols,y=ARI,fill=Protocols)) + labs(x = '',y='') + geom_boxplot() + 
scale_fill_manual(values=col_type) + scale_y_continuous(limits=c(0.6,0.9), breaks = seq(0.6,0.9,by=0.1)) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5),
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.position = 'bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


## gain cnvs ##

input <- './Results/10x/all_cell/copykat/v1/sens_spec_summary/'
input2 <- './Results/c1_fda/all_cell/copykat/v1/sens_spec_summary/'
input3 <- './Results/c1_llu/all_cell/copykat/v1/sens_spec_summary/'
input4 <- './Results/WaferGen/all_cell/copykat/v1/sens_spec_summary/'
input5 <- './Results/bulk/all_cell/copykat/v2/sens_spec_summary/'

cnv_anno <- read.csv('./Anno/cyto/cnv_cyto_brca.csv')
cnv_10x <- read.csv(paste0(input,'cnv_gainv1.csv'))
cnv_10x[,1] <- as.character(cnv_10x[,1])
cnv_c1_ht <- read.csv(paste0(input2,'cnv_gainv1.csv'))
cnv_c1_ht[,1] <- as.character(cnv_c1_ht[,1])
cnv_c1 <- read.csv(paste0(input3,'cnv_gainv1.csv'))
cnv_c1[,1] <- as.character(cnv_c1[,1])
cnv_icell <- read.csv(paste0(input4,'cnv_gainv1.csv'))
cnv_icell[,1] <- as.character(cnv_icell[,1])
cnv_bulk <- read.csv(paste0(input5,'cnv_gainv1.csv'))
cnv_bulk[,1] <- as.character(cnv_bulk[,1])

index <- match(cnv_anno[,1],cnv_10x[,1])
cnv_anno <- cnv_anno[!is.na(index),]
cnv_10x <- cnv_10x[match(cnv_anno[,1],cnv_10x[,1]),]
cnv_c1_ht <- cnv_c1_ht[match(cnv_anno[,1],cnv_c1_ht[,1]),]
cnv_c1 <- cnv_c1[match(cnv_anno[,1],cnv_c1[,1]),]
cnv_icell <- cnv_icell[match(cnv_anno[,1],cnv_icell[,1]),]
cnv_bulk <- cnv_bulk[match(cnv_anno[,1],cnv_bulk[,1]),]

score_wgs <- as.integer(cnv_10x[,2])
score_10x <- as.integer(cnv_10x[,5])
score_c1_ht <- as.integer(cnv_c1_ht[,5])
score_c1 <- as.integer(cnv_c1[,5])
score_icell <- as.integer(cnv_icell[,5])
score_bulk <- as.integer(cnv_bulk[,5])

ari_wgs_10x <- adj.rand.index(score_wgs,score_10x)
ari_wgs_c1_ht <- adj.rand.index(score_wgs,score_c1_ht)
ari_wgs_c1 <- adj.rand.index(score_wgs,score_c1)
ari_wgs_icell <- adj.rand.index(score_wgs,score_icell)
ari_wgs_bulk <- adj.rand.index(score_wgs,score_bulk)
ari_10x_c1_ht <- adj.rand.index(score_10x,score_c1_ht)
ari_10x_c1 <- adj.rand.index(score_10x,score_c1)
ari_10x_icell <- adj.rand.index(score_10x,score_icell)
ari_10x_bulk <- adj.rand.index(score_10x,score_bulk)
ari_c1_ht_c1 <- adj.rand.index(score_c1_ht,score_c1)
ari_c1_ht_icell <- adj.rand.index(score_c1_ht,score_icell)
ari_c1_ht_bulk <- adj.rand.index(score_c1_ht,score_bulk)
ari_c1_icell <- adj.rand.index(score_c1,score_icell)
ari_c1_bulk <- adj.rand.index(score_c1,score_bulk)
ari_icell_bulk <- adj.rand.index(score_icell,score_bulk)

ARI <- c(ari_wgs_10x,ari_wgs_c1_ht,ari_wgs_c1,ari_wgs_icell,ari_wgs_bulk,ari_10x_c1_ht,ari_10x_c1,ari_10x_icell,ari_10x_bulk, ari_c1_ht_c1,ari_c1_ht_icell,ari_c1_ht_bulk,ari_c1_icell,ari_c1_bulk,ari_icell_bulk)
ARI <- rep(ARI,2)
Protocols <- c(rep('WGS',5),rep('10x',4),rep('C1-HT',3),rep('C1',2),'ICELL8','10x','C1-HT','C1','ICELL8','Bulk RNA-seq',
'C1-HT','C1','ICELL8','Bulk RNA-seq','C1','ICELL8','Bulk RNA-seq','ICELL8',rep('Bulk RNA-seq',2))
da_ari_gain2 <- cbind.data.frame(ARI,Protocols)
da_ari_gain2[,2] <- factor(da_ari_gain2[,2],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq','WGS'))

p7_gain <-  ggplot(da_ari_gain2,aes(x=Protocols,y=ARI,fill=Protocols)) + labs(x = '',y='') + geom_boxplot() + 
scale_fill_manual(values=col_type) + scale_y_continuous(limits=c(0,0.7), breaks = seq(0,0.7,by=0.1)) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5),
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.position = 'bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))



## loss cnvs ##

input <- './Results/10x/all_cell/copykat/v1/sens_spec_summary/'
input2 <- './Results/c1_fda/all_cell/copykat/v1/sens_spec_summary/'
input3 <- './Results/c1_llu/all_cell/copykat/v1/sens_spec_summary/'
input4 <- './Results/WaferGen/all_cell/copykat/v1/sens_spec_summary/'
input5 <- './Results/bulk/all_cell/copykat/v2/sens_spec_summary/'

cnv_anno <- read.csv('./Anno/cyto/cnv_cyto_brca.csv')
cnv_10x <- read.csv(paste0(input,'cnv_lossv1.csv'))
cnv_10x[,1] <- as.character(cnv_10x[,1])
cnv_10x[30,1] <- paste0(cnv_10x[30,1],'.2')
cnv_c1_ht <- read.csv(paste0(input2,'cnv_lossv1.csv'))
cnv_c1_ht[,1] <- as.character(cnv_c1_ht[,1])
cnv_c1_ht[30,1] <- paste0(cnv_c1_ht[30,1],'.2')
cnv_c1 <- read.csv(paste0(input3,'cnv_lossv1.csv'))
cnv_c1[,1] <- as.character(cnv_c1[,1])
cnv_c1[30,1] <- paste0(cnv_c1[30,1],'.2')
cnv_icell <- read.csv(paste0(input4,'cnv_lossv1.csv'))
cnv_icell[,1] <- as.character(cnv_icell[,1])
cnv_icell[30,1] <- paste0(cnv_icell[30,1],'.2')
cnv_bulk <- read.csv(paste0(input5,'cnv_lossv1.csv'))
cnv_bulk[,1] <- as.character(cnv_bulk[,1])
cnv_bulk[30,1] <- paste0(cnv_bulk[30,1],'.2')

index <- match(cnv_anno[,1],cnv_10x[,1])
cnv_anno <- cnv_anno[!is.na(index),]
cnv_10x <- cnv_10x[match(cnv_anno[,1],cnv_10x[,1]),]
cnv_c1_ht <- cnv_c1_ht[match(cnv_anno[,1],cnv_c1_ht[,1]),]
cnv_c1 <- cnv_c1[match(cnv_anno[,1],cnv_c1[,1]),]
cnv_icell <- cnv_icell[match(cnv_anno[,1],cnv_icell[,1]),]
cnv_bulk <- cnv_bulk[match(cnv_anno[,1],cnv_bulk[,1]),]

score_wgs <- as.integer(cnv_10x[,2])
score_10x <- as.integer(cnv_10x[,5])
score_c1_ht <- as.integer(cnv_c1_ht[,5])
score_c1 <- as.integer(cnv_c1[,5])
score_icell <- as.integer(cnv_icell[,5])
score_bulk <- as.integer(cnv_bulk[,5])

ari_wgs_10x <- adj.rand.index(score_wgs,score_10x)
ari_wgs_c1_ht <- adj.rand.index(score_wgs,score_c1_ht)
ari_wgs_c1 <- adj.rand.index(score_wgs,score_c1)
ari_wgs_icell <- adj.rand.index(score_wgs,score_icell)
ari_wgs_bulk <- adj.rand.index(score_wgs,score_bulk)
ari_10x_c1_ht <- adj.rand.index(score_10x,score_c1_ht)
ari_10x_c1 <- adj.rand.index(score_10x,score_c1)
ari_10x_icell <- adj.rand.index(score_10x,score_icell)
ari_10x_bulk <- adj.rand.index(score_10x,score_bulk)
ari_c1_ht_c1 <- adj.rand.index(score_c1_ht,score_c1)
ari_c1_ht_icell <- adj.rand.index(score_c1_ht,score_icell)
ari_c1_ht_bulk <- adj.rand.index(score_c1_ht,score_bulk)
ari_c1_icell <- adj.rand.index(score_c1,score_icell)
ari_c1_bulk <- adj.rand.index(score_c1,score_bulk)
ari_icell_bulk <- adj.rand.index(score_icell,score_bulk)

ARI <- c(ari_wgs_10x,ari_wgs_c1_ht,ari_wgs_c1,ari_wgs_icell,ari_wgs_bulk,ari_10x_c1_ht,ari_10x_c1,ari_10x_icell,ari_10x_bulk, ari_c1_ht_c1,ari_c1_ht_icell,ari_c1_ht_bulk,ari_c1_icell,ari_c1_bulk,ari_icell_bulk)
ARI <- rep(ARI,2)
Protocols <- c(rep('WGS',5),rep('10x',4),rep('C1-HT',3),rep('C1',2),'ICELL8','10x','C1-HT','C1','ICELL8','Bulk RNA-seq',
'C1-HT','C1','ICELL8','Bulk RNA-seq','C1','ICELL8','Bulk RNA-seq','ICELL8',rep('Bulk RNA-seq',2))
da_ari_loss2 <- cbind.data.frame(ARI,Protocols)
da_ari_loss2[,2] <- factor(da_ari_loss2[,2],levels=c('C1','ICELL8','C1-HT','10x','Bulk RNA-seq','WGS'))

p7_loss <-  ggplot(da_ari_loss2,aes(x=Protocols,y=ARI,fill=Protocols)) + labs(x = '',y='') + geom_boxplot() + 
scale_fill_manual(values=col_type) + scale_y_continuous(limits=c(0,0.7), breaks = seq(0,0.7,by=0.1)) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5),
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.position = 'bottom',
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


##########################################
## ROC for casper and copykat ##
##########################################

input <- paste0('./Results/sens_spec_eva/cyto/all_cell/ref_v1/v1/')
input2 <- paste0('./Results/sens_spec_eva/cyto/all_cell/ref_v2/v1/')
ref_type <- c('scRNA-seq','Bulk RNA-seq','GTex')

method <- c('CaSpER','CopyKAT')
method2 <- c('casper','copykat')
da <- c()
da2 <- c()
for(i in 1:length(method))
{
temp <- read.csv(paste0(input,'sens_spec_all_',method2[i],'.csv'))
colnames(temp)[1:2] <- c('Sensitivity','Specificity')
temp <- cbind(temp,method[i])
da <- rbind(da,temp)

temp2 <- read.csv(paste0(input2,'sens_spec_all_',method2[i],'.csv'))
colnames(temp2)[1:2] <- c('Sensitivity','Specificity')
temp2 <- temp2[temp2[,4]=='Bulk',]
temp2 <- cbind(temp2,method[i])
da2 <- rbind(da2,temp2)
}
colnames(da)[5] <- colnames(da2)[5] <- 'Methods'
da <- rbind(da,da2)
da[,2] <- 1-da[,2]
levels(da[,4])[5] <- 'Bulk RNA-seq'
da[,4] <- factor(da[,4],levels=c('C1','ICELL8','C1-HT','10x', 'Bulk RNA-seq'))

p8 <- ggplot(da,aes(x=Specificity,y=Sensitivity,col=Protocols,linetype=Methods)) + labs(x = '1-Specificity') + geom_path() +
scale_colour_manual(values=col_type) + scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


p1_leg <- get_legend(p1)
p1 <- p1 + theme(legend.position='none')
p2 <- p2 + theme(legend.position='none')
p3 <- p3 + theme(legend.position='none')
p4_leg <- get_legend(p4)
p4 <- p4 + theme(legend.position='none')
p5 <- p5 + theme(legend.position='none')
p6 <- p6 + theme(legend.position='none')
p7 <- p7 + theme(legend.position='none')
p8_leg <- get_legend(p8)
p8 <- p8 + theme(legend.position='none')

p9 <- plot_grid(p1,p2,p3,p1_leg, ncol=4, labels=c('a','b','c'), rel_widths=c(2,2,2,1.1))
p10 <- plot_grid(plot_grid(p4,p6,nrow=2,labels=c('d','e')),p4_leg, rel_widths=c(6,1.1))
p11 <- plot_grid(p5, p7, p8, p8_leg, ncol=4, labels=c('f','g','h'), rel_widths=c(2,2,2,1.1))

p12 <- plot_grid(p9,p10,p11,nrow=3,rel_heights=c(1,1.5,1))

jpeg('./Manuscript/figure_table/figure2_v7.jpeg',width=12,height=12,res=300,units='in')
p12
dev.off()

pdf('./Manuscript/figure_table/figure2_v7.pdf',width=12,height=12)
p12
dev.off()






























