
library(ggpubr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(RColorBrewer)

col_type <- brewer.pal(9,'Set1')


#####################################################
## reference for infercnv, casper, scicnv, copykat ##
#####################################################

input <- paste0('./Results/sens_spec_eva/cyto/all_cell/ref_v',1:3,'/v1/')
ref_type <- c('scRNA-seq','Bulk RNA-seq','GTex')

da <- c()
for(i in 1:length(input))
{
infercnv_summary <- read.csv(paste0(input[i],'sens_spec_all_infercnv.csv'))
infercnv_summary <- cbind.data.frame(infercnv_summary[infercnv_summary[,3]=='sc_0.1',-3],'inferCNV')
casper_summary <- read.csv(paste0(input[i],'sens_spec_all_casper.csv'))
casper_summary <- cbind.data.frame(casper_summary[casper_summary[,3]=='sc_0.1',-3],'CaSpER')
scicnv_summary <- read.csv(paste0(input[i],'sens_spec_all_scicnv.csv'))
scicnv_summary <- cbind.data.frame(scicnv_summary[scicnv_summary[,3]=='sc_0.1',-3],'sciCNV')
copykat_summary <- read.csv(paste0(input[i],'sens_spec_all_copykat.csv'))
copykat_summary <- cbind.data.frame(copykat_summary[copykat_summary[,3]=='sc_0.1',-3],'CopyKAT')
colnames(infercnv_summary)[-3] <- colnames(casper_summary)[-3] <- colnames(scicnv_summary)[-3] <- colnames(copykat_summary)[-3] <- c('Sensitivity','Specificity','Methods')
temp <- rbind(infercnv_summary,casper_summary,scicnv_summary,copykat_summary)
temp <- cbind(temp,ref_type[i])
da <- rbind(da,temp)
}
colnames(da)[5] <-'References'
da <- da[da[,3]!='Bulk',]
da[,3] <- factor(da[,3],levels=c('C1','ICELL8','C1-HT','10x'))
da_sens <- cbind.data.frame(da[,-2],'Sensitivity')
da_spec <- cbind.data.frame(da[,-1],'Specificity')
colnames(da_sens)[c(1,5)] <- colnames(da_spec)[c(1,5)] <- c('Value','Metrics')
da2 <- rbind(da_sens,da_spec)

p1 <- ggplot(da2,aes(x=Methods,y=Value,fill=References)) + labs(x = 'Methods',y='') + geom_boxplot() + 
facet_grid(cols=vars(Metrics)) + scale_fill_manual(values=col_type) + 
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.text.x = element_text(angle=20, hjust=1), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))



##########################################
## reference ROC for casper and copykat, infercnv, and scicnv ##
##########################################

## casper ##

input <- paste0('./Results/sens_spec_eva/cyto/all_cell/ref_v',1:3,'/v1/')
ref_type <- c('scRNA-seq','Bulk RNA-seq','GTex')

da <- c()
for(i in 1:length(input))
{
temp <- read.csv(paste0(input[i],'sens_spec_all_casper.csv'))
colnames(temp)[1:2] <- c('Sensitivity','Specificity')
temp <- cbind(temp,ref_type[i])
da <- rbind(da,temp)
}
colnames(da)[5] <-'References'
da[,2] <- 1-da[,2]
da[,4] <- factor(da[,4],levels=c('C1','ICELL8','C1-HT','10x'))
da <- da[!is.na(da[,4]),]
da <- da[da[,3]!='sc_0',]

p2 <- ggplot(da,aes(x=Specificity,y=Sensitivity,col=References,linetype=Protocols)) + labs(x = '1-Specificity') + geom_path() +
scale_colour_manual(values=col_type) + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + geom_abline(slope = 1) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


## copykat ##

input <- paste0('./Results/sens_spec_eva/cyto/all_cell/ref_v',1:3,'/v1/')
ref_type <- c('scRNA-seq','Bulk RNA-seq','GTex')

da2 <- c()
for(i in 1:length(input))
{
temp <- read.csv(paste0(input[i],'sens_spec_all_copykat.csv'))
colnames(temp)[1:2] <- c('Sensitivity','Specificity')
temp <- cbind(temp,ref_type[i])
da2 <- rbind(da2,temp)
}
colnames(da2)[5] <-'References'
da2[,2] <- 1-da2[,2]
da2[,4] <- factor(da2[,4],levels=c('C1','ICELL8','C1-HT','10x'))
da2 <- da2[!is.na(da2[,4]),]
da2 <- da2[da2[,3]!='sc_0',]

p3 <- ggplot(da2,aes(x=Specificity,y=Sensitivity,col=References,linetype=Protocols)) + labs(x = '1-Specificity') + geom_path() +
scale_colour_manual(values=col_type)  + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + geom_abline(slope = 1) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))



## infercnv ##

input <- paste0('./Results/sens_spec_eva/cyto/all_cell/ref_v',1:3,'/v1/')
ref_type <- c('scRNA-seq','Bulk RNA-seq','GTex')

da2 <- c()
for(i in 1:length(input))
{
temp <- read.csv(paste0(input[i],'sens_spec_all_infercnv.csv'))
colnames(temp)[1:2] <- c('Sensitivity','Specificity')
temp <- cbind(temp,ref_type[i])
da2 <- rbind(da2,temp)
}
colnames(da2)[5] <-'References'
da2[,2] <- 1-da2[,2]
da2[,4] <- factor(da2[,4],levels=c('C1','ICELL8','C1-HT','10x'))
da2 <- da2[!is.na(da2[,4]),]

p4 <- ggplot(da2,aes(x=Specificity,y=Sensitivity,col=References,linetype=Protocols)) + labs(x = '1-Specificity') + geom_path() +
scale_colour_manual(values=col_type)  + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + geom_abline(slope = 1) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


## scicnv ##

input <- paste0('./Results/sens_spec_eva/cyto/all_cell/ref_v',1:3,'/v1/')
ref_type <- c('scRNA-seq','Bulk RNA-seq','GTex')

da2 <- c()
for(i in 1:length(input))
{
temp <- read.csv(paste0(input[i],'sens_spec_all_scicnv.csv'))
colnames(temp)[1:2] <- c('Sensitivity','Specificity')
temp <- cbind(temp,ref_type[i])
da2 <- rbind(da2,temp)
}
colnames(da2)[5] <-'References'
da2[,2] <- 1-da2[,2]
da2[,4] <- factor(da2[,4],levels=c('C1','ICELL8','C1-HT','10x'))
da2 <- da2[!is.na(da2[,4]),]

p5 <- ggplot(da2,aes(x=Specificity,y=Sensitivity,col=References,linetype=Protocols)) + labs(x = '1-Specificity') + geom_path() +
scale_colour_manual(values=col_type)  + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + geom_abline(slope = 1) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


## compare bulk rnaseq between casper and copykat ## 

da_bulk <- rbind(cbind(da[da[,5]=='Bulk RNA-seq',],Methods='CaSpER'),cbind(da2[da2[,5]=='Bulk RNA-seq',],Methods='CopyKAT'))

p6 <- ggplot(da_bulk,aes(x=Specificity,y=Sensitivity,col=Methods,linetype=Protocols)) + labs(x = '1-Specificity') + geom_path() +
scale_colour_manual(values=col_type)  + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + geom_abline(slope = 1) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))


da_gtex <- rbind(cbind(da[da[,5]=='GTex',],Methods='CaSpER'),cbind(da2[da2[,5]=='GTex',],Methods='CopyKAT'))

p7 <- ggplot(da_gtex,aes(x=Specificity,y=Sensitivity,col=Methods,linetype=Protocols)) + labs(x = '1-Specificity') + geom_path() +
scale_colour_manual(values=col_type) + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + geom_abline(slope = 1) +
theme(axis.text = element_text(size = 9,family='Times'), plot.title = element_text(size = 9,family='Times', hjust=0.5), 
axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), legend.justification = c("left", "top"),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill='white',colour='black'))



##############################################
## read length and read depth effect ##
##############################################

## get all data ##

input <- './Results/sens_spec_eva/cyto/read_len_subsampling_SE_subCell/cell_1/'
rl <- dir(input)
rd <- dir(paste0(input,rl[1]))

cnv_mds <- c('CaSpER','CopyKAT')
da <- c()
for(i in 1:length(rl))
{
	for(j in 1:length(rd))
	{
	temp_input <- paste0(input,rl[i],'/',rd[j],'/ref_v1/v1/')
	files <- dir(temp_input)[1:2]
	for(k in 1:length(files))
	{
	temp <- read.csv(paste0(temp_input,files[k]))
	temp <- cbind.data.frame(temp,rl[i],rd[j],cnv_mds[k])
	colnames(temp)[-(3:4)] <- c('Sensitivity','Specificity','Read_length','Read_depth','Methods')
	da <- rbind(da,temp)
	}
	}
}
da[,2] <- 1-da[,2]
da[,6] <- factor(da[,6],levels=rd[c(5,1,3,4,2)])
da[,5] <- factor(da[,5],levels=rl[c(4,5,1:3)])
levels(da[,4])[4] <- 'C1'

protocols <- names(table(da[,4]))
protocols <- protocols[c(4,5,1,2,3)]
p <- vector('list',length(protocols)*length(cnv_mds))

for(i in 1:length(protocols))
{
	for(j in 1:length(cnv_mds))
	{
	temp_da <- da[da[,4]==protocols[i],]
	temp_da <- temp_da[temp_da[,7]==cnv_mds[j],]
	if(i == 3)
	{temp_da <- temp_da[temp_da[,5]=='100bp',]}
	else if(i == 4)
	{temp_da <- temp_da[temp_da[,5]=='125bp',]}
	else
	{temp_da <- temp_da[temp_da[,5]=='150bp',]}
	p[[(i-1)*2+j]] <- ggplot(temp_da,aes(x = Specificity, y = Sensitivity, color=Read_depth)) + 
	labs(x = '1-Specificity') + geom_path() + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
	scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
	scale_colour_manual(values=col_type) + labs(title = paste('Protocol:',protocols[i],'  Method:',cnv_mds[j]), 
	x = '1-Specificity', color='Read depth') + theme(axis.text = element_text(size = 9,family='Times'), 
	plot.title = element_text(size = 9,family='Times',hjust=0.5), 
	axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
	legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), 
	panel.grid.major = element_blank(), legend.justification = c("left", "top"), panel.grid.minor = element_blank(), 
	panel.background = element_rect(fill='white',colour='black'))

	}
}

p_legend <- get_legend(p[[1]])
for(i in 1:length(p))
{
	p[[i]] <- p[[i]] + theme(legend.position='none')
}

p6_legend <- get_legend(p6)
p6 <- p6 + theme(legend.position='none')
p7 <- p7 + theme(legend.position='none')
p8 <- plot_grid(p1,p6,p7,p6_legend,ncol=4,rel_widths=c(4,2,2,0.7),labels=c('a','b','c'))
p9 <- plot_grid(p[[1]],p[[3]],p[[5]],p[[7]],p[[2]],p[[4]],p[[6]],p[[8]],nrow=2, labels=c('d','e','f','g','h','i','j','k'))
p10 <- plot_grid(p9,p_legend,nrow=1,rel_widths=c(8,0.7))
p11 <- plot_grid(p8,p10,nrow=2,rel_heights=c(1,2))

jpeg('./Manuscript/figure_table/figure3_v3.jpeg',width=12,height=12,res=300,units='in')
p11
dev.off()

pdf('./Manuscript/figure_table/figure3_v3.pdf',width=12,height=12)
p11
dev.off()


p2_legend <- get_legend(p2)
p2 <- p2 + theme(legend.position='none')
p3 <- p3 + theme(legend.position='none')
p8 <- plot_grid(p1,p2,p3,p2_legend,ncol=4,rel_widths=c(4,2,2,0.8),labels=c('a','b','c'))
p9 <- plot_grid(p[[1]],p[[3]],p[[5]],p[[7]],p[[2]],p[[4]],p[[6]],p[[8]],nrow=2, labels=c('d','e','f','g','h','i','j','k'))
p10 <- plot_grid(p9,p_legend,nrow=1,rel_widths=c(8,0.7))
p11 <- plot_grid(p8,p10,nrow=2,rel_heights=c(1,2))

jpeg('./Manuscript/figure_table/figure3_v4.jpeg',width=12,height=12,res=300,units='in')
p11
dev.off()

pdf('./Manuscript/figure_table/figure3_v4.pdf',width=12,height=12)
p11
dev.off()






protocols <- names(table(da[,4]))
protocols <- protocols[c(4,5,1,2)]
cnv_mds <- c('CaSpER','CopyKAT')
n <- 1
p_rl <- c()
for(j in 1:length(cnv_mds))
{
	for(i in 1:length(protocols))
	{
	temp_da <- da[da[,4]==protocols[i],]
	temp_da <- temp_da[temp_da[,7]==cnv_mds[j],]
	temp_rl <- table(temp_da[,5])
	temp_rl <- temp_rl[temp_rl!=0]
	rl <- names(temp_rl)
	for(k in 1:length(rl))
	{
	temp_da2 <- temp_da[temp_da[,5]==rl[k],]
	p_rl[[n]] <- ggplot(temp_da2,aes(x = Specificity, y = Sensitivity, color=Read_depth)) + 
	labs(x = '1-Specificity') + geom_path() + scale_x_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
	scale_y_continuous(limits=c(0,1), breaks = seq(0,1,by=0.1)) + 
	scale_colour_manual(values=col_type) + labs(title = paste('Protocol:',protocols[i],' Method:',cnv_mds[j],' rl:',rl[k]), 
	x = '1-Specificity', color='Read depth') + theme(axis.text = element_text(size = 9,family='Times'), 
	plot.title = element_text(size = 9,family='Times',hjust=0.5), 
	axis.title.x = element_text(size = 9,family='Times'), axis.title.y = element_text(size = 9,family='Times'), 
	legend.title = element_text(size = 9,family='Times'), legend.text = element_text(size = 9,family='Times'), 
	panel.grid.major = element_blank(), legend.justification = c("left", "top"), panel.grid.minor = element_blank(), 
	panel.background = element_rect(fill='white',colour='black'))
	n <- n+1
	}
	}
}

p_rl_legend <- get_legend(p_rl[[1]])
for(i in 1:length(p_rl))
{
	p_rl[[i]] <- p_rl[[i]] + theme(legend.position='none')
}

p_rl_casper1 <- plot_grid(p_rl[[1]],p_rl[[2]],p_rl[[3]],p_rl[[4]],p_rl[[5]],nrow=1,labels=c('a','b','c','d','e'))
p_rl_casper2 <- plot_grid(p_rl[[6]],p_rl[[7]],p_rl[[8]],p_rl[[9]],p_rl[[10]],nrow=1,labels=c('f','g','h','i','j'))
p_rl_casper3 <- plot_grid(p_rl[[14]],p_rl[[15]],p_rl[[16]],p_rl[[17]],nrow=1,labels=c('k','l','m','n','o'))
p_rl_casper4 <- plot_grid(p_rl[[11]],p_rl[[12]],p_rl[[13]],p_rl_legend,nrow=1,labels=c('p','q'))
p_rl_casper5 <- plot_grid(p_rl_casper1,p_rl_casper2,p_rl_casper3,p_rl_casper4,nrow=4)


p_rl_copykat1 <- plot_grid(p_rl[[18]],p_rl[[19]],p_rl[[20]],p_rl[[21]],p_rl[[22]],nrow=1,labels=c('a','b','c','d','e'))
p_rl_copykat2 <- plot_grid(p_rl[[23]],p_rl[[24]],p_rl[[25]],p_rl[[26]],p_rl[[27]],nrow=1,labels=c('f','g','h','i','j'))
p_rl_copykat3 <- plot_grid(p_rl[[31]],p_rl[[32]],p_rl[[33]],p_rl[[34]],nrow=1,labels=c('k','l','m','n','o'))
p_rl_copykat4 <- plot_grid(p_rl[[28]],p_rl[[29]],p_rl[[30]],p_rl_legend,nrow=1,labels=c('p','q'))
p_rl_copykat5 <- plot_grid(p_rl_copykat1,p_rl_copykat2,p_rl_copykat3,p_rl_copykat4,nrow=4)


jpeg('./Manuscript/figure_table/sup_figure/sup_fig3_v2.jpeg',width=15,height=12,res=300,units='in')
p_rl_casper5
dev.off()

pdf('./Manuscript/figure_table/sup_figure/sup_fig3_v2.pdf',width=15,height=12)
p_rl_casper5
dev.off()


jpeg('./Manuscript/figure_table/sup_figure/sup_fig4_v2.jpeg',width=15,height=12,res=300,units='in')
p_rl_copykat5
dev.off()

pdf('./Manuscript/figure_table/sup_figure/sup_fig4_v2.pdf',width=15,height=12)
p_rl_copykat5
dev.off()


