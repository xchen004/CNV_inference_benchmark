
#################
## 10x_3cl_310 ##
#################


library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)
library(RColorBrewer)

input <- './Results/Tian/infercnv/10x_3cl_310/results/'
output <- './Results/Tian/infercnv/10x_3cl_310/summary/'
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_310.csv',row.names=1)
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,0,1)))
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
col <- brewer.pal(8,"Set1")
cls_all <- cbind(cell_cls,cls)
dir.create(output,recursive=T)
cls3_table <- table(cell_cls[,2],cls[,1])
cls4_table <- table(cell_cls[,2],cls[,2])
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)
write.csv(cls3_table,file=paste0(output,'cls3_table.csv'),quote=F)
write.csv(cls4_table,file=paste0(output,'cls4_table.csv'),quote=F)

bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(dend)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=3, value=col[1:3]) %>% plot
colored_bars(colors=bar_col,dend=dend,y_shift=-30,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975'),col=col[1:3],pch=15,bty='n')
dev.off()

## cluster accuracy ##

H1975:	336/346=97.11%
H2228:	325/327=99.39%
HCC827:	618/927=66.67%




#################
## 10x_3cl_220 ##
#################


library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)
library(RColorBrewer)

input <- './Results/Tian/infercnv/10x_3cl_220/results/'
output <- './Results/Tian/infercnv/10x_3cl_220/summary/'
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_220.csv',row.names=1)
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
cell_cls <- cell_cls[rownames(cls),]
colnames(cls) <- paste0('pred_cls',colnames(cls))
col <- brewer.pal(8,"Set1")
cls_all <- cbind(cell_cls,cls)
dir.create(output,recursive=T)
cls3_table <- table(cell_cls[,2],cls[,1])
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)
write.csv(cls3_table,file=paste0(output,'cls3_table.csv'),quote=F)

bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(dend)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=3, value=col[1:3]) %>% plot
colored_bars(colors=bar_col,dend=dend,y_shift=-30,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975'),col=col[1:3],pch=15,bty='n')
dev.off()

## cluster accuracy ##

H1975:	336/343=97.96%
H2228:	325/325=100%
HCC827:	296/296=100%







#################
## 10x_5cl_310 ##
#################


library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)
library(RColorBrewer)

input <- './Results/Tian/infercnv/10x_5cl_310/results/'
output <- './Results/Tian/infercnv/10x_5cl_310/summary/'
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_5cl_310.csv',row.names=1)
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,3,4,0,1)))
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,5:14)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
col <- brewer.pal(8,"Set1")
cls_all <- cbind(cell_cls,cls)
dir.create(output,recursive=T)
cls5_table <- table(cell_cls[,2],cls[,1])
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)
write.csv(cls5_table,file=paste0(output,'cls5_table.csv'),quote=F)

bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(dend)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=5, value=col[c(4,1,5,3,2)]) %>% plot
colored_bars(colors=bar_col,dend=dend,y_shift=-30,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975','A549','H838'),col=col[1:5],pch=15,bty='n')
dev.off()

## cluster accuracy ##

H1975:	425/456=93.20%
H2228:	751/794=94.58%
HCC827:	600/600=100%
A549:	1257/1265=99.37%
H838:	836/919=90.97%





#################
## 10x_5cl_201 ##
#################


library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)
library(RColorBrewer)

input <- './Results/Tian/infercnv/10x_5cl_201/results/'
output <- './Results/Tian/infercnv/10x_5cl_201/summary/'
cell_cls <- read.csv('./Results/Tian/cluster/sc_10x_5cl_201.csv',row.names=1)
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,3,4,0,1)))
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,5:14)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
col <- brewer.pal(8,"Set1")
cls_all <- cbind(cell_cls,cls)
dir.create(output,recursive=T)
cls5_table <- table(cell_cls[,2],cls[,1])
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)
write.csv(cls5_table,file=paste0(output,'cls5_table.csv'),quote=F)

bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(dend)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=5, value=col[c(1,4,5,2,3)]) %>% plot
colored_bars(colors=bar_col,dend=dend,y_shift=-30,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975','A549','H838'),col=col[1:5],pch=15,bty='n')
dev.off()

## cluster accuracy ##

H1975:	423/429=98.60%
H2228:	743/746=99.60%
HCC827:	568/576=98.61%
A549:	1187/1191=99.66%
H838:	852/852=100%



#############
## celseq2 ##
#############


library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)
library(RColorBrewer)

input <- './Results/Tian/infercnv/celseq2/results/'
output <- './Results/Tian/infercnv/celseq2/summary/'
cell_cls <- read.csv('./Results/Tian/cluster/celseq2.csv',row.names=1)
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(1,2,0)))
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
colnames(cls) <- paste0('pred_cls',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
col <- brewer.pal(8,"Set1")
cls_all <- cbind(cell_cls,cls)
dir.create(output,recursive=T)
cls3_table <- table(cell_cls[,2],cls[,1])
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)
write.csv(cls3_table,file=paste0(output,'cls3_table.csv'),quote=F)

bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(dend)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=3, value=col[1:3]) %>% plot
colored_bars(colors=bar_col,dend=dend,y_shift=-30,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975'),col=col[1:3],pch=15,bty='n')
dev.off()


## cluster accuracy ##

H1975:	89/91=97.8%
H2228:	75/76=98.68%
HCC827:	75/75=100%



#############
## dropseq ##
#############

library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(dendextend)
library(RColorBrewer)

input <- './Results/Tian/infercnv/dropseq/results/'
output <- './Results/Tian/infercnv/dropseq/summary/'
cell_cls <- read.csv('./Results/Tian/cluster/dropseq.csv',row.names=1)
cell_cls[,1] <- as.integer(factor(cell_cls[,1],levels=c(2,1,0)))
rownames(cell_cls) <- paste0(rownames(cell_cls),'_1')

dend <- read.newick(paste0(input,'infercnv.observations_dendrogram.txt')) %>% as.hclust
cls <- cutree(dend,3:12)
colnames(cls) <- paste0('pred_cls_',colnames(cls))
cell_cls <- cell_cls[rownames(cls),]
col <- brewer.pal(8,"Set1")
cls_all <- cbind(cell_cls,cls)
dir.create(output,recursive=T)
cls3_table <- table(cell_cls[,2],cls[,1])
write.csv(cls_all,file=paste0(output,'cls_all.csv'),quote=F)
write.csv(cls3_table,file=paste0(output,'cls3_table.csv'),quote=F)

bar_col <- col[cell_cls[,1]]
dend <- as.dendrogram(dend)

jpeg(paste0(output,'cls.jpeg'),units='in',width=8,height=8,res=300)
dend %>% set('labels','') %>% set('branches_k_color',k=3, value=col[1:3]) %>% plot
colored_bars(colors=bar_col,dend=dend,y_shift=-30,rowLabels='Cell lines')
legend('topright',c('H2228','HCC827','H1975'),col=col[1:3],pch=15,bty='n')
dev.off()


## cluster accuracy ##

H1975:	82/84=97.62%
H2228:	54/56=96.43%
HCC827:	61/61=100%














