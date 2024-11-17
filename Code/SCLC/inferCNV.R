#########################
## generate raw counts ##
#########################

load('./Results/SCLC/Star_featureCounts/Ensembl_GRCh38/Gtf/gene_counts.rdata')
t_counts <- gene_counts
t_counts <- as.matrix(t_counts)
colnames(t_counts) <- paste(colnames(t_counts),1,sep="_")

load('./Anno/GTEx/lung_gene_counts.Rdata')
n_bulk <- gene_counts
overlap_gene <- intersect(rownames(t_counts),rownames(n_bulk))
t_counts <- t_counts[overlap_gene,]
n_bulk <- n_bulk[overlap_gene,]
colnames(n_bulk) <- paste(colnames(n_bulk),2,sep='_')
raw_counts <- cbind(t_counts,n_bulk)
sample_anno <- cbind(colnames(raw_counts),c(rep('Tumor',ncol(t_counts)),rep('Normal',ncol(n_bulk))))

output <- './Results/SCLC/infercnv/'
dir.create(output,recursive=T)
write.table(as.matrix(raw_counts),file=paste0(output,'raw_counts.txt'),sep='\t',row.names=T,col.names=T,quote=F)
write.table(sample_anno,file=paste0(output,'sample_anno.txt'),sep='\t',row.names=F,col.names=F,quote=F)


##############
## infercnv ##
##############

library(infercnv)
library(biomaRt)

input <- './Results/SCLC/infercnv/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts.txt.gz')
sample_input <- paste0(input,'/sample_anno.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- paste0(input,'/results')
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 20)






