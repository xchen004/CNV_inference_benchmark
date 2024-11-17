
#################
## 10x ##
#################

library(infercnv)
library(biomaRt)

input <- './Data/10x/gene_counts/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts_v3.txt.gz')
sample_input <- paste0(input,'/sample_anno_v3.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- './Results/10x/all_cell/infercnv/v3/subcls/results/'
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
#infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=output, cluster_by_groups=T, denoise=T, HMM=T, num_threads = 16)

infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 16)


#################
## c1 ##
#################

library(infercnv)
library(biomaRt)

input <- './Data/C1/gene_counts/read_len_SE/150bp/'
input2 <- './Anno/'
raw_input <- paste0(input,'raw_counts_v3.txt.gz')
sample_input <- paste0(input,'sample_anno_v3.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- './Results/c1_llu/all_cell/infercnv/v3/subcls/results/'
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
#infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=output, cluster_by_groups=T, denoise=T, HMM=T, num_threads = 16)

infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 16)


#################
## c1-fda ##
#################

library(infercnv)
library(biomaRt)

input <- './Data/C1_ceber/gene_counts/'
input2 <- './Anno/'
raw_input <- paste0(input,'raw_counts_v3.txt.gz')
sample_input <- paste0(input,'sample_anno_v3.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- './Results/c1_fda/all_cell/infercnv/v3/subcls/results/'
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
#infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=output, cluster_by_groups=T, denoise=T, HMM=T, num_threads = 16)

infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=output, denoise=T, HMM=T, HMM_type='i3', analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 16)



#################
## icell8 ##
#################

library(infercnv)
library(biomaRt)

input <- './Data/WaferGen/gene_counts/read_len_SE/150bp/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts_v3.txt.gz')
sample_input <- paste0(input,'/sample_anno_v3.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- './Results/WaferGen/all_cell/infercnv/v3/subcls/results/'
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
#infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, cluster_by_groups=T, denoise=T, HMM=T, num_threads = 16)

infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 16)



#################
## bulk rnaseq ##
#################

library(infercnv)
library(biomaRt)

input <- './Data/RNAseq/merge/Star/NCBI_GRCh38/Gtf/HCC1395/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts_v3.txt.gz')
sample_input <- paste0(input,'/sample_anno_v3.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- './Results/bulk/infercnv/v3/results/'
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, cluster_by_groups=T, denoise=T, HMM=T, HMM_type='i3', num_threads = 16)



