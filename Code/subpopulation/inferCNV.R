
#################
## 10x_3cl_310 ##
#################

library(infercnv)
library(biomaRt)

input <- './Results/Tian/infercnv/10x_3cl_310/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts.txt.gz')
sample_input <- paste0(input,'/sample_anno.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- paste0(input,'/results')
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 10)



#################
## 10x_3cl_220 ##
#################

library(infercnv)
library(biomaRt)

input <- './Results/Tian/infercnv/10x_3cl_220/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts.txt.gz')
sample_input <- paste0(input,'/sample_anno.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- paste0(input,'/results')
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 10)



#################
## 10x_5cl_310 ##
#################

library(infercnv)
library(biomaRt)

input <- './Results/Tian/infercnv/10x_5cl_310/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts.txt.gz')
sample_input <- paste0(input,'/sample_anno.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- paste0(input,'/results')
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 10)



#################
## 10x_5cl_201 ##
#################

library(infercnv)
library(biomaRt)

input <- './Results/Tian/infercnv/10x_5cl_201/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts.txt.gz')
sample_input <- paste0(input,'/sample_anno.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- paste0(input,'/results')
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=0.1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 40)



#############
## celseq2 ##
#############

library(infercnv)
library(biomaRt)

input <- './Results/Tian/infercnv/celseq2/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts.txt.gz')
sample_input <- paste0(input,'/sample_anno.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- paste0(input,'/results')
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 20)



#############
## dropseq ##
#############

library(infercnv)
library(biomaRt)

input <- './Results/Tian/infercnv/dropseq/'
input2 <- './Anno/'
raw_input <- paste0(input,'/raw_counts.txt.gz')
sample_input <- paste0(input,'/sample_anno.txt.gz')
gene_input <- paste0(input2,'gene_anno_inferCNV.txt')

output <- paste0(input,'/results')
dir.create(output,recursive=T)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=raw_input,annotations_file=sample_input,delim="\t",gene_order_file=gene_input, ref_group_names="Normal")
infercnv_obj <- infercnv::run(infercnv_obj, cutoff=1, out_dir=output, denoise=T, HMM=T, analysis_mode='subclusters', tumor_subcluster_partition_method = 'random_trees', tumor_subcluster_pval=0.05, num_threads = 20)





