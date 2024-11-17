
## version 1 use single cell data as reference ##

#################
## 10x ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/10x/gene_counts/read_len_subsampling_SE/100bp/'
output <- './Results/10x/all_cell/casper/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/10x/gene_counts/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)

colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Data/10x/gene_counts/HCC1395BL/gene_counts_featureCounts.rdata')
gene_hcc1395bl <- gene_counts
colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395_2.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
maf <- read.table(paste0(input,'/down_sample_merge_hcc1395bl_2.baf'),sep='\t')
loh[[2]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
names(loh) <- c('HCC1395','HCC1395BL')
loh.name.mapping <- rbind.data.frame(cbind('HCC1395',colnames(gene_hcc1395)),cbind('HCC1395BL',colnames(gene_hcc1395bl)))
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=colnames(gene_hcc1395bl),cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))



#################
## c1_fda ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/C1_ceber/gene_counts/read_len_subsampling_SE/125bp/500K/'
output <- './Results/c1_fda/all_cell/casper/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/C1_ceber/gene_counts/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)

colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Data/C1_ceber/gene_counts/HCC1395BL/gene_counts_featureCounts.rdata')
gene_hcc1395bl <- gene_counts
colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
maf <- read.table(paste0(input,'/down_sample_merge_hcc1395bl.baf'),sep='\t')
loh[[2]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
names(loh) <- c('HCC1395','HCC1395BL')
loh.name.mapping <- rbind.data.frame(cbind('HCC1395',colnames(gene_hcc1395)),cbind('HCC1395BL',colnames(gene_hcc1395bl)))
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=colnames(gene_hcc1395bl),cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))




#################
## C1 ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/C1/gene_counts/read_len_subsampling_SE/150bp/1M/'
output <- './Results/c1_llu/all_cell/casper/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/C1/gene_counts/read_len_SE/150bp/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)

colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Data/C1/gene_counts/read_len_SE/150bp/HCC1395BL/gene_counts_featureCounts.rdata')
gene_hcc1395bl <- gene_counts
colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
maf <- read.table(paste0(input,'/down_sample_merge_hcc1395bl.baf'),sep='\t')
loh[[2]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
names(loh) <- c('HCC1395','HCC1395BL')
loh.name.mapping <- rbind.data.frame(cbind('HCC1395',colnames(gene_hcc1395)),cbind('HCC1395BL',colnames(gene_hcc1395bl)))
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
annotation=annotation2,method='iterative',loh=loh,control.sample.ids=colnames(gene_hcc1395bl),cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))

gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))


#################
## icell8 ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/WaferGen/gene_counts/read_len_subsampling_SE/150bp/1M/'
output <- './Results/WaferGen/all_cell/casper/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/WaferGen/gene_counts/read_len_SE/150bp/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)

colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')
load('./Data/WaferGen/gene_counts/read_len_SE/150bp/HCC1395BL/gene_counts_featureCounts.rdata')
gene_hcc1395bl <- gene_counts
colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
maf <- read.table(paste0(input,'/down_sample_merge_hcc1395bl.baf'),sep='\t')
loh[[2]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
names(loh) <- c('HCC1395','HCC1395BL')
loh.name.mapping <- rbind.data.frame(cbind('HCC1395',colnames(gene_hcc1395)),cbind('HCC1395BL',colnames(gene_hcc1395bl)))
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
annotation=annotation2,method='iterative',loh=loh,control.sample.ids=colnames(gene_hcc1395bl),cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))

gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
rna.matrix <- rna.matrix[,1:598]

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))





## version 2 use bulk data as reference ##

#################
## 10x ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/10x/gene_counts/read_len_subsampling_SE/100bp/'
input2 <- './Data/RNAseq/merge/Star/NCBI_GRCh38/Gtf/'
output <- './Results/10x/all_cell/casper/v2/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/10x/gene_counts/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input2,'txi_star_rsem_hcc1395bl.rdata'))
gene_hcc1395bl <- txi_rsem_hcc1395bl$counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]
gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395_2.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))




#################
## c1-fda ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/C1_ceber/gene_counts/read_len_subsampling_SE/125bp/500K/'
input2 <- './Data/RNAseq/merge/Star/NCBI_GRCh38/Gtf/'
output <- './Results/c1_fda/all_cell/casper/v2/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/C1_ceber/gene_counts/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input2,'txi_star_rsem_hcc1395bl.rdata'))
gene_hcc1395bl <- txi_rsem_hcc1395bl$counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))




#################
## C1 ##
#################


library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/C1/gene_counts/read_len_subsampling_SE/150bp/1M/'
input2 <- './Data/RNAseq/merge/Star/NCBI_GRCh38/Gtf/'
output <- './Results/c1_llu/all_cell/casper/v2/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/C1/gene_counts/read_len_SE/150bp/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input2,'txi_star_rsem_hcc1395bl.rdata'))
gene_hcc1395bl <- txi_rsem_hcc1395bl$counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))




#################
## icell8 ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/WaferGen/gene_counts/read_len_subsampling_SE/150bp/1M/'
input2 <- './Data/RNAseq/merge/Star/NCBI_GRCh38/Gtf/'
output <- './Results/WaferGen/all_cell/casper/v2/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/WaferGen/gene_counts/read_len_SE/150bp/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input2,'txi_star_rsem_hcc1395bl.rdata'))
gene_hcc1395bl <- txi_rsem_hcc1395bl$counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))



#################
## bulk rnaseq ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/RNAseq/merge/Star/NCBI_GRCh38/Gtf/'
output <- './Results/bulk/all_cell/casper/v2/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load(paste0(input,'txi_star_rsem_hcc1395.rdata'))
gene_hcc1395 <- txi_rsem_hcc1395$counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input,'txi_star_rsem_hcc1395bl.rdata'))
gene_hcc1395bl <- txi_rsem_hcc1395bl$counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/HCC1395/merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
maf <- read.table(paste0(input,'/HCC1395BL/merge_hcc1395bl.baf'),sep='\t')
loh[[2]] <- data.frame(chr = maf[, 1], position = maf[, 2], 
            alt = maf[, 5], ref = maf[, 6] - maf[, 5], coverage = maf[, 
                6], baf = maf[, 5]/maf[, 6], dev = abs(as.numeric(maf[, 
                5]/maf[, 6]) - 0.5))
names(loh) <- c('HCC1395','HCC1395BL')
loh.name.mapping <- rbind.data.frame(cbind('HCC1395',colnames(gene_hcc1395)),cbind('HCC1395BL',colnames(gene_hcc1395bl)))
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='bulk',cnv.scale=3,loh.scale=3,
annotation=annotation2,method='iterative',loh=loh,control.sample.ids=colnames(gene_hcc1395bl),cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))



## version 3 use GTex bulk data as reference ##

#################
## 10x ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/10x/gene_counts/read_len_subsampling_SE/100bp/'
input2 <- './Anno/GTEx/'
output <- './Results/10x/all_cell/casper/v3/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/10x/gene_counts/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input2,'breast_gene_counts.Rdata'))
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]
gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395_2.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))




#################
## c1-fda ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/C1_ceber/gene_counts/read_len_subsampling_SE/125bp/500K/'
input2 <- './Anno/GTEx/'
output <- './Results/c1_fda/all_cell/casper/v3/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/C1_ceber/gene_counts/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input2,'breast_gene_counts.Rdata'))
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))




#################
## C1 ##
#################


library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/C1/gene_counts/read_len_subsampling_SE/150bp/1M/'
input2 <- './Anno/GTEx/'
output <- './Results/c1_llu/all_cell/casper/v3/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/C1/gene_counts/read_len_SE/150bp/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input2,'breast_gene_counts.Rdata'))
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]
gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[grep('_1',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))




#################
## icell8 ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Data/WaferGen/gene_counts/read_len_subsampling_SE/150bp/1M/'
input2 <- './Anno/GTEx/'
output <- './Results/WaferGen/all_cell/casper/v3/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/WaferGen/gene_counts/read_len_SE/150bp/HCC1395/gene_counts_featureCounts.rdata')
gene_hcc1395 <- gene_counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input2,'breast_gene_counts.Rdata'))
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[-grep('_2',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))



#################
## bulk rnaseq ##
#################

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Anno/GTEx/'
output <- './Results/bulk/all_cell/casper/v3/'
dir.create(output,recursive=T)

cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
load('./Data/RNAseq/merge/Star/NCBI_GRCh38/Gtf/txi_star_rsem_hcc1395.rdata')
gene_hcc1395 <- txi_rsem_hcc1395$counts
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_hcc1395), ishg19=F, centromere)
colnames(gene_hcc1395) <- paste0(colnames(gene_hcc1395),'_1')

load(paste0(input,'breast_gene_counts.Rdata'))
gene_hcc1395bl <- gene_counts
control.sample.ids <- colnames(gene_hcc1395bl) <- paste0(colnames(gene_hcc1395bl),'_2')
overlap_gene <- intersect(rownames(gene_hcc1395),rownames(gene_hcc1395bl))
index <- match(annotation$Gene,overlap_gene)
annotation2 <- annotation[!is.na(index),]

gene_exprs <- t(cbind(gene_hcc1395[annotation2$Gene,],gene_hcc1395bl[annotation2$Gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
gene_exprs <- as.matrix(gene_exprs)

maf <- read.table(paste0(input,'/down_sample_merge_hcc1395.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('HCC1395')
loh.name.mapping <- cbind.data.frame('HCC1395',colnames(gene_hcc1395))
levels(loh.name.mapping[,1])=c('HCC1395',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='bulk',cnv.scale=3,loh.scale=3,
annotation=annotation2,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))
l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))


gamma <- 7
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loss.final <- loss[loss$count>=gamma, ]
gain.final <- gain[gain$count>=gamma, ]
all.summary<- rbind(loss.final, gain.final)
all.summary <- all.summary[-grep('_2',all.summary$ID),]
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "Gene")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,1])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
all.samples <- all.samples[grep('_1',all.samples)]
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

save(rna.matrix,file=paste0(output,'rna.matrix.rdata'))






