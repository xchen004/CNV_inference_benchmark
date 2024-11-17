
#############
## casper ##
#############

library(CaSpER)
library(biomaRt)
library(GenomicFeatures)

input <- './Results/SCLC/casper/baf/'
output <- './Results/SCLC/casper/results/'
dir.create(output,recursive=T)

load('./Results/SCLC/Star_featureCounts/Ensembl_GRCh38/Gtf/gene_counts.rdata')
gene_sclc <- gene_counts
colnames(gene_sclc) <- paste0(colnames(gene_sclc),'_1')
cytoband <- read.table('./Anno/cyto/cytoBand_hg38.txt')
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(gene_sclc), ishg19=F, centromere)

load('./Anno/GTEx/lung_gene_counts.Rdata')
gene_ctrl <- gene_counts
control.sample.ids <- colnames(gene_ctrl) <- paste0(colnames(gene_ctrl),'_2')
overlap_gene <- intersect(rownames(gene_sclc),rownames(gene_ctrl))
gene_exprs <- t(cbind(gene_sclc[overlap_gene,],gene_ctrl[overlap_gene,]))
gene_exprs <- t(log(gene_exprs/apply(gene_exprs,1,sum)*1e6+1))
index <- match(annotation$Gene,rownames(gene_exprs))
gene_exprs <- gene_exprs[index[!is.na(index)],]
annotation <- annotation[!is.na(index),]

maf <- read.table(paste0(input,'sample_merge.baf'),sep='\t')
loh <- list()
loh[[1]] <- data.frame(chr = maf[,1], position = maf[,2], alt = maf[,5], ref = maf[,6] - maf[,5], 
		coverage = maf[,6], baf = maf[,5]/maf[,6], dev = abs(as.numeric(maf[,5]/maf[,6]) - 0.5))
names(loh) <- c('sclc')
loh.name.mapping <- cbind.data.frame('sclc',colnames(gene_sclc))
levels(loh.name.mapping[,1])=c('sclc',control.sample.ids)
colnames(loh.name.mapping) <- c('loh.name','sample.name')

object <- CreateCasperObject(raw.data=gene_exprs,loh.name.mapping=loh.name.mapping,sequencing.type='single-cell',cnv.scale=3,loh.scale=3,
expr.cutoff=0.1,annotation=annotation,method='iterative',loh=loh,control.sample.ids=control.sample.ids,cytoband=cytoband)
final.objects <- runCaSpER(object, removeCentromere=F, cytoband=cytoband, method="iterative")
save(final.objects,file=paste0(output,'/casper_final.rdata'))

all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects)
save(segment.summary,file=paste0(output,'/casper_segment.rdata'))

l_cnv <- extractLargeScaleEvents(final.objects, thr=0.75)
save(l_cnv,file=paste0(output,'/large_cnv.rdata'))



