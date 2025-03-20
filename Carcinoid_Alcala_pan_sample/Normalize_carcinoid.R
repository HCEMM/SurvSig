library(rtracklayer)
library(edgeR)
library(DESeq2)
library(ComplexHeatmap)
library(umap)
library(ggplot2)

# paper sources: 
#       https://academic.oup.com/gigascience/article/9/11/giaa112/5943495?login=false#supplementary-data
#       https://www.nature.com/articles/s41467-019-11276-9#MOESM1

#data from https://github.com/IARCbioinfo/DRMetrics/tree/NextJournalH

#gtf: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
gtf = rtracklayer::import("gencode.v33.annotation.gtf.gz")
gtf.exon = gtf[gtf$type == "exon"]
grl <- reduce(split(gtf.exon, elementMetadata(gtf.exon)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
reducedGTF$width = width(reducedGTF)
reducedGTF.df = reducedGTF
names(reducedGTF.df) = NULL

reducedGTF.df = as.data.frame(reducedGTF.df)
reducedGTF.df$ENS = names(reducedGTF)
aggreg_lengths = aggregate(reducedGTF.df$width, by=list(Category=reducedGTF.df$ENS), FUN=sum)
aggreg_lengths = as.data.frame(aggreg_lengths)
colnames(aggreg_lengths) = c("ENSEMBL", "width")
rownames(aggreg_lengths) = aggreg_lengths$ENSEMBL
gene_maps = data.frame(
  ENS = gtf.exon$gene_id,
  SYMBOL = gtf.exon$gene_name
)

gene_maps = gene_maps[complete.cases(gene_maps),]
gene_maps = gene_maps[!duplicated(gene_maps$ENS),]
rownames(gene_maps) = gene_maps$ENS
aggreg_lengths$SYMBOL = gene_maps[aggreg_lengths$ENSEMBL, "SYMBOL"]

raw_counts = read.table("read_counts_all.txt", header = T, sep = " ", row.names = 1, check.names = F)
attrib = read.table("Attributes.txt", header = T, check.names = F, sep = "\t")

carc_types = c("Atypical", "Typical")
attrib.carc = attrib[attrib$Histopathology %in% carc_types,]

raw_counts.carc = raw_counts[,attrib.carc$Sample_ID]

aggreg_lengths = aggreg_lengths[rownames(aggreg_lengths) %in% rownames(raw_counts.carc),]
raw_counts.carc = raw_counts.carc[rownames(aggreg_lengths),]

gene_lengths = data.frame(
  GeneID = rownames(raw_counts.carc),
  Length = aggreg_lengths[rownames(raw_counts.carc), "width"]
)

GeneDF_EdgeR <- edgeR::DGEList(counts = raw_counts.carc, genes = as.data.frame(gene_lengths))
GeneDF_Norm  <- edgeR::calcNormFactors(GeneDF_EdgeR, method = 'TMM')
tmm_fpkm <- edgeR::rpkm(GeneDF_Norm, normalized.lib.sizes = TRUE, log = F)
tmm_fpkm_log = log2(tmm_fpkm + 1)

aggreg_lengths$sum = rowSums(tmm_fpkm_log)
aggreg_lengths = aggreg_lengths[order(aggreg_lengths$sum, decreasing = T),]
aggreg_lengths.uniq = aggreg_lengths[!duplicated(aggreg_lengths$SYMBOL),]

tmm_fpkm_log = tmm_fpkm_log[aggreg_lengths.uniq$ENSEMBL,]
rownames(tmm_fpkm_log) = aggreg_lengths.uniq$SYMBOL

genes = c("ASCL1", "NEUROD1", "POU2F3", "YAP1")
Heatmap(tmm_fpkm_log[genes,])

sds = data.frame(gene = rownames(tmm_fpkm_log),
                 sd = apply(tmm_fpkm_log, 1, sd))
sds = sds[order(sds$sd, decreasing = T),]

Heatmap(tmm_fpkm_log[sds$gene[1:1000],], show_row_names = F, name = "Log2")

tmm_fpkm_log.z = as.data.frame(t(scale(t(tmm_fpkm_log))))
Heatmap(tmm_fpkm_log.z[sds$gene[1:1000],], show_row_names = F, name = "Zscored")

umap.config = umap.defaults
umap.config$n_neighbors = 15

umap_res = umap(t(tmm_fpkm_log[sds$gene[1:5000],]), config = umap.config)
pldf = as.data.frame(umap_res$layout)
colnames(pldf) = c("UMAP_1", "UMAP_2")
pldf$source = "LNEN"
pldf$source = ifelse(grepl(pattern = "SRR", x = rownames(pldf)), "SRA", pldf$source)
pldf$source = ifelse(grepl(pattern = "S0", x = rownames(pldf)), "S0", pldf$source)

ggplot(pldf,aes(x = UMAP_1, y = UMAP_2, fill = source)) +
  geom_point(colour="black",pch=21, size=5) +
  theme_bw()

pcares = prcomp(t(tmm_fpkm_log[sds$gene[1:5000],]),
                center = TRUE,
                scale. = F)
pldf = as.data.frame(pcares$x[,1:2])
pldf$source = "LNEN"
pldf$source = ifelse(grepl(pattern = "SRR", x = rownames(pldf)), "SRA", pldf$source)
pldf$source = ifelse(grepl(pattern = "S0", x = rownames(pldf)), "S0", pldf$source)

ggplot(pldf,aes(x = PC1, y = PC2, fill = source)) +
  geom_point(colour="black",pch=21, size=5) +
  theme_bw()

#write.table(tmm_fpkm_log, "Carcinoid_pan_sample.tsv", col.names = NA, row.names = T, sep = "\t", quote = F)

sample_ann = attrib.carc
rownames(sample_ann) = sample_ann$Sample_ID
sample_ann = sample_ann[,c("Sample_ID", "Histopathology",
                           "Age", "Sex", "Smoking_status",
                           "Histopathology_simplified", "Molecular_clusters",
                           "Histopathology_MachineLearning_prediction")]
#write.table(sample_ann, "Sample_annotation.tsv", col.names = T, row.names = F,
#            sep = "\t", quote = F)
