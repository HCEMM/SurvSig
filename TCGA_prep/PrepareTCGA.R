library(TCGAbiolinks)

# Data is downloaded from the pan-cancer atlas:
# https://gdc.cancer.gov/about-data/publications/pancanatlas

# Manifest file: https://gdc.cancer.gov/files/public/file/PanCan-General_Open_GDC-Manifest_3.txt
#id	filename	md5	size
#7d4c0344-f018-4ab0-949a-09815f483480	merge_merged_reals.tar.gz	ff8bf50dcd3a314162af71d1b8e538b6	388646603
#0f4f5701-7b61-41ae-bda9-2805d1ca9781	TCGA_mastercalls.abs_segtabs.fixed.txt	585c8793730f294d7bf0144566bb37fa	253061161
#1a7d7be8-675d-4e60-a105-19d4121bdebf	merged_sample_quality_annotations.tsv	05ddd2270fb1fb24fbdc2fe9bf7384e5	8463670
#55d9bf6f-0712-4315-b588-e6f8e295018e	PanCanAtlas_miRNA_sample_information_list.txt	02bb56712be34bcd58c50d90387aebde	553408
#d82e2c44-89eb-43d9-b6d3-712732bf6a53	jhu-usc.edu_PANCAN_merged_HumanMethylation27_HumanMethylation450.betaValue_whitelisted.tsv	5cec086f0b002d17befef76a3241e73b	5022150019
#4f277128-f793-4354-a13d-30cc7fe9f6b5	TCGA_mastercalls.abs_tables_JSedit.fixed.txt	8ea2ca92c8ae58350538999dfa1174da	901812
#00a32f7a-c85f-4f86-850d-be53973cbc4d	broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg	7fafc537807d5b3ddf0bb89665279a9d	166238394
#0fc78496-818b-4896-bd83-52db1f533c5c	clinical_PANCAN_patient_with_followup.tsv	ffcb35edda305dd8d615497f9214eb92	18633685
#1c8cfe5f-e52d-41ba-94da-f15ea1337efc	mc3.v0.2.8.PUBLIC.maf.gz	639ad8f8386e98dacc22e439188aa8fa	753339089
#1b5f413e-a8d1-4d10-92eb-7c4ae739ed81	TCGA-CDR-SupplementalTableS1.xlsx	a4591b2dcee39591f59e5e25a6ce75fa	2945129
#fcbb373e-28d4-4818-92f3-601ede3da5e1	TCGA-RPPA-pancan-clean.txt	e2b914c7ecd369589275d546d9555b05	18901234
#3586c0da-64d0-4b74-a449-5ff4d9136611	EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv	02e72c33071307ff6570621480d3c90b	1882540959
#1c6174d9-8ffb-466e-b5ee-07b204c15cf8	pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.csv	c7501dc3c505ca172a6a05b611bd11c3	67167640
#99b0c493-9e94-4d99-af9f-151e46bab989	jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv	a92f50490cf4eca98b0d19e10927de9d	41541692788
#210c0e4d-de0b-4e00-b323-79991cac8f66	PANCAN_ISAR.seg.txt	ee8d16b55290754d6e7ac76d25989b47	104660421

expression_manifest_id = "3586c0da-64d0-4b74-a449-5ff4d9136611"
desination_file = "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
system(paste0("curl --remote-name https://api.gdc.cancer.gov/data/", 
              expression_manifest_id))
file.rename(expression_manifest_id, desination_file)


# Clinical data: downloaded, as suggested on GDC site, from the TCGA-CDR paper
# https://www.cell.com/cell/fulltext/S0092-8674(18)30229-0
# https://www.cell.com/cms/10.1016/j.cell.2018.02.052/attachment/bbf46a06-1fb0-417a-a259-fd47591180e4/mmc1.xlsx
clinical_manifest_id = "1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"
clinical_xls_file = "TCGA-CDR-SupplementalTableS1.xlsx"

system(paste0("curl --remote-name https://api.gdc.cancer.gov/data/", 
              clinical_manifest_id))
file.rename(clinical_manifest_id, clinical_xls_file)
# import manually downloaded excel file
clinical_xls = readxl::read_xlsx(clinical_xls_file)
# remove redacted samples
clinical_df = as.data.frame(clinical_xls[is.na(clinical_xls$Redaction),])
rownames(clinical_df) = clinical_df$bcr_patient_barcode
clinical_df[,1] = NULL


# import expression data
TCGA_expression = read.table(desination_file, header = T, sep = "\t", 
                             check.names = F, comment.char = "")
rownames(TCGA_expression) = TCGA_expression$gene_id
tcga_genes = data.frame(original_name = TCGA_expression$gene_id,
                        SYMBOL = NA,
                        ENTREZ = NA)
tcga_genes$SYMBOL = sapply(1:nrow(tcga_genes), function(i) {
  strsplit(x = tcga_genes$original_name[i], split = "\\|")[[1]][1]
})
tcga_genes$ENTREZ = sapply(1:nrow(tcga_genes), function(i) {
  strsplit(x = tcga_genes$original_name[i], split = "\\|")[[1]][2]
})
tcga_genes$mean = rowMeans(TCGA_expression)
TCGA_expression$gene_id = NULL

tcga_genes = tcga_genes[tcga_genes$SYMBOL != "?",]
tcga_genes = tcga_genes[order(tcga_genes$mean, decreasing = T),]
tcga_genes = tcga_genes[!duplicated(tcga_genes$SYMBOL),]
TCGA_expression = TCGA_expression[tcga_genes$original_name,]

sample_df = data.frame(name = colnames(TCGA_expression),
                       sample_type = substr(colnames(TCGA_expression), 14, 15),
                       short_name = substr(colnames(TCGA_expression), 1, 12))

sample_df$type = ifelse(as.numeric(sample_df$sample_type) > 10, "normal", "tumor")
sample_df = sample_df[order(as.numeric(sample_df$sample_type), decreasing = F),]
sample_df = sample_df[sample_df$type == "tumor",]
sample_df = sample_df[!duplicated(sample_df$short_name),]

TCGA_expression = TCGA_expression[,sample_df$name]
colnames(TCGA_expression) = sample_df$short_name

clinical_df = clinical_df[clinical_df$bcr_patient_barcode %in% colnames(TCGA_expression),]
TCGA_expression = TCGA_expression[,clinical_df$bcr_patient_barcode]

write.table(TCGA_expression, 
            "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.renamed.tsv",
            col.names = NA, row.names = T, sep = "\t", quote = F)

write.table(clinical_df, 
            "TCGA-CDR_clinical_data.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

# Prepare TCGA subtype information
# http://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/subtypes.html
TCGA_molecular_subtypes <- as.data.frame(PanCancerAtlas_subtypes())
rownames(TCGA_molecular_subtypes) = TCGA_molecular_subtypes$pan.samplesID

sample_df = data.frame(name = rownames(TCGA_molecular_subtypes),
                       sample_type = 0,
                       short_name = substr(rownames(TCGA_molecular_subtypes), 1, 12))

sample_df$sample_type = sapply(1:nrow(sample_df), function(i) {
  if(nchar(sample_df$name[i]) > 12) {
    as.numeric(substr(sample_df$name[i], 14, 15))
  } else {
    0
  }
})

sample_df = sample_df[order(sample_df$sample_type, decreasing = F),]
sample_df = sample_df[!duplicated(sample_df$short_name),]
TCGA_molecular_subtypes = TCGA_molecular_subtypes[sample_df$name,]
rownames(TCGA_molecular_subtypes) = sample_df$short_name

write.table(TCGA_molecular_subtypes, 
            "TCGAbiolinks_subtypes.tsv",
            col.names = T, row.names = F, sep = "\t", quote = F)

