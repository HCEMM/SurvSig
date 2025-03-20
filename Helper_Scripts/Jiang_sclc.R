
fpkm = read.table("./source_data/GSE60052_79tumor.7normal.normalized.log2.data.Rda.tsv", header = T, sep = "\t", row.names = 1, check.names = F)
colnames(fpkm) = gsub(pattern = "A", replacement = "", x = colnames(fpkm))
colnames(fpkm) = gsub(pattern = "B", replacement = "", x = colnames(fpkm))
colnames(fpkm) = gsub(pattern = "T", replacement = "", x = colnames(fpkm))

clinical = as.data.frame(readxl::read_xlsx("./source_data/pgen.1005895.s012.xlsx"))
colnames(clinical) = clinical[1,]
clinical = clinical[2:nrow(clinical),]
rownames(clinical) = clinical[,1]
clinical[,1] = NULL
clinical = clinical[1:(nrow(clinical)-1),]
clinical = clinical[!is.na(clinical$`Tumor RNAseq data`),]
clinical = clinical[order(rownames(clinical)),]
colnames(fpkm)[!colnames(fpkm) %in% rownames(clinical)]
clinical$Gender = ifelse(as.numeric(clinical$Gender) > 1, "Female", "Male")
clinical$Smoke = ifelse(as.numeric(clinical$Smoke) > 1, "No", "Yes")
clinical[,c(1,2,3,4,5,6,7,8,13)] = NULL

clinical = clinical[rownames(clinical) %in% colnames(fpkm),]
fpkm = fpkm[,rownames(clinical)]

clustering = clinical[,c("Months", "Survival Status")]
colnames(clustering) = c("Months", "status")
clustering$Months = as.numeric(clustering$Months)
clustering$status = as.numeric(clustering$status)
clustering = as.data.frame(clustering)

fpkm = as.matrix(fpkm)
clustering$cluster = ifelse(as.numeric(fpkm["SRSF1",]) > quantile(fpkm["SRSF1",],.75)[[1]], "high", "down")

fit <- survfit(Surv(Months ,status)~cluster, data = clustering)
ggsurvplot(fit, risk.table = T, pval = TRUE)

#################### EMT Score ####################

load(file='source_data/EMT_gene_list.rda')
gene_list = EMT_gene_list
rownames(gene_list) = gene_list$genes
gene_list$genes = NULL

filtered_log = fpkm[rownames(fpkm) %in% rownames(gene_list), ]

# Convert filtered data to a matrix
filtered_log <- as.matrix(filtered_log)

# Define a function to calculate EMT scores - https://github.com/korkutlab/imogimap/
get_emt_score <- function(pdf) {
  pdf <- as.matrix(pdf)
  missing_EMT <- EMT_gene_list[-which(EMT_gene_list$genes %in% rownames(pdf)),]$genes
  if(length(missing_EMT) > 0) {
    warning(length(missing_EMT)," missing EMT signature genes: ", 
            paste(missing_EMT, collapse = " "), "\nCheck EMT_gene_list for all signature genes.\n")
  }
  
  sub_list <- EMT_gene_list[EMT_gene_list$genes %in% rownames(pdf),]
  pdf_sub <- pdf[rownames(pdf) %in% sub_list$genes,]
  
  if(nrow(pdf_sub) == 0) {
    warning("No signature genes found!")
    return(tibble(Tumor_Sample_ID = character(), score = numeric()))
  } else {
    #pdf_sub <- log2(pdf_sub + 1)
    pdf_sub <- t(scale(t(pdf_sub), center = TRUE, scale = TRUE))
    pdf_sub <- sweep(pdf_sub, 1, sub_list$sign, "*")
    pscore <- colMeans(pdf_sub, na.rm = TRUE)
    pscore <- data.frame(Tumor_Sample_ID = colnames(pdf_sub), score = pscore)
  }
  return(pscore)
}

# Calculate EMT score using the function
emt_score_result <- get_emt_score(filtered_log)

# Extract the 'score' column
emt_scores <- emt_score_result[, "score", drop = FALSE]

# Define common_samples if not already defined (assuming it is a vector of sample names)
common_samples <- intersect(colnames(filtered_log), rownames(emt_scores))

# Align the EMT scores with the common samples
emt_scores_aligned <- emt_scores[common_samples, , drop = FALSE]
clinical$EMT = ifelse(as.numeric(emt_scores_aligned$score) > 0,  "Mesenchymal", "Epithelial")

#################### NE Score ####################

# Function to calculate NE Score
ne_Score <- function(input_data, ne50_data) {
  nb <- ncol(input_data)
  cat("Starting number of samples:", nb, "\n")
  
  # Find common genes between input and NE50 data
  cnames <- intersect(rownames(input_data), rownames(ne50_data))
  if (length(cnames) < 50) {
    cat("Number of common genes:", length(cnames), "\n")
    cat("Missing genes:", setdiff(rownames(ne50_data), cnames), "\n")
  }
  ne50_data <- ne50_data[cnames, ]
  input_data <- input_data[cnames, ]
  stopifnot(identical(rownames(input_data), rownames(ne50_data)))
  
  # Remove columns (cell lines) with all NA values
  nba <- colSums(is.na(input_data))
  if (any(nba == nrow(input_data))) {
    input_data <- input_data[, nba != nrow(input_data)]
    nb <- ncol(input_data)
  }
  cat("Final number of samples (no NA):", nb, "\n")
  
  # Initialize result matrix
  res <- matrix(NA, nb, 3)
  colnames(res) <- c("cor.NE", "cor.Non.NE", "score")
  rownames(res) <- colnames(input_data)
  
  # Calculate correlations and scores
  for (k in 1:nb) {
    res[k, 1] <- cor(ne50_data$NE.mean, input_data[, k], use = "pairwise.complete.obs")
    res[k, 2] <- cor(ne50_data$Non.NE.mean, input_data[, k], use = "pairwise.complete.obs")
    res[k, 3] <- (res[k, 1] - res[k, 2]) / 2
  }
  
  return(res)
}

# Load NE50 signatures
gene_list <- read.table("source_data/NE50-signatures.txt", header = TRUE, sep = "\t", comment.char = "")
gene_names <- gene_list$Symbol

ne50signature_file <- "source_data/NE50-signatures.txt"
ne50sign <- read.delim(ne50signature_file, stringsAsFactors = FALSE, row.names = 1)
filtered_data <- fpkm[rownames(fpkm) %in% gene_names, ]

# Calculate NE scores
ne_score_result <- ne_Score(filtered_data, ne50sign)

# Extract 'score' column
ne_scores <- ne_score_result[, "score", drop = FALSE]

clinical$NE = ifelse(as.numeric(emt_scores_aligned$score) > 0,  "NE", "non-NE")

#################### NAPY Classification ####################

napy_genes = c("ASCL1", "NEUROD1", "POU2F3", "YAP1")
ext_genes = c(napy_genes, c("CHGA", "SYP"), c("REST", "NOTCH1", "NOTCH2", "NOTCH3"))
gene_split = c("LTF", "LTF", "LTF", "LTF", "NE marker", "NE marker", "NOTCH", "NOTCH", "NOTCH", "NOTCH")
gene_split = factor(gene_split, levels = c("NE marker", "LTF", "NOTCH"))
clust_array = c()
tmp_exp_df = fpkm[napy_genes,]

for(i in 1:ncol(tmp_exp_df)) {
  clust_array[i] = rownames(tmp_exp_df)[as.numeric(tmp_exp_df[,i]) == max(as.numeric(tmp_exp_df[,i]))]
}

clinical$NAPY = clust_array
