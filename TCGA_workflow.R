
library(TCGAbiolinks)
library(DT)
library(openxlsx)
library(ELMER)
library(MultiAssayExperiment)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(readr)
library(conflicted)




# This workflow is designed to perform Transcriptome and methylation analysis of TCGA data
# Here any one of the four BRCA subtypes can be selected for differential expression, miRNA-gene correlation, and methylation analysis

############################################# Section 1: Download data ########################################################


################## Subtypes notes ####################

# Luminal A tumors are characterized by the presence of ER and/or PR and the absence of HER2,
# and have a low expression of cell proliferation marker Ki-67 (less than 20%) (Figure 1).
# Clinically they are low grade, slow growing, and have the best prognosis with less incidence
# of relapse and higher survival rate.

# Luminal B tumors are of higher grade and worse prognosis compared to Luminal A.
# They are ER positive and can be PR negative and have a high expression of Ki67
# (greater than 20%) (Figure 2). They are generally of intermediate/high histologic grade.
# These tumors may benefit from hormonal therapy along with chemotherapy. The elevated Ki67 makes
# them grow faster than luminal A and worse prognosis (32).

# The HER2-positive group constitutes 10–15% of breast cancers and is characterized by high HER2
# expression with absence of ER and PR (Figure 3). They grow faster than the luminal ones and the
# prognosis has improved after the introduction of HER2-targeted therapies. The HER2-positive subtype
# is more aggressive and fast-growing. Within this, two subgroups can be distinguished: luminal HER2
# (E+, PR+, HER2+ and Ki-67:15–30%) and HER2-enriched (HER2+, E-, PR-, Ki-67>30%) (36). They have a worse
# prognosis compared to luminal tumors

# Triple-negative breast cancer is ER-negative, PR-negative, and HER2-negative (Figure 4). They constitute
# about 20% of all breast cancers. It is most common among women under 40 years of age, and in African-American women.
# The TNBC subtype is further classified into several additional subgroups including basal-like (BL1 and BL2),
# claudin-low, mesenchymal (MES), luminal androgen receptor (LAR), and immunomodulatory (IM), the first two being
# the most frequent with 50–70% and 20–30% of cases (41). Moreover, each of these has unique clinical outcomes,
# phenotypes, and pharmacological sensitivities. TNBC presents an aggressive behavior and 80% of breast cancer tumors
# (tumor suppressor gene BRCA1 and BRCA2) belong to this group (28).


cancer_type = "TCGA-BRCA"

################## Section 2: Gene Preprocessing and Differential Expression ##################################################################


# Download tumor sample Gene expression
query_exp_tumor <- GDCquery(
  project = cancer_type,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor")
)

GDCdownload(query_exp_tumor)
exp_tumor <- GDCprepare(
  query = query_exp_tumor
)

# Download normal sample Gene expression
query_exp_normal <- GDCquery(
  project = cancer_type,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Solid Tissue Normal")
)

GDCdownload(query_exp_normal)
exp_normal <- GDCprepare(
  query = query_exp_normal
)

################ BRCA subtype selection #####################################
clinical_data <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

metadataObj <- colData(exp_tumor) # contains dataframe with PAM50 tumor classifications and other clinical data for each sample
#luminal_A_samples <- metadata[!is.na(metadata$paper_BRCA_Subtype_PAM50) & metadata$paper_BRCA_Subtype_PAM50 == "LumA",]
luminal_B_samples <- metadata[!is.na(metadata$paper_BRCA_Subtype_PAM50) & metadata$paper_BRCA_Subtype_PAM50 == "LumB",]
#HER2_samples <- metadata[!is.na(metadata$paper_BRCA_Subtype_PAM50) & metadata$paper_BRCA_Subtype_PAM50 == "Her2",]
#Basal_samples <- metadata[!is.na(metadata$paper_BRCA_Subtype_PAM50) & metadata$paper_BRCA_Subtype_PAM50 == "Basal",]

# It is important to note that sample ids in exp_tumor is large ids with extra characters and samples$bcr_patient_barcode
# has tructued id format. Therefore to match these and pull subtype sample expression data, we need following
matched_columns <- sapply(colnames(exp_tumor), function(column_name) {
  any(sapply(luminal_B_samples$bcr_patient_barcode, function(truncated_id) {
    grepl(truncated_id, column_name)
  }))
})

exp_tumor_filtered <- exp_tumor[, matched_columns]
# Data preprocessing
exp_preprocessed_tumor <- TCGAanalyze_Preprocessing(
  object = exp_tumor_filtered,
  cor.cut = 0.9,
  datatype = "unstranded",
  filename = "IlluminaHiSeq_RNASeqV2.png"
)

exp_preprocessed_normal <- TCGAanalyze_Preprocessing(
  object = exp_normal,
  cor.cut = 0.9,
  datatype = "unstranded",
  filename = "IlluminaHiSeq_RNASeqV2_normal.png"
)

exp_normalized <- TCGAanalyze_Normalization(
  tabDF = cbind(exp_preprocessed_tumor, exp_preprocessed_normal),
  geneInfo = TCGAbiolinks::geneInfoHT,
  method = "gcContent"
)

exp_filtered <- TCGAanalyze_Filtering(
  tabDF = exp_normalized,
  method = "quantile",
  qnt.cut =  0.25
)

exp_filtered_tumor <- exp_filtered[
  ,colnames(exp_filtered) %in% colnames(exp_tumor_filtered)
]

exp_filtered_normal <-   exp_filtered[
  ,colnames(exp_filtered) %in% colnames(exp_normal)
]

# differential expression
diff_expressed_genes <- TCGAanalyze_DEA(
  mat1 = exp_filtered_normal,
  mat2 = exp_filtered_tumor,
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.05 ,
  logFC.cut = 0,
  method = "glmLRT"
)

write.xlsx(diff_expressed_genes, file = "diff_expressed_genes_HER2.xlsx")


################## Section 3: miRNA Preprocessing and Differential Expression ##################################################################

# Download mirna expression
query_miRNA_tumor <- GDCquery(
  project = cancer_type,
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type = c("Primary Tumor")
)

GDCdownload(query_miRNA_tumor)
miRNA_tumor <- GDCprepare(
  query = query_miRNA_tumor
)

query_miRNA_normal <- GDCquery(
  project = cancer_type,
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type = c("Solid Tissue Normal")
)

GDCdownload(query_miRNA_normal)
miRNA_normal <- GDCprepare(
  query = query_miRNA_normal
)

################ BRCA subtype selection (miRNA data) #####################################

# steps to find matching TCGA patient sample ids between gene and mirna data
extract_base_id_mirna <- function(full_id) {
  # Assuming the format is always "read_count_" followed by the TCGA ID
  id_parts <- strsplit(full_id, "_")[[1]]
  base_id <- id_parts[length(id_parts)] # Get the last part after splitting
  base_id <- unlist(strsplit(base_id, "-"))[1:3]
  return(paste(base_id, collapse = "-"))
}

extract_base_id <- function(full_id) {
  # Split the ID and extract the base part
  base_id <- unlist(strsplit(full_id, "-"))[1:3]
  return(paste(base_id, collapse = "-"))
}

contains_cross_mapped <- function(column_name) {
  return(grepl("cross-mapped", column_name))
}

contains_reads_per_million <- function(column_name) {
  return(grepl("reads_per_million", column_name))
}

rownames(miRNA_tumor) <- miRNA_tumor$miRNA_ID
rownames(miRNA_normal) <- miRNA_normal$miRNA_ID

# Apply to miRNA data
mirna_base_ids_tumor <- sapply(colnames(miRNA_tumor), extract_base_id_mirna)
mirna_base_ids_normal <- sapply(colnames(miRNA_normal), extract_base_id_mirna)
gene_expr_base_ids_tumor <- sapply(colnames(exp_filtered_tumor), extract_base_id)
gene_expr_base_ids_normal <- sapply(colnames(exp_filtered_normal), extract_base_id)

#Filter
miRNA_filtered_tumor <- miRNA_tumor[, !sapply(colnames(miRNA_tumor), contains_cross_mapped)]
miRNA_filtered_normal <- miRNA_normal[, !sapply(colnames(miRNA_normal), contains_cross_mapped)]
miRNA_filtered_tumor <- miRNA_tumor[, !sapply(colnames(miRNA_tumor), contains_reads_per_million)]
miRNA_filtered_normal <- miRNA_normal[, !sapply(colnames(miRNA_normal), contains_reads_per_million)]

# Function to remove "read_count_" prefix from the IDs
remove_read_count_prefix <- function(id) {
  gsub("read_count_", "", id)
}

colnames(miRNA_filtered_tumor) <- sapply(colnames(miRNA_filtered_tumor), remove_read_count_prefix)
colnames(miRNA_filtered_normal) <- sapply(colnames(miRNA_filtered_normal), remove_read_count_prefix)

# Remove the 'miRNA_ID' column from both as its a non numeric column and it will create issue with TCGAanalyze_DEA
miRNA_filtered_normal <- miRNA_filtered_normal[, !names(miRNA_filtered_normal) %in% "miRNA_ID"]
miRNA_filtered_tumor <- miRNA_filtered_tumor[, !names(miRNA_filtered_tumor) %in% "miRNA_ID"]

extract_base_sample_id <- function(id) {
  return(substr(id, 1, nchar(id) -9))
}

############################ Refine and remove duplicates ##########################################

# Apply the function to the column names of both data frames
exp_filtered_tumor_base_ids <- sapply(colnames(exp_filtered_tumor), extract_base_sample_id)
miRNA_filtered_tumor_base_ids <- sapply(colnames(miRNA_filtered_tumor), extract_base_sample_id)

# Find common base sample IDs
common_base_samples <- dplyr::intersect(exp_filtered_tumor_base_ids, miRNA_filtered_tumor_base_ids)

# Filter both data frames to only include columns with common base sample IDs
exp_filtered_tumor_common <- exp_filtered_tumor[, exp_filtered_tumor_base_ids %in% common_base_samples]
miRNA_filtered_tumor_common <- miRNA_filtered_tumor[, miRNA_filtered_tumor_base_ids %in% common_base_samples]

# Ensure the columns are in the same order for both data frames
exp_filtered_tumor_common <- exp_filtered_tumor_common[, order(colnames(exp_filtered_tumor_common))]
miRNA_filtered_tumor_common <- miRNA_filtered_tumor_common[, order(colnames(miRNA_filtered_tumor_common))]

# Count frequency of each common base sample ID in miRNA_filtered_tumor_base_ids
miRNA_count <- table(miRNA_filtered_tumor_base_ids[miRNA_filtered_tumor_base_ids %in% common_base_samples])
# Count frequency of each common base sample ID in exp_filtered_tumor_base_ids
exp_count <- table(exp_filtered_tumor_base_ids[exp_filtered_tumor_base_ids %in% common_base_samples])
# Combine the counts into a data frame for easy viewing
combined_counts <- data.frame(
  SampleID = names(miRNA_count),
  miRNA_Count = as.integer(miRNA_count),
  exp_Count = as.integer(exp_count[names(miRNA_count)])
)

duplicates <- combined_counts$SampleID[combined_counts$miRNA_Count > 1 | combined_counts$exp_Count > 1]
exp_filtered_tumor_common <- exp_filtered_tumor_common[, !colnames(exp_filtered_tumor_common) %in% duplicates]
miRNA_filtered_tumor_common <- miRNA_filtered_tumor_common[, !colnames(miRNA_filtered_tumor_common) %in% duplicates]

# Function to check if any of the duplicates are in a given column name
has_duplicate <- function(column_name, duplicates) {
  any(sapply(duplicates, function(duplicate) grepl(duplicate, column_name)))
}

# Apply the function to each column name and negate the result to filter out duplicates
exp_filtered_tumor_common <- exp_filtered_tumor_common[, !sapply(colnames(exp_filtered_tumor_common), has_duplicate, duplicates=duplicates)]
miRNA_filtered_tumor_common <- miRNA_filtered_tumor_common[, !sapply(colnames(miRNA_filtered_tumor_common), has_duplicate, duplicates=duplicates)]

dim(exp_filtered_tumor_common)
dim(miRNA_filtered_tumor_common)

# differential expression
diff_expressed_miRNA <- TCGAanalyze_DEA(
  mat1 = miRNA_filtered_normal,
  mat2 = miRNA_filtered_tumor_common,
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.05 ,
  logFC.cut = 0,
  method = "glmLRT"
)

write.xlsx(diff_expressed_miRNA, file = "diff_expressed_miRNA_LUMB.xlsx", rowNames = TRUE)
gc()

############################# Section 4: Anti correlated gene-mirna pairs ############################

# Extract upregulated genes and downregulated miRNAs
upregulated_genes <- diff_expressed_genes %>% dplyr::filter(logFC > 2)
downregulated_miRNA <- diff_expressed_miRNA %>% dplyr::filter(logFC < -2)

# Extract upregulated genes and downregulated miRNAs
downregulated_genes <- diff_expressed_genes %>% dplyr::filter(logFC < -2)
upregulated_miRNA <- diff_expressed_miRNA %>% dplyr::filter(logFC > 2)

extract_base_sample_id <- function(id) {
  # Assuming the format is always "TCGA-XX-XXXX-01A-XX", ignoring the last 4 characters
  return(substr(id, 1, nchar(id) -9))
}

# ------------------------------------ upregulated genes - Downregulated miRNAs ---------------------------
# Perform correlation analysis
correlated_pairs <- list()
total_operations <- length(rownames(upregulated_genes)) * length(rownames(downregulated_miRNA))
progress_bar <- txtProgressBar(min = 0, max = total_operations, style = 3)
operation_count = 0
for (gene in rownames(upregulated_genes)) {
  for (mirna in rownames(downregulated_miRNA)) {
    operation_count <- operation_count + 1
    setTxtProgressBar(progress_bar, operation_count)
    gene_name <- upregulated_genes[gene, "gene_name"]
    # Run Spearman correlation
    correlation <- cor.test(as.numeric(exp_filtered_tumor_common[gene,]), as.numeric(miRNA_filtered_tumor_common[mirna,]), method = "spearman")
    # Check for the correlation threshold
    if (!is.na(correlation$estimate)) {
      if (correlation$estimate <= -0.4) {
        correlated_pairs[[paste(gene_name, mirna, sep = "_")]] <- c(correlation$estimate, correlation$p.value)
      }
    }
  }
}

close(progress_bar)

# Convert the list to a dataframe
correlated_pairs_df <- do.call(rbind, lapply(names(correlated_pairs), function(x) {
  parts <- unlist(strsplit(x, "_"))
  data.frame(
    gene = parts[1],
    miRNA = parts[2],
    correlation = correlated_pairs[[x]][1],
    p_value = correlated_pairs[[x]][2]
  )
}))

# Write the correlated pairs to an Excel file
write.xlsx(correlated_pairs_df, file = "correlated_pairs_upgene_downmirna_LUMB.xlsx", rowNames = FALSE)

# ------------------------------------ Downregulated genes - Upregulated miRNAs ---------------------------

correlated_pairs <- list()
total_operations <- length(rownames(downregulated_genes)) * length(rownames(upregulated_miRNA))
progress_bar <- txtProgressBar(min = 0, max = total_operations, style = 3)
operation_count = 0
for (gene in rownames(downregulated_genes)) {
  for (mirna in rownames(upregulated_miRNA)) {
    operation_count <- operation_count + 1
    setTxtProgressBar(progress_bar, operation_count)
    gene_name <- downregulated_genes[gene, "gene_name"]
    # Run Spearman correlation
    correlation <- cor.test(as.numeric(exp_filtered_tumor_common[gene,]), as.numeric(miRNA_filtered_tumor_common[mirna,]), method = "spearman")
    # Check for the correlation threshold
    if (!is.na(correlation$estimate)) {
      if (correlation$estimate <= -0.4) {
        correlated_pairs[[paste(gene_name, mirna, sep = "_")]] <- c(correlation$estimate, correlation$p.value)
      }
    }
  }
}

close(progress_bar)

# Convert the list to a dataframe
correlated_pairs_df <- do.call(rbind, lapply(names(correlated_pairs), function(x) {
  parts <- unlist(strsplit(x, "_"))
  data.frame(
    gene = parts[1],
    miRNA = parts[2],
    correlation = correlated_pairs[[x]][1],
    p_value = correlated_pairs[[x]][2]
  )
}))

# Write the correlated pairs to an Excel file
write.xlsx(correlated_pairs_df, file = "correlated_pairs_downgene_upmirna_LUMB.xlsx", rowNames = FALSE)
gc()
############################  Section 5: Enrichment Analysis ###############################################

ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Tumor Vs Normal",
  RegulonList = diff_expressed_genes$gene_name
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = "Biological_Process.pdf",
  GOBPTab = ansEA$ResBP,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)


TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = "Biological_Process_2.pdf",
  GOCCTab = ansEA$ResCC,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = "Molecular_Functions.pdf",
  GOMFTab = ansEA$ResMF,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = "Canonical_pathways.pdf",
  PathTab = ansEA$ResPat,
  nRGTab = rownames(diff_expressed_genes),
  nBar = 20
)

# confirm pathway enrichment

library(SummarizedExperiment)


exp_filtered = cbind(exp_filtered_tumor, exp_filtered_normal)
# DEGs TopTable
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = diff_expressed_genes,
  typeCond1 = "Tumor",
  typeCond2 = "Normal",
  TableCond1 = exp_filtered[,colnames(exp_filtered_tumor)],
  TableCond2 = exp_filtered[,colnames(exp_filtered_normal)]
)

dataDEGsFiltLevel$GeneID <- 0

saveRDS(dataDEGsFiltLevel, "dataDEGsFiltLevel.rds")

library(clusterProfiler)
# Converting Gene symbol to geneID
eg = as.data.frame(
  bitr(
    dataDEGsFiltLevel$mRNA,
    fromType = "ENSEMBL",
    toType = c("ENTREZID","SYMBOL"),
    OrgDb = "org.Hs.eg.db"
  )
)
eg <- eg[!duplicated(eg$SYMBOL),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]

dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$ENSEMBL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[eg$ENSEMBL,]
rownames(dataDEGsFiltLevel) <- eg$SYMBOL

all(eg$SYMBOL == rownames(dataDEGsFiltLevel))
dataDEGsFiltLevel$GeneID <- eg$ENTREZID
dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID

library(pathview)

hsa04670 <- pathview::pathview(
  gene.data  = genelistDEGs,
  pathway.id = "hsa04670",
  species    = "hsa"
)
gc()
############################# Section 6: Methylation - Gene targets ####################################

group.col <- "sample_type"
group1 <- "Primary Tumor"
group2 <- "Solid Tissue Normal"

query_met_tumor <- GDCquery(
  project = cancer_type,
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value",
  sample.type = c("Primary Tumor")
)

GDCdownload(query_met_tumor)
met_tumor <- GDCprepare(query_met_tumor, save = FALSE)

query_met_normal <- GDCquery(
  project = cancer_type,
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value",
  sample.type = c("Solid Tissue Normal")
)

GDCdownload(query_met_normal)
met_normal <- GDCprepare(query_met_normal, save = FALSE)

################ subtype selection #####################################
clinical_data <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
metadataObj <- colData(met_tumor) # contains dataframe with PAM50 tumor classifications and other clinical data for each sample
LumB_samples <- metadataObj[!is.na(metadataObj$paper_BRCA_Subtype_PAM50) & metadataObj$paper_BRCA_Subtype_PAM50 == "LumB",]
# It is important to note that sample ids in exp_tumor is large ids with extra characters and samples$bcr_patient_barcode
# has trunctued id format. Therefore to match these and pull subtype sample expression data, we need following

matched_columns <- sapply(colnames(met_tumor), function(column_name) {
  any(sapply(LumB_samples$bcr_patient_barcode, function(truncated_id) {
    grepl(truncated_id, column_name)
  }))
})

met_tumor_filtered <- met_tumor[, matched_columns]
common_met <- dplyr::intersect(rownames(met_tumor_filtered), rownames(met_normal))

# Subset each dataset to include only common genes
met_filtered_tumor_common <- met_tumor_filtered[common_met, ]
met_filtered_normal_common <- met_normal[common_met, ]

# Subset tumor dataset for common CpG islands
met_filtered_tumor_common <- met_filtered_tumor_common[rownames(met_filtered_tumor_common) %in% rownames(met_filtered_normal_common), ]

# Subset normal dataset for common CpG islands
met_filtered_normal_common <- met_filtered_normal_common[rownames(met_filtered_normal_common) %in% rownames(met_filtered_tumor_common), ]
exp_combined <- cbind(exp_filtered_tumor_common, exp_filtered_normal_common)
assay_tumor <- assay(met_filtered_tumor_common)
assay_normal <- assay(met_filtered_normal_common)
met <- cbind(assay_tumor, assay_normal)
distal.probes <- get.feature.probe(
  genome = "hg19",
  met.platform = "450K"
)

mae <- createMAE(
  exp = assay(exp_combined),
  met = met,
  save = TRUE,
  linearize.exp = TRUE,
  save.filename = "mae.rda",
  filter.probes = distal.probes,
  met.platform = "450K",
  genome = "hg19",
  TCGA = TRUE
)

saveRDS(mae, file = "mae.rds")

# Note: Check colData(mae_tumor)$sample_type <- contain Tumor and Solid tissue normal. Useful to supply in diff analysis

message("Get diff methylated probes")
cores <- 8
group.col <- "sample_type"
group1 <- "Primary Tumor"
group2 <- "Solid Tissue Normal"


direction <- "hypo" #change this for hyper methylated case <---------
dir.out <- paste0("elmer/",direction)
dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)
Sig.probes_hypo <- get.diff.meth(
  data = mae,
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  minSubgroupFrac = 1.0,
  sig.dif = 0.2,
  diff.dir = direction, # Search for hypomethylated probes in group 1
  cores = cores,
  dir.out = dir.out,
  pvalue = 0.05
)

message("Get nearby genes")
nearGenes_hypo <- GetNearGenes(
  data = mae,
  numFlankingGenes = 4,
  probes = Sig.probes_hypo$probe
)

write.xlsx(nearGenes_hypo, file = "nearGenes_hypo.xlsx")

message("Get anti correlated probes-genes")
pair_hypo <- get.pair(
  data = mae,
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  nearGenes = nearGenes_hypo,
  mode = "supervised",
  minSubgroupFrac = 1, # % of samples to use in to create groups U/M
  raw.pvalue = 0.05,   # defualt is 0.001
  Pe = 0.05, # Please set to 0.001 to get significant results
  filter.probes = TRUE, # See preAssociationProbeFiltering function
  filter.percentage = 0.05,
  save = FALSE, # Create CVS file
  filter.portion = 0.2,
  dir.out = dir.out,
  diff.dir = direction,
  cores = cores,
  label = direction
)

write.xlsx(pair_hypo, file = "pair_hypo.xlsx")

#-------------------------------------------------------------
# Step 3.3: Motif enrichment analysis on the selected probes |
#-------------------------------------------------------------
enriched.motif_hypo <- get.enriched.motif(
  data = mae,
  probes = pair_hypo$Probe,
  dir.out = dir.out,
  label = direction,
  pvalue = 0.05, # default is FDR < 0.05
  min.incidence = 10,
  lower.OR = 1.1
)

motif.enrichment_hypo <- read.csv(paste0(dir.out,"/getMotif.",direction, ".motif.enrichment.csv"))
write.xlsx(motif.enrichment_hypo %>% gt::gt(), file = "motif.enrichment_hypo.xlsx")

#-------------------------------------------------------------
# Step 3.4: Identifying regulatory TFs                        |
#-------------------------------------------------------------
TF_hypo <- get.TFs(
  data = mae,
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  mode = "supervised",
  enriched.motif = enriched.motif_hypo,
  dir.out = dir.out,
  cores = cores,
  save.plots = FALSE,
  diff.dir = direction,
  label = direction
)

TF.meth.cor_hypo <- get(load(paste0(dir.out,"/getTF.",direction,".TFs.with.motif.pvalue.rda")))
saveRDS(TF.meth.cor_hypo, file = "TF.meth.cor_hypo.rds")

# Step 1: Filter out rows where all values are greater than 0.05
filtered_rows_TF_meth_cor_hypo <- TF.meth.cor_hypo[apply(TF.meth.cor_hypo, 1, function(x) any(x <= 0.05)), ]

# Step 2: Filter out columns where all values are greater than 0.05
final_filtered_TF_meth_cor_hypo <- filtered_rows_TF_meth_cor_hypo[, colSums(filtered_rows_TF_meth_cor_hypo <= 0.05) > 0]
final_filtered_TF_meth_cor_hypo <- t(final_filtered_TF_meth_cor_hypo)

write.csv(final_filtered_TF_meth_cor_hypo, "filtered_TF.meth.cor_hyper.csv", row.names = TRUE)

# Step 1: Create a data frame from enriched.motif_hypo
motif_data_hypo <- stack(enriched.motif_hypo)
colnames(motif_data_hypo) <- c("Probe", "Motif")

# Step 2: Join pair_hypo with motif_data based on Probe
combined_data_hypo <- pair_hypo %>%
  dplyr::select(Probe, Symbol) %>%
  inner_join(motif_data_hypo, by = "Probe")

result_hypo <- combined_data_hypo %>%
  group_by(Symbol) %>%
  summarize(Motifs = paste(unique(Motif), collapse = ", "))

write.csv(result_hypo, "gene_Motif_hypo.csv", row.names = TRUE) # Final gene- differentially enriched TF methylation site pairs

motif.enrichment.plot(
  motif.enrichment = motif.enrichment_hypo,
  save = FALSE,
  significant = list(lowerOR = 1.2)
)

pdf("TF_hypo.pdf")
scatter.plot(
  data = mae,
  category = group.col,
  save = FALSE,
  lm_line = TRUE,
  byTF = list(
    TF = unlist(stringr::str_split(TF_hypo[1,"top_5percent_TFs"],";"))[1:4],
    probe = enriched.motif_hypo[[TF_hypo$motif[1]]]
  )
)
dev.off()

gc()

direction <- "hyper" # <--------- ##################################
dir.out <- paste0("elmer/",direction)
dir.create(dir.out, showWarnings = FALSE, recursive = TRUE)
message("Get diff methylated probes")
Sig.probes_hyper <- get.diff.meth(
  data = mae,
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  minSubgroupFrac = 1.0,
  sig.dif = 0.2,
  diff.dir = direction, # Search for hypermethylated probes in group 1
  cores = cores,
  dir.out = dir.out,
  pvalue = 0.05
)

datatable(
  data = Sig.probes_hyper[1:10,],
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
  rownames = TRUE
)
message("Get nearby genes")
nearGenes_hyper <- GetNearGenes(
  data = mae,
  numFlankingGenes = 4, # default is 20 genes
  probes = Sig.probes_hyper$probe
)

message("Get anti correlated probes-genes")
pair_hyper <- get.pair(
  data = mae,
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  nearGenes = nearGenes_hyper,
  mode = "supervised",
  minSubgroupFrac = 1,
  raw.pvalue = 0.05,
  Pe = 0.05,
  filter.probes = TRUE,
  filter.percentage = 0.05,
  save = FALSE,
  filter.portion = 0.2,
  dir.out = dir.out,
  diff.dir = direction,
  cores = cores,
  label = direction
)

write.xlsx(pair_hyper, file = "pair_hyper.xlsx")

enriched.motif_hyper <- get.enriched.motif(
  data = mae,
  probes = pair_hyper$Probe,
  dir.out = dir.out,
  label = direction,
  pvalue = 0.05, # default is FDR < 0.05
  min.incidence = 10,
  lower.OR = 1.1
)

motif.enrichment_hyper <- read.csv(paste0(dir.out,"/getMotif.",direction, ".motif.enrichment.csv"))
write.xlsx(motif.enrichment_hyper %>% gt::gt(), file = "motif.enrichment_hyper.xlsx")

TF_hyper <- get.TFs(
  data = mae,
  group.col = group.col,
  group1 = group1,
  group2 =  group2,
  mode = "supervised",
  enriched.motif = enriched.motif_hyper,
  dir.out = dir.out,
  cores = cores,
  save.plots = FALSE,
  diff.dir = direction,
  label = direction
)

TF.meth.cor_hyper <- get(load(paste0(dir.out,"/getTF.",direction,".TFs.with.motif.pvalue.rda")))
write.csv(TF.meth.cor_hyper, "TF_meth_cor_hyper.csv", row.names = TRUE)

# Step 1: Filter out rows where all values are greater than 0.05
filtered_rows_TF_meth_cor_hyper <- TF.meth.cor_hyper[apply(TF.meth.cor_hyper, 1, function(x) any(x <= 0.05)), ]

dim(filtered_rows_TF_meth_cor_hyper)
# Step 2: Filter out columns where all values are greater than 0.05
final_filtered_TF_meth_cor_hyper <- filtered_rows_TF_meth_cor_hyper[, colSums(filtered_rows_TF_meth_cor_hyper <= 0.05) > 0]
final_filtered_TF_meth_cor_hyper <- t(final_filtered_TF_meth_cor_hyper)
write.csv(final_filtered_TF_meth_cor_hyper, "filtered_TF.meth.cor_hyper.csv", row.names = TRUE)

# Step 1: Create a data frame from enriched.motif_hypo
motif_data_hyper <- stack(enriched.motif_hyper )
colnames(motif_data_hyper) <- c("Probe", "Motif")

# Step 2: Join pair_hypo with motif_data based on Probe
combined_data_hyper <- pair_hyper %>%
  dplyr::select(Probe, Symbol) %>%
  inner_join(motif_data_hyper, by = "Probe")

result_hyper <- combined_data_hyper %>%
  group_by(Symbol) %>%
  summarize(Motifs = paste(unique(Motif), collapse = ", "))


write.csv(result_hyper, "gene_Motif_hyper.csv", row.names = TRUE)

