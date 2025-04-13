# importing libraries to R script
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(GEOquery)
library(apeglm)


# read in counts data
counts_data <- read.table(file = "datasets/deseq2/GSE273142_Processed_data.tab")

counts_data <- counts_data %>%
  rename(GSM8422516 = JG.27250_S1.rsem.genes.results) %>%
  rename(GSM8422517 = JG.27251_S2.rsem.genes.results) %>%
  rename(GSM8422518 = JG.27255_S6.rsem.genes.results) %>%
  rename(GSM8422519 = JG.27256_S7.rsem.genes.results) %>%
  rename(GSM8422520 = JG.27260_S11.rsem.genes.results) %>%
  rename(GSM8422521 = JG.27261_S12.rsem.genes.results)


gse <- getGEO(GEO = "GSE273142", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))   #dim(metadata)

colData <- metadata %>%
  select(12) %>%
  rename(condition = characteristics_ch1.2) %>%
  mutate(condition = gsub("genotype: ", "", condition)) %>%
  mutate(condition = gsub(" ", "_", condition))

# read in sample info
colData <- read.csv(file = "datasets/deseq2/sample_info.csv")
head(colData)

# making sure the row names in colData matches to column names in Counts_data
all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data) == rownames(colData))

# Construct a DESeq Dataset Object
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_data),
  colData = colData,
  design = ~ condition
)

# pre-filtering : removing rows with low gene count
# keeping rows that have atleast 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the factor level
dds$condition <- factor(dds$condition, levels = c("DAPK3_WT", "DAPK3_KO"))
dds$condition <- droplevels(dds$condition)

# run DESeq : perform statistical test to identify the differentially expressed genes
dds <- DESeq(dds)
deseq_result <- results(dds)


# explore result
# default resut with with p-value < 0.1
deseq_result
summary(deseq_result)

# Result as dataframe of R
deseq_result_data <- as.data.frame(deseq_result)
head(deseq_result_data)

# re-order by p-value
deseq_result_reordered <- deseq_result_data[order(deseq_result_data$padj),]
head(deseq_result_reordered)


deseq_result_data["ENSG00000167657_DAPK3",]


# p-value of interest
result0.01 <- results(dds, alpha = 0.01)
summary(result0.01)




# constrasts
resultsNames(dds)
result_contrast <- results(dds, contrast = c("condition", "DAPK3_KO", "DAPK3_WT"))
result_contrast



# Shrinkage of log2foldchange
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_DAPK3_KO_vs_DAPK3_WT", type="apeglm")
resLFC_reordered <- resLFC[order(resLFC$padj),]

# Data Visualization
# MA plot
ma_plot <- plotMA(resLFC)

# Plot count of each gene
plotCounts(dds, gene="ENSG00000080824_HSP90AA1", intgroup="condition")



