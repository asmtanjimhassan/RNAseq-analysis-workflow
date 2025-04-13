# importing libraries to R script
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(dplyr)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations


# importing dataset
result <- read.csv(file = "datasets/pea/GSE247894.csv")

# summarizing the data
result <- result %>%
  dplyr::select(2,5,6,7) %>%
  dplyr::rename(log2fc = log2FoldChange) %>%
  dplyr::rename(gene = GeneName)

# preparing the path variable
list.files
input <- "datasets/pea/"
output <- "output/pea/"
bg_path <- "output/pea/"

# Functions ------------------------------------------------
## Function: Adjacency matrix to list ####
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}



# Preparing the background genes ######################################

# Get the genes that are present in your dataframe
genes_in_data <- result$gene

# Read in the .gmt file
file <- "datasets/pea/msigdb_v2024_1_Hs_GMTs/c2.cp.reactome.v2024.1.Hs.symbols.gmt"
pwl2 <- read.gmt(file) 

# Subset to the genes that are present in our dataset
pwl2 <- pwl2[pwl2$gene %in% genes_in_data,] 

# Save the filtered background gene set
saveRDS(pwl2, "output/pea/msigdb_Hs.RDS")



# annotate the differentially expressed gene ##########################
result <- result %>%
  mutate(diffexp = case_when(
    log2fc > 0 & padj < 0.05 ~ 'up',
    log2fc < 0 & padj < 0.05 ~ 'down',
    padj > 0.05 ~ 'no'
    )
  )

# removing the non-significant genes
result <- result[result$diffexp !='no',]

# Substitute names so they are annotated nicely in the heatmap later
result$diffexp <- gsub('down', 'Healthy', gsub('up', 'Severe', result$diffexp))
unique(result$diffexp)

# split the dataset into a list of sub-dataframe 
# i.e. upregulated and downregulated genes
deg_result_list <- split(result, result$diffexp)


## Run ClusterProfiler ###########################################


# Settings
name_of_comparison <- 'severevshealthy'
background_genes <- 'msigdb'
bg_genes <- readRDS(paste0(bg_path, 'msigdb_Hs.RDS')) # read in the background genes
padj_cutoff <- 0.05 # p-adjusted threshold, used to filter out pathways
genecount_cutoff <- 5 # minimum number of genes in the pathway, used to filter out pathways



# Run clusterProfiler on each sub-dataframe
res <- lapply(names(deg_result_list),
              function(x) enricher(gene = deg_result_list[[x]]$gene,
                                   TERM2GENE = bg_genes))

names(res) <- names(deg_result_list)


#Convert the list of enrichResults for each sample_pattern to a dataframe with the pathways
res_df <- lapply(names(res), function(x) rbind(res[[x]]@result))
names(res_df) <- names(res)
res_df <- do.call(rbind, res_df)
head(res_df)


res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                            diffexp = gsub('\\.GOBP.*$|\\.KEGG.*$|\\.REACTOME.*$', '', rownames(res_df)))


# Subset to those pathways that have p adj < cutoff and gene count > cutoff 
# (you can also do this in the enricher function)

# select only target pathways have p adjusted < 0.05 and at least 6 genes
target_pws <- unique(res_df$ID[res_df$p.adjust < padj_cutoff & res_df$Count > genecount_cutoff]) 
res_df <- res_df[res_df$ID %in% target_pws,]

print('Saving clusterprofiler results')
write.csv(res_df, 'output/pea/enr_analysis_result.csv', row.names = FALSE)
