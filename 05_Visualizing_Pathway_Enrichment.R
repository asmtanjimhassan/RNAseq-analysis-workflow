# importing libraries to R script
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(dplyr)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library(ggplot2)

# import dataset
res_df <- read.csv(file = "output/pea/enr_analysis_result.csv")
bg_genes <- readRDS("output/pea/msigdb_Hs.RDS")

# and the deseq result data frame

# importing dataset
deseq.result <- read.csv(file = "datasets/pea/GSE247894.csv")

# summarizing the data
deseq.result <- deseq.result %>%
  dplyr::select(2,5,6,7) %>%
  dplyr::rename(log2fc = log2FoldChange) %>%
  dplyr::rename(gene = GeneName)

# Convert clusterProfiler object to a new "enrichResult" object
# Select only upregulated genes in Severe
res_df <- res_df %>%
  filter(diffexp == 'Healthy') %>% 
  dplyr::select(!c('minuslog10padj', 'diffexp')) %>%
  mutate(ID = gsub("REACTOME_", "", ID)) %>%
  mutate(Description = gsub("REACTOME_", "", Description))

rownames(res_df) <- res_df$ID

enrichres <- new("enrichResult",
                 readable = FALSE,
                 result = res_df,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "human",
                 ontology = "UNKNOWN",
                 gene = deseq.result$gene,
                 keytype = "UNKNOWN",
                 universe = unique(bg_genes$gene),
                 gene2Symbol = character(0),
                 geneSets = bg_genes)
class(enrichres)


# Barplot
barplot(enrichres, showCategory = 20, font.size = 3)
mutate(enrichres, qscore = -log(p.adjust, base = 10)) %>% 
  barplot(x = "qscore")

# Dotplot
dotplot(enrichres, showCategory = 15) + ggtitle("Severe vs Healthy")

# Cnetplot
cnetplot(enrichres)

# Heatplot
heatplot(enrichres, showCategory = 3)

# Treeplot
# calculate pairwise similarities of the enriched terms using Jaccard’s similarity index
enrichres2 <- pairwise_termsim(enrichres)
treeplot(enrichres2)

# Enrichment map 
emapplot(enrichres2,
         cex_label_category = 0.4)

# Upsetplot
upsetplot(enrichres)
