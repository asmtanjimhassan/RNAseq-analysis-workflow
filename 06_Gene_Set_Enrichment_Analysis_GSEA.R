# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(fgsea)


# importing dataset
result <- read.csv(file = "datasets/pea/GSE247894.csv")

# summarizing the data
result <- result %>%
  dplyr::select(2,5,6,7) %>%
  dplyr::rename(log2fc = log2FoldChange) %>%
  dplyr::rename(gene = GeneName)



# Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# Function: prepare_gmt
prepare_gmt <- function(gmt_file, gene_in_data, savefile = FALSE){
  
  #read in gmt file
  gmt <- gmtPathways("datasets/pea/msigdb_v2024_1_Hs_GMTs/c2.cp.reactome.v2024.1.Hs.symbols.gmt")
  hidden <- unique(unlist(gmt))
  
  # convert gmt file to matrix with the genes as rows and
  # for each GO annotations the value are either 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  # subset to the gene that are present in our data to avoid bias
  hidden1 <- intersect(gene_in_data, hidden)
  
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,]) > 5)]]
  
  final_list <- matrix_to_list(mat)

  print("Shit is done")
  return(final_list)
}


# prepare background genes
bg_path = "datasets/pea/msigdb_v2024_1_Hs_GMTs/"

gene_in_data <- result$gene
list.files(bg_path)
gmt_file <- list.files(path = bg_path, pattern = ".gmt", full.names = TRUE)
bg_genes <- prepare_gmt(gmt_file[16], gene_in_data, savefile = FALSE)


# 2. Prepare your ranked list of genes ##########################
# we will use the signed p values from spatial DGE as ranking
rankings <- sign(result$log2fc)*(-log10(result$pvalue))
names(rankings) <- result$gene # genes as names#

rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

max(rankings)
min(rankings)

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking

ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings.
                 # if only pos or neg, set to 'pos', 'neg'
                 minSize = 10,
                 maxSize = 500,
                 nproc = 1) # for parallelisation


## 6. Check results ------------------------------------------------------
# Top 6 enriched pathways (ordered by p-val)
head(GSEAres[order(pval), ])

sum(GSEAres[, padj < 0.01])
sum(GSEAres[, pval < 0.01])

number_of_top_pathways_up <- 10
number_of_top_pathways_down <- 10
topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(bg_genes[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()


# Select only independent pathways, removing redundancies/similar pathways
collapsedPathways <- collapsePathways(GSEAres[order(pval)][pval < 0.05], bg_genes, rankings)
mainPathways <- GSEAres[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway]
#pdf(file = paste0('GSEA/Selected_pathways/', paste0(filename, background_genes, '_gsea_mainpathways.pdf')), width = 20, height = 15)
plotGseaTable(bg_genes[mainPathways], rankings, GSEAres, gseaParam = 0.5)
#dev.off()


# plot the most significantly enriched pathway
plotEnrichment(bg_genes[[head(GSEAres[order(padj), ], 1)$pathway]],
               rankings) + 
  labs(title = head(GSEAres[order(padj), ], 1)$pathway)



