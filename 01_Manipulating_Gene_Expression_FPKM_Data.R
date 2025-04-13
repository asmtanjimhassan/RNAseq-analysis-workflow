# importing libraries to R script
library(tidyverse)
library(dplyr)
library(GEOquery)

# importing dataset to R script
data <- read.csv(file = "datasets/GSE183947_fpkm.csv")


# get metadata
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
metadata <- pData(phenoData(gse[[1]]))   #dim(metadata)
view(metadata)

# summarizing the metadata
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

# summarizing the GSE data
data.long <- data %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)
  #head()

# merging both dataframe
data.long <- data.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

# manipulating dataset
data_summarized <- data.long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  group_by(gene, tissue) %>%
  summarize(
    mean_FPKM = mean(FPKM),
    median_FPKM = median(FPKM),
    std_FPKM = sd(FPKM)
  )

view(data_summarized)

# Export to CSV File
write.csv(data.long,"output/gene_expression_data_summarized.csv", row.names = FALSE)
