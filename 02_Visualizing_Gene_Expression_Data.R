# importing libraries to R script
library(tidyverse)
library(dplyr)
library(ggplot2)

# importing dataset to R script
data <- read.csv(file = "output/gene_expression_data_summarized.csv")


# bar chart
barchart <- data %>%
  filter(gene == "BRCA2") %>%
  ggplot(.,
         aes(
           x = samples,
           y = FPKM,
           fill = tissue)
         ) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Exporting the Bar Chart
ggsave(barchart, filename="output/02_Bar_Chart.pdf", width = 14, height = 6)


# density plot
density_plot <- data %>%
  filter(gene == "BRCA1") %>%
  ggplot(.,
         aes(
           x = FPKM,
           fill = tissue))+
  geom_density(alpha = 0.6)

# Exporting the Density Plot
ggsave(density_plot, filename="output/02_Density_Plot.pdf", width = 8, height = 6)


# Box Plot
box_plot <- data %>%
  filter(gene == "BRCA1") %>%
  ggplot(.,
         aes(
           x = metastasis,
           y = FPKM))+
  geom_boxplot()

# Exporting the Box Plot
ggsave(box_plot, filename="output/02_Box_Plot.pdf", width = 8, height = 6)

#Scatter Plot
scatter_plot <- data %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  tidyr::spread(key = gene, value = FPKM) %>%
  ggplot(.,
         aes(
           x = BRCA1,
           y = BRCA2
           
         )
  ) +
  geom_point(aes(color = tissue)) +
  geom_smooth(method = "lm", se = FALSE)

# Exporting the Box Plot
ggsave(scatter_plot, filename="output/02_Scatter_Plot.pdf", width = 8, height = 6)


# Heatmap
heatmap <- data %>%
  filter(gene %in% c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')) %>%
  ggplot(.,
         aes(
           x = samples,
           y = gene,
           fill = FPKM
         )
  )+
  geom_tile() +
  scale_fill_gradient(low="white", high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(heatmap, filename="Output/lec02_Heat_Map.pdf", width = 14, height = 6)
