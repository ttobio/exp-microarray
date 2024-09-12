#Heatmap
library(tidyverse)
exp <- read.csv("agg.csv")
degs <- read.csv("DEGS_lfc1.5.csv")
ex_degs <- exp %>%
  filter(exp$X %in% degs$X)
rownames(ex_degs) <- ex_degs$X
ex_degs <- ex_degs[,-1]
head(ex_degs)

library(pheatmap)
p <- read.csv("pheno_data.csv")
my_sample_col <- data.frame(sample = p$status)
row.names(my_sample_col) <- colnames(ex_degs)
pheatmap(ex_degs, annotation_col = my_sample_col)


my_hclust_gene <- hclust(dist(ex_degs), method = "complete")

# install if necessary
#install.packages("dendextend")

library(dendextend)

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)

my_gene_col

my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))

head(my_gene_col)

my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster1", no = "cluster2"))
row.names(my_sample_col) <- colnames(ex_degs)

my_colour = list(
  sample = c(healthy = "#e89829", HCC = "#cc4ee0"),
  cluster = c(cluster1 = "blue", cluster2 = "red")
)
p_color <- colorRampPalette((c("red", "black", "green")))(50)
p <- pheatmap(ex_degs,
              color = p_color,
              annotation_colors = my_colour,
              annotation_row = my_gene_col,
              annotation_col = my_sample_col,
              cutree_rows = 4,
              cutree_cols = 2)

save_pheatmap_png <- function(x, filename, width=1200, height=1700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(p, "heatmap_new.png")
#===============================================================================
#Volcano plot
library(EnhancedVolcano)
RES <- read.csv("RES.csv")
volcano_plot <- EnhancedVolcano(RES, x = "logFC", y = "adj.P.Val",title = "HCC VS Healthy", lab= RES$X, pCutoff = 1e-2, FCcutoff = 0.5,
                                pointSize = 3.0,
                                labSize = 6.0)

print(volcano_plot)
png(filename = "volcano_plot_new.png", width = 2400, height = 2400, res = 200)  # Set output details

plot(volcano_plot)  

dev.off()