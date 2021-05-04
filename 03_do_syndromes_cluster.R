library(plyr)
library(ggplot2)
library(ggfortify)
library(cluster)
library(dendextend)
library(factoextra)

rm(list=ls())
dev.off()

# Load and merge groups
files <- list.files(path = "data/simulated_group/", pattern= "*.Rds", full.names = TRUE)
groups <- lapply(files, readRDS)
all_data <- ldply(groups, rbind)

# Simple PCA with all individual per syndrome
filename_pca <- "figures/pca_all_groups.png"

pca <- prcomp(all_data[,-1:-2],
              center = TRUE,
              scale = TRUE)

autoplot(pca, 
         data = all_data, 
         colour = "syndrome",
         label = TRUE, 
         label.size = 2,
         size = 3) +
  labs(title = "All groups",
       color = "Syndrome") +
  theme(panel.grid.major = element_line(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#ffffff",color = NA),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 0, color = "#4e4d47"),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size=5)))

dev.copy(png, filename_pca, width=6, height=6, units="in", res=500)
dev.off()

# PCA with variables' vectors
filename_pca_vectors <- "figures/pca_all_groups_vectors.png"
fviz_pca_biplot(pca,
                repel=FALSE,
                geom.ind = "point",
                fill.ind = all_data$syndrome,
                pointshape = 21 ,
                pointsize = 2,
                alpha.ind=0.5,
                alpha.var=1, 
                col.var="contrib",
                legend.title = list(fill = "Syndrome"),
                addEllipses = TRUE) +
  labs(title = "PCA All Groups") +
  scale_color_continuous(guide = 'none') +
  theme(panel.grid.major = element_line(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "#ffffff",color = NA),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14, hjust = 0, color = "#4e4d47"),
        legend.position = "right")

dev.copy(png, filename_pca_vectors, width=6, height=6, units="in", res=500)
dev.off()

# Testing the influence of the clustering method
## Method 1: Ward's agglomerative hierarchical clustering
test01 <- agnes(all_data[,-1:-2], method = "ward", diss = FALSE, metric = "euclidean")
plot(test01, which.plot = 2, main = "Ward???s agglomerative hierarchical clustering")


## Method 2: Ward's agglomerative hierarchical clustering
test02_diss <- dist(all_data[,-1:-2], method = "euclidean")
  
test02_ward <- hclust(test02_diss, method = "ward.D")
test02_average <- hclust(test02_diss, method = "average")

dend_02_ward <- as.dendrogram (test02_ward)
dend_02_average <- as.dendrogram (test02_average)

tanglegram(dend_02_ward, dend_02_average)

