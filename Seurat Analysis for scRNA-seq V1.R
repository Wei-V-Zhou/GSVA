
# Clear the environment
rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# Load libraries
{
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(tidyverse)
  library(reticulate)
}

# Load data
if (!file.exists("m4T1.scset.Rdata")){
  dat <- read.table("m4T1.raw.featurecounts.txt", header = T)
  ann <- read.csv("m4T1.cell.annotation.csv", header = T)
  counts <- dat[ , 5:ncol(dat)]
  colnames(counts) <- ann[ , 1]
  save(counts, ann, file = "m4T1.scset.Rdata")
}
load("m4T1.scset.Rdata")

# Choose early-stage cells for analysis
ann.D4 <- ann[ann$CellType == "BoneMet_D4", ]
{
  ann.D10 <- ann[ann$CellType == "BoneMet_D10", ]
  ann.D16 <- ann[ann$CellType == "BoneMet_D16", ]
  ann.CL <- ann[ann$CellType == "CellLine", ]
  ann.PT <- ann[ann$CellType == "PrimaryTumor", ]
}
counts_D4 <- counts[ , ann.D4[ , 1]]

# Remove mCherry = 0
count <- counts_D4[ , - which((counts["mCherry", ] != 0) == FALSE)]
ann.D4 <- ann.D4[- which((counts["mCherry", ] != 0) == FALSE), ]
count <- count[- c(which(rownames(count) == "mCherry"), which(rownames(count) == "luciGFP")), ]

# Create a Seurat Object
scset <- CreateSeuratObject(
  counts = as.matrix(count),
  meta.data = ann.D4,
  project = "4T1.scRNA",
  min.cells = 3,
  min.features = 200)

# Cell QC
scset[["percent.mt"]] <- PercentageFeatureSet(scset, pattern = "^mt-")
VlnPlot(scset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(scset, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
scset <- subset(scset, subset = nFeature_RNA > 5000 & nFeature_RNA < 12500 & percent.mt < 5)
dim(scset)

# Normalizing the data
scset <- NormalizeData(scset, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
scset <- FindVariableFeatures(scset, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(scset), 10)
# Plot varible features with and without labels
plot1 <- VariableFeaturePlot(scset)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))

# Scaling the data
all.genes <- rownames(scset)
scset <- ScaleData(scset, features = all.genes)

# Perform linear dimensional reduction
scset <- RunPCA(scset, features = VariableFeatures(object = scset))
print(scset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scset, dims = 1:2, reduction = "pca")
DimPlot(scset, reduction = "pca")

# Determine the dimensionality of the dataset roughly
# Way 1: DimHeatmap function: Link the GSEA analysis
DimHeatmap(scset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(scset, dims = 1:15, cells = 500, balanced = TRUE)
# Way 2: JackStrawPlot function: Significant PCs show a strong enrichment of features.
# Find the sharp drop-off in significance after the first X-Y PCs
scset <- JackStraw(scset, num.replicate = 100)
scset <- ScoreJackStraw(scset, dims = 1:20)
JackStrawPlot(scset, dims = 1:15)
# Way 3: ElbowPlot function: A heuristic method generates an "ElbowPlot"
# Observe the elbow around PCX-Y and the majority of true signal is captured at pcY
ElbowPlot(scset)

# The gold standard to decide the dimensionality:
# 1. the culmutative contribution of PCs > 90%
# 2. the variance contribution of PCs < 5%
# 3. the differience variance of neighborhood PCs < 0.1%
# The method to get the accurate PCs:
# Determine the percent of variation associated with each PC
pct <- scset[["pca"]]@stdev / sum(scset[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
# Determine the PCs exihibit cumulative percent greater than 90% and variance less than 5%
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variance of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] +1
# Conclude the number of PCs
pcs <- min(co1, co2)
# Create a dataframe with values
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
# ElbowPlot for visualization
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() +
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


# Cluster the cells
scset <- FindNeighbors(scset, dims = 1:13)
scset <- FindClusters(scset, resolution = 0.5)
head(Idents(scset), 6)

# Run non-linear dimensional reduction(UMAP/tSNE)
scset <- RunUMAP(scset, dims = 1:13)
DimPlot(scset, reduction = "umap", pt.size = 1)
scset <- RunTSNE(scset, dims = 1:13)
DimPlot(scset, reduction = "tsne")

# Assign cell type to cluster
current.cluster.ids <- c(0, 1, 2, 3, 4)
new.cluster.ids <- c("Stromal cell", "B cell", "Neutrophil", "Macrophage", "Erythroid cell")
scset@active.ident <- plyr::mapvalues(x = scset@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(object = scset, reduction = "umap", do.label = TRUE, pt.size = 1)

# Save the data
save(scset, file = "m4T1.new.featurecounts.RData")

#============================#
#       Musician: Resonance  #
#           Date: 2019/09/06 #
# Revised author: Resonance  #
#           Time: 2019/11/02 #
#============================#