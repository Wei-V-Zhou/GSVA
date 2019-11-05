
rm(list=ls())
gc()
options(stringsAsFactors=FALSE)

# Load libraries
{
  library(ggplot2)
  library(GSEABase)
  library(GSVA)
  library(gplots)
  library(pheatmap)
  library(gdata)
  library(geneplotter)
}

# Load the method of phylogenetic tree
{
  cor.dist <- function (x) 
  {
    as.dist(1-cor(t(x), use="pairwise.complete.obs",method="kendall"))
  }
  cor.dist.spearman <- function (x) 
  {
    as.dist(1-cor(t(x), use="pairwise.complete.obs",method="spearman"))
  }
  dist.manhattan <- function(x){dist(x, method="manhattan")}
  hclust.ward.d <- function(x){hclust(x, method="ward.D")}
  hclust.ward.d2 <- function(x){hclust(x, method="ward.D2")}
}


# Read Hallmark Gene Sets, please obtain it from original authors' website: MSigDB
hm.sets <- getGmt("h.all.v7.0.symbols.gmt")

# Load primary dataset
if(!file.exists("m4T1.primary.Rdata")){
  load("m4T1.D4.sceset.qc.new.Rdata")
  exprs <- log2(calculateCPM(counts(sce.all)) +1)
  exprs_norm <- t(scale(t(exprs)))
  save(ann, exprs_norm, file = "m4T1.primary.Rdata")
}
load("m4T1.primary.Rdata")

# Isolate D4 breast cancer cells
ann.cellname_D4 <- ann[ann$CellType=="BoneMet_D4", ]
exprs_D4 <- exprs_norm[ , colnames(exprs_norm)==row.names(ann.cellname_D4)]
rownames(exprs_D4) <- toupper(rownames(exprs_D4))

{
  # Choose SDC4 cluster
  marker.group <- "SDC4"
  marker.exprs <- exprs_D4[rownames(exprs_D4)==marker.group, ]
  marker.data <- rbind(marker.exprs, t(ann.cellname_D4))
  colnames(marker.data) <- marker.data[nrow(marker.data), ]
  expr.D4 <- as.matrix(marker.data[ , colnames(marker.data)=="BoneMet_D4"])
  {
    colnames(expr.D4) <- expr.D4[(nrow(expr.D4)-1), ]
    expr.D4.1 <- expr.D4[c(-nrow(expr.D4),-(nrow(expr.D4)-1)), ]
  }
  Stromal.PT_D4 <- t(expr.D4.1)
  row.names(Stromal.PT_D4) <- marker.group
  anno <- ann.cellname_D4[colnames(Stromal.PT_D4), ]
  anno$CellType <- "Stromal"
  anno_SDC4 <- colnames(Stromal.PT_D4)[Stromal.PT_D4 > 0]
  exprs_D4 <- exprs_D4[ , anno_SDC4]
}

# Run GSVA on the Breast Cancer data using Hallmark Gene Sets
datasets.hm <- gsva(as.matrix(exprs_D4), hm.sets, method="gsva")
rownames(datasets.hm) <- gsub("HALLMARK_","",rownames(datasets.hm))

# Generate GSVA analysis
pheatmap(datasets.hm, treeheight_col = 30, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2", 
         scale = "column",color = colorRampPalette(c("white", "white", "red"))(50), fontsize_row = 7)

