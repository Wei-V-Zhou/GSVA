
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

# Choose the Cd19 and Stromal groups for annotation colorbar
{
  # Choose Cd19 cluster
  marker.group <- "CD19"
  marker.exprs <- exprs_D4[rownames(exprs_D4)==marker.group, ]
  marker.data <- rbind(marker.exprs, t(ann.cellname_D4))
  colnames(marker.data) <- marker.data[nrow(marker.data), ]
  expr.D4 <- as.matrix(marker.data[ , colnames(marker.data)=="BoneMet_D4"])
  {
    colnames(expr.D4) <- expr.D4[(nrow(expr.D4)-1), ]
    expr.D4.1 <- expr.D4[c(-nrow(expr.D4),-(nrow(expr.D4)-1)), ]
  }
  Cd19.PT_D4 <- t(expr.D4.1)
  row.names(Cd19.PT_D4) <- marker.group
  anno <- ann.cellname_D4[colnames(Cd19.PT_D4), ]
  anno$CellType <- "CD19"
  anno_Cd19 <- colnames(Cd19.PT_D4)[Cd19.PT_D4 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Cd19, 2] = "Others"
  annotation_col_Cd19 = data.frame(anno_df$CellType)
  row.names(annotation_col_Cd19) <- anno[ , 1]
  colnames(annotation_col_Cd19)<-"CD19"
}
{
  # Choose Jag1 cluster
  rm(marker.group, marker.exprs, marker.data, expr.D4, expr.D4.1 )
  marker.group <- "JAG1"
  marker.exprs <- exprs_D4[rownames(exprs_D4)==marker.group, ]
  marker.data <- rbind(marker.exprs, t(ann.cellname_D4))
  colnames(marker.data) <- marker.data[nrow(marker.data), ]
  expr.D4 <- as.matrix(marker.data[ , colnames(marker.data)=="BoneMet_D4"])
  {
    colnames(expr.D4) <- expr.D4[(nrow(expr.D4)-1), ]
    expr.D4.1 <- expr.D4[c(-nrow(expr.D4),-(nrow(expr.D4)-1)), ]
  }
  Jag1.PT_D4 <- t(expr.D4.1)
  row.names(Jag1.PT_D4) <- marker.group
  anno <- ann.cellname_D4[colnames(Jag1.PT_D4), ]
  anno$CellType <- "JAG1"
  anno_Jag1 <- colnames(Jag1.PT_D4)[Jag1.PT_D4 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Jag1, 2] = "Others"
  annotation_col_Jag1 = data.frame(anno_df$CellType)
  row.names(annotation_col_Jag1) <- anno[ , 1]
  colnames(annotation_col_Jag1)<-"JAG1"
}
{
  # Choose Ly6g cluster
  rm(marker.group, marker.exprs, marker.data, expr.D4, expr.D4.1 )
  marker.group <- "LY6G"
  marker.exprs <- exprs_D4[rownames(exprs_D4)==marker.group, ]
  marker.data <- rbind(marker.exprs, t(ann.cellname_D4))
  colnames(marker.data) <- marker.data[nrow(marker.data), ]
  expr.D4 <- as.matrix(marker.data[ , colnames(marker.data)=="BoneMet_D4"])
  {
    colnames(expr.D4) <- expr.D4[(nrow(expr.D4)-1), ]
    expr.D4.1 <- expr.D4[c(-nrow(expr.D4),-(nrow(expr.D4)-1)), ]
  }
  Ly6g.PT_D4 <- t(expr.D4.1)
  row.names(Ly6g.PT_D4) <- marker.group
  anno <- ann.cellname_D4[colnames(Ly6g.PT_D4), ]
  anno$CellType <- "LY6G"
  anno_Ly6g <- colnames(Ly6g.PT_D4)[Ly6g.PT_D4 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Ly6g, 2] = "Others"
  annotation_col_Ly6g = data.frame(anno_df$CellType)
  row.names(annotation_col_Ly6g) <- anno[ , 1]
  colnames(annotation_col_Ly6g)<-"LY6G"
}
{
  # Choose Cd68 cluster
  rm(marker.group, marker.exprs, marker.data, expr.D4, expr.D4.1 )
  marker.group <- "CD68"
  marker.exprs <- exprs_D4[rownames(exprs_D4)==marker.group, ]
  marker.data <- rbind(marker.exprs, t(ann.cellname_D4))
  colnames(marker.data) <- marker.data[nrow(marker.data), ]
  expr.D4 <- as.matrix(marker.data[ , colnames(marker.data)=="BoneMet_D4"])
  {
    colnames(expr.D4) <- expr.D4[(nrow(expr.D4)-1), ]
    expr.D4.1 <- expr.D4[c(-nrow(expr.D4),-(nrow(expr.D4)-1)), ]
  }
  Cd68.PT_D4 <- t(expr.D4.1)
  row.names(Cd68.PT_D4) <- marker.group
  anno <- ann.cellname_D4[colnames(Cd68.PT_D4), ]
  anno$CellType <- "CD68"
  anno_Cd68 <- colnames(Cd68.PT_D4)[Cd68.PT_D4 < 0]
  anno_df <- data.frame(anno)
  anno_df[anno_Cd68, 2] = "Others"
  annotation_col_Cd68 = data.frame(anno_df$CellType)
  row.names(annotation_col_Cd68) <- anno[ , 1]
  colnames(annotation_col_Cd68)<-"CD68"
}

annotat_col = data.frame(c(annotation_col_Cd19, annotation_col_Jag1, annotation_col_Ly6g, annotation_col_Cd68))
row.names(annotat_col) = anno_df[ , 1]
annot_color = list(CD19 = c(CD19 = "red", Others = "white"), JAG1 = c(JAG1 = "green", Others = "white"),
                   LY6G = c(LY6G = "blue", Others = "white"), CD68 = c(CD68 = "yellow", Others = "white"))

# Run GSVA on the Breast Cancer data using Hallmark Gene Sets
datasets.hm <- gsva(as.matrix(exprs_D4), hm.sets, method="gsva")
rownames(datasets.hm) <- gsub("HALLMARK_","",rownames(datasets.hm))
oncogene <- c("TGF_BETA_SIGNALING", "G2M_CHECKPOINT", "OXIDATIVE_PHOSPHORYLATION", "INFLAMMATORY_RESPONSE",
              "EPITHELIAL_MESENCHYMAL_TRANSITION", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "E2F_TARGETS", "UNFOLDED_PROTEIN_RESPONSE",
              "COMPLEMENT", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "IL6_JAK_STAT3_SIGNALING")
dat.hm <- datasets.hm[oncogene, ]
# hot_genes <- c("KLF4")
# Add_genes <- exprs_D4[hot_genes, ]
# DAT.HM <- rbind(datasets.hm, Add_genes)

# Generate GSVA analysis
pheatmap(datasets.hm, treeheight_col = 30, treeheight_row = 0, show_colnames = F, clustering_method="ward.D2", scale = "column",
         color = colorRampPalette(c("blue", "white", "red"))(50), fontsize_row = 7, annotation_col = annotat_col, annotation_colors = annot_color)

# heatmap.2(datasets.hm, dendrogram = "col", lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 ),
# col=c(rep("blue",20), dChip.colors(100), rep("red",20)), trace="none", keysize = 1,
# scale=("row"), dist=cor.dist, hclust=hclust.ward.d)$colInd


#============================#
#       Musician: Resonance  #
#           Date: 2019/10/05 #
# Revised author: Resonance  #
#           Time: 2019/10/26 #
#============================#