library(Seurat)
library(ggplot2)
library(reticulate)
library(readxl)
use_python("~/.conda/envs/jupyter/bin/python")

# Load data
l3 <- readRDS("larvaOL.integrated150.rds")

# Quick glimpse
DimPlot(l3, label = TRUE) + NoLegend()
DefaultAssay(l3) <- "RNA"

# Use NE markers to locate neuroepithelial cells
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3714590/
#Idents(l3) <- "old.ident"

VlnPlot(l3, c("shg", "Brd", "Tom", "dpn"), assay = "RNA",
        ncol = 1, pt.size = 0)

# Average expression of the marker genes
avg <- AverageExpression(l3, c("shg", "Tom", "Brd"),
                         assay = "RNA")

# Setting threshold for Tom and Brd
tom_thres <- quantile(FetchData(l3, "rna_Tom")[, 1], 0.95)
brd_thres <- quantile(FetchData(l3, "rna_Brd")[, 1], 0.95)

# Identifying candidate NE clusters that are Tom+/Brd+
tom_cluster <- colnames(avg$RNA)[avg$RNA["Tom", ] > tom_thres]
brd_cluster <- colnames(avg$RNA)[avg$RNA["Brd", ] > brd_thres]

ne_c<- intersect(brd_cluster, tom_cluster)

# Examining DE-Cad expression in the candidate NE clusters versus the
# whole dataset
dif <- FindMarkers(l3, ident.1 = ne_c)

# Based on marker expression
# Cluster 9 are NBs
# Cluster 1, 57, 10, 19, 74, 37 could be GMCs
nes <- subset(l3, idents = ne_c)
DefaultAssay(nes) <- "integrated"

nes <- FindVariableFeatures(nes, assay = "RNA",
                            selection.method = "vst",
                            nfeatures = 2000,
                            verbose = FALSE)

nes <- RunPCA(nes, npcs = 250)
ElbowPlot(nes, ndims = 250)
nes <- RunUMAP(nes, dims = 1:100)
nes <- FindNeighbors(nes, dims = 1:100)
nes <- FindClusters(nes, resolution = 0.6)

DimPlot(nes, label = TRUE)

FeaturePlot(nes, c("shg", "Brd", "Tom", "dpn"),
            min.cutoff = "q5", max.cutoff = "q95")



# All possible NEs
# Contains a lot of lamina NE (labeled by Tll, Dac, and Eya)
VlnPlot(nes, c("hth", "tll"), assay = "RNA", pt.size = 0)

# Subsetting OPC NEs by removing Tll+ cells
opcc <- c(0, 2, 4, 6, 8, 9, 10, 11, 12)
opc <- subset(nes, idents = opcc)
opc <- FindVariableFeatures(opc, assay = "RNA",
                            selection.method = "vst", nfeatures = 2000,
                            verbose = FALSE)
opc <- RunPCA(opc, npcs = 250)
ElbowPlot(opc, ndims = 250)

opc <- RunUMAP(opc, dims = 1:30)
opc <- FindNeighbors(opc, dims = 1:30)
opc <- FindClusters(opc, resolution = 1)
#saveRDS(opc, "dataset/opc.rds")

# Clusters with low Hth expression was examined
# and discarded based on the expression of mira, ase, and pros
# This might contains differentiating cells (GMCs?)

opcne <- subset(opc, idents = c(2, 3, 8, 7, 9, 11))
opcne <- FindVariableFeatures(opcne, , assay = "RNA",
                              selection.method = "vst",
                              nfeatures = 2000, verbose = FALSE)
opcne <- RunPCA(opcne, npcs = 250)
ElbowPlot(opcne, ndims = 30)
opcne <- RunUMAP(opcne, dims = 1:20)
opcne <- FindNeighbors(opcne, dims = 1:20)
opcne <- FindClusters(opcne, resolution = 2)
DimPlot(opcne, label = TRUE) + NoLegend()

# This could still contain some lamina precursors
# By the expression os, cv-c, and tll

mne <- subset(opcne, idents = c(0, 2, 3, 4, 6, 7, 9, 10))
mne <- FindVariableFeatures(mne, , assay = "RNA",
                            selection.method = "vst", nfeatures = 2000,
                            verbose = FALSE)
mne <- RunPCA(mne, npcs = 250)
ElbowPlot(mne, ndims = 50)
mne <- RunUMAP(mne, dims = 1:20)
mne <- FindNeighbors(mne, dims = 1:20)
mne <- FindClusters(mne, resolution = 2)
DimPlot(mne, label = TRUE) + NoLegend()
saveRDS(mne, "medulla_ne.rds")

# Supervised PC selection
# Hoping some PCs are more informative than others, 
# I'll try focusing on PCs that contain known spatial factors to see if cells are better separated.
# Get feature embeddings
pc_loading <- mne@reductions$pca@feature.loadings

# Select top and bottom 10 features for each PC
top_features <- apply(pc_loading, 2, function(x) {
  order_by_weight <- order(x)
  results <- row.names(pc_loading)[c(head(order_by_weight, 20), tail(order_by_weight, 20))]
  return(results)
})
head(top_features)

# List PCs with spatial factors in the top 20

stf <- c("dpp", "Vsx1", "Rx", "Optix", "wg", "brk", "cg", "bi", "ci")
stf_pc <- apply(top_features, 2, function(x) {
  stf_present <- length(intersect(x, stf)) > 0
  return(stf_present)
})

# Keep top 30 PCs that contains the spatial factors in the top and bottom 10 features for each PC
pc_use <- seq(length(stf_pc))[stf_pc]
length(pc_use)

# Using spatial PCs to cluster (aiming for over-clustering for manual merge)
mne <- FindNeighbors(mne, reduction = "pca", dims = pc_use, graph.name = "stf_snn")
mne <- RunUMAP(mne, dims = pc_use)
mne <- FindClusters(mne, resolution = 1.5, graph.name = "stf_snn")

DimPlot(mne, reduction = "umap", label = TRUE) + NoLegend()

FeaturePlot(mne, features = c("dpp", "Vsx1", "Rx", "Optix", "wg"), 
            min.cutoff = "q5", max.cutoff = "q95", 
            reduction = "umap")

# Check the size of each clusters 
# to make sure the size is large enough for a classifier to run upon
table(mne@active.ident)

# Build tree and get nodes of low confidence
mne <- BuildClusterTree(mne, assay = "integrated", verbose = TRUE)
PlotClusterTree(mne)

oob_node <- AssessAllNodes(mne, assay = "integrated")
oob_node$low_confidence

# Deal with the terminal low-confidence nodes
merge_temp <- IdentRename(mne@active.ident, old.ident = c(3, 8), new.ident = 3)
mne@meta.data$merged <- merge_temp
Idents(mne) <- mne@meta.data$merged

# Re-build tree and merge again
mne <- BuildClusterTree(mne, assay = "integrated", verbose = TRUE)
oob_node <- AssessAllNodes(mne, assay = "integrated")
oob_node$low_confidence

# Plot semi-supervised merge
DimPlot(mne, label = TRUE, reduction = "umap") + NoLegend()

# Supervised merging with prior knowledge
VlnPlot(mne, features = c("dpp", "Vsx1", "Rx", "Optix", "wg"),
        assay = "RNA", ncol = 2, pt.size = 0.1)


merge_temp <- IdentRename(mne@active.ident, old.ident = 0, new.ident = "Wg")
merge_temp <- IdentRename(merge_temp, old.ident = c(9, 10), new.ident = "Dpp")
merge_temp <- IdentRename(merge_temp, old.ident = 7, new.ident = "Rx")
merge_temp <- IdentRename(merge_temp, old.ident = c(1, 4, 6), new.ident = "Vsx")
merge_temp <- IdentRename(merge_temp, old.ident = c(2, 3, 5), 
                          new.ident = "Optix")
mne@meta.data$annot <- merge_temp
Idents(mne) <- mne@meta.data$annot

VlnPlot(mne, features = c("dpp", "Vsx1", "Rx", "Optix", "wg"), assay = "RNA", ncol = 2, pt.size = 0.1)

# Save temporay result to conserve computation effort
saveRDS(mne, file = "medulla_ne_reanalysis.rds")
