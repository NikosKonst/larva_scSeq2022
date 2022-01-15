library(Seurat)

OL.combined <- readRDS("OL.combined_published.annotations.rds")

### Figure 4 ###

png("classic_Hth.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_hth", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_Ey.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_ey", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_Erm.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_erm", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_D.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_D", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_Scro.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_scro", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_Hbn.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_hbn", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_Opa.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_opa", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_Oaz.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_Oaz", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_B-H1.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_B-H1", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_Slp1.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_slp1", pt.size = 1,  reduction = "umap")
dev.off()

png("classic_tll.png",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_tll", pt.size = 1,  reduction = "umap")
dev.off()

### Figure 5 ###
lab<-read.table("colors.txt", comment.char="")
lab2<-unlist(lab)
names(lab2)<-levels(OL.combined)

png("dimplot_colors.tiff", width = 10000, height=5000, res=200)
DimPlot(OL.combined, label = T, pt.size = 1, cols = lab2, na.value = "red")
dev.off()

### Figure S9 ###
OL.combined<-BuildClusterTree(OL.combined)

lab3<-read.table("colors2.txt", comment.char="")
lab4<-unlist(lab3)
names(lab4)<-levels(OL.combined)

png("Circ_Cluster_tree_all_color2.tiff", width = 10000, height=10000, res=200)
PlotClusterTree(OL.combined, type  = "fan", tip.color = lab4, cex  = 2, font = 2)
dev.off()

