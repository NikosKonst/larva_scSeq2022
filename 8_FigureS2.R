library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(monocle)
library(reshape)

source("LayeredFeaturePlot_white_background_pt.brightness.R")

OL.combined<-readRDS(file="OL.combined_published.annotations.rds")


#### Figure S2A #####

tiff("Figure S2_dpn.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_dpn", reduction = "umap")
dev.off()

tiff("Figure S2_shg.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_shg", reduction = "umap")
dev.off()

tiff("Figure S2_ase.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_ase", reduction = "umap")
dev.off()

tiff("Figure S2_elav.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_elav", reduction = "umap")
dev.off()

tiff("Figure S2_repo.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_repo", reduction = "umap")
dev.off()


#### Figure S2B ####

tiff("Figure S2_gcm.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_gcm", reduction = "umap", pt.size = 0.4)
dev.off()

tiff("Figure S2_eya.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_eya", reduction = "umap", pt.size = 0.4)
dev.off()

tiff("Figure S2_sim.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_sim", reduction = "umap", pt.size = 0.4)
dev.off()

tiff("Figure S2_tll.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_tll", reduction = "umap", pt.size = 0.4)
dev.off()

tiff("Figure S2_dac.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_dac", reduction = "umap", pt.size = 0.4)
dev.off()

tiff("Figure S2_acj6.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_acj6", reduction = "umap", pt.size = 0.4)
dev.off()

### Figure S2C ###

tiff("Figure S2_Optix_Vsx_Rx.tiff",width = 6000, height=2000, res=300)
LayeredFeaturePlot(OL.combined, features = c("rna_Vsx1","rna_Optix","rna_Rx"), reduction = "umap", pt.size = 0.4)
dev.off()

tiff("Figure S2_Vsx.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_Vsx1", reduction = "umap", pt.size = 0.4)
dev.off()

tiff("Figure S2_Optix.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_Optix", reduction = "umap", pt.size = 0.4)
dev.off()

tiff("Figure S2_Rx.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "rna_Rx", reduction = "umap", pt.size = 0.4)
dev.off()

