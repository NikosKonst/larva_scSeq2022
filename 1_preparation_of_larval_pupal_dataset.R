library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

Lib1.data<-Read10X(data.dir="Lib1/Lib_10X_larva_2018_12_20_Lib1/outs/filtered_feature_bc_matrix/")
Lib1 <-CreateSeuratObject(counts =Lib1.data, min.cells = 3, min.features = 200,project="Lib1")
Lib1 <- NormalizeData(object = Lib1, verbose = FALSE)
Lib1 <- FindVariableFeatures(object = Lib1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib2.data<-Read10X(data.dir="Lib2/Lib_10X_larva_2018_12_20_Lib2/outs/filtered_feature_bc_matrix/")
Lib2 <-CreateSeuratObject(counts =Lib2.data, min.cells = 3, min.features = 200,project="Lib2")
Lib2 <- NormalizeData(object = Lib2, verbose = FALSE)
Lib2 <- FindVariableFeatures(object = Lib2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib3.data<-Read10X(data.dir="Lib3/Lib_10X_larva_2018_12_20_Lib3/outs/filtered_feature_bc_matrix/")
Lib3 <-CreateSeuratObject(counts =Lib3.data, min.cells = 3, min.features = 200,project="Lib3")
Lib3 <- NormalizeData(object = Lib3, verbose = FALSE)
Lib3 <- FindVariableFeatures(object = Lib3, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib4.data<-Read10X(data.dir="Lib4/Lib_10X_larva_2019_04_12_Lib4/outs/filtered_feature_bc_matrix")
Lib4 <-CreateSeuratObject(counts =Lib4.data, min.cells = 3, min.features = 200,project="Lib4")
Lib4 <- NormalizeData(object = Lib4, verbose = FALSE)
Lib4 <- FindVariableFeatures(object = Lib4, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib5.data<-Read10X(data.dir="Lib5/Lib_10X_larva_2019_04_12_Lib5/outs/filtered_feature_bc_matrix")
Lib5 <-CreateSeuratObject(counts =Lib5.data, min.cells = 3, min.features = 200,project="Lib5")
Lib5 <- NormalizeData(object = Lib5, verbose = FALSE)
Lib5 <- FindVariableFeatures(object = Lib5, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib6.data<-Read10X(data.dir="Lib6/Lib_10X_larva_2019_04_12_Lib6/outs/filtered_feature_bc_matrix")
Lib6 <-CreateSeuratObject(counts =Lib6.data, min.cells = 3, min.features = 200,project="Lib6")
Lib6 <- NormalizeData(object = Lib6, verbose = FALSE)
Lib6 <- FindVariableFeatures(object = Lib6, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib7.data<-Read10X(data.dir="Lib7/Lib_10X_larva_2019_04_12_Lib7/outs/filtered_feature_bc_matrix")
Lib7 <-CreateSeuratObject(counts =Lib7.data, min.cells = 3, min.features = 200,project="Lib7")
Lib7 <- NormalizeData(object = Lib7, verbose = FALSE)
Lib7 <- FindVariableFeatures(object = Lib7, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib8.data<-Read10X(data.dir="Lib8/Lib_10X_larva_2019_04_12_Lib8/outs/filtered_feature_bc_matrix")
Lib8 <-CreateSeuratObject(counts =Lib8.data, min.cells = 3, min.features = 200,project="Lib8")
Lib8 <- NormalizeData(object = Lib8, verbose = FALSE)
Lib8 <- FindVariableFeatures(object = Lib8, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib9.data<-Read10X(data.dir="Lib9/Lib_10X_larva_2019_04_12_Lib9/outs/filtered_feature_bc_matrix")
Lib9 <-CreateSeuratObject(counts =Lib9.data, min.cells = 3, min.features = 200,project="Lib9")
Lib9 <- NormalizeData(object = Lib9, verbose = FALSE)
Lib9 <- FindVariableFeatures(object = Lib9, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

Lib10.data<-Read10X(data.dir="Lib10/Lib_10X_larva_2019_04_12_Lib10/outs/filtered_feature_bc_matrix")
Lib10 <-CreateSeuratObject(counts =Lib10.data, min.cells = 3, min.features = 200,project="Lib10")
Lib10 <- NormalizeData(object = Lib10, verbose = FALSE)
Lib10 <- FindVariableFeatures(object = Lib10, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

reference.list <- c(Lib1, Lib2, Lib3,Lib4,Lib5,Lib6, Lib7, Lib8,Lib9,Lib10)
larvaOL.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:150)
larvaOL.integrated150 <- IntegrateData(anchorset = larvaOL.anchors, dims = 1:150)
larvaOL.integrated150 <- ScaleData(object = larvaOL.integrated150, verbose = FALSE)
larvaOL.integrated150 <- RunPCA(object = larvaOL.integrated150, npcs = 150, verbose = FALSE)

A=MixingMetric(object = larvaOL.integrated150, grouping.var = "orig.ident",reduction = "pca",dims=1:150)
#save(A, file="MixMet150.Rdata")
B=LocalStruct(object=larvaOL.integrated150, grouping.var = "orig.ident")
#save(B,file="LocStruc150.Rdata")


larvaOL.integrated150<-RunTSNE(object = larvaOL.integrated150, dims = 1:150)
larvaOL.integrated150 <- RunUMAP(object = larvaOL.integrated150, dims = 1:150)

larvaOL.integrated150<-FindNeighbors(object = larvaOL.integrated150, dims = 1:150)
larvaOL.integrated150<-FindClusters(object = larvaOL.integrated150, resolution = 2)

saveRDS(larvaOL.integrated150, file="larvaOL.integrated150.rds")


#### Merge with the P15 dataset ###

P15_new<-readRDS("P15_10Xv3Adjusted.rds")
P15<-UpdateSeuratObject(P15_new)

OL.combined <- merge(larva, y = P15, add.cell.ids = c("L3", "P15"), project = "OL")
DefaultAssay(OL.combined)<-"RNA"
OL.combined <-NormalizeData(OL.combined)
OL.combined <- FindVariableFeatures(object = OL.combined, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
OL.combined <-ScaleData(object = OL.combined, verbose = FALSE)
OL.combined <- RunPCA(object = OL.combined, npcs = 150, verbose = FALSE)
OL.combined<-RunTSNE(object = OL.combined, dims = 1:150)
OL.combined <- RunUMAP(object = OL.combined, dims = 1:150)

saveRDS(OL.combined, file="OL.combined150.rds")

