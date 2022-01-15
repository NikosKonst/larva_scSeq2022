library(Seurat)
library(Matrix)

### cortex ###

tab<-read.table("GSE118953_raw_count.tsv", h=T)

row.names(tab)<-tab[,1]
tab<-tab[,-1]

cds<-CreateSeuratObject(counts = tab,project="MouseCortex")
VlnPlot(cds, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)


cds <- NormalizeData(cds)
cds <- FindVariableFeatures(cds, selection.method = "vst", nfeatures = 2000)
cds <- ScaleData(cds)
cds <- RunPCA(cds, features = VariableFeatures(object = cds))

cds <- JackStraw(cds, num.replicate = 100)
cds <- ScoreJackStraw(cds, dims = 1:20)
JackStrawPlot(cds, dims = 1:20)

ElbowPlot(cds, ndims = 25)


cds <- FindNeighbors(cds, dims = 1:25)
cds <- FindClusters(cds, resolution = 0.5)

cds <- RunUMAP(cds, dims = 1:25)


# radial glia is clusters 2 and 3 #
RG<-subset(cds,idents=c(2,3))
RG <- RunUMAP(RG, dims = 1:25)

Idents(RG)<-RG@meta.data$orig.ident
RG_sub<-subset(RG,idents=c("E12.1H", "E13.1H", "E14.1H", "E15.1H"))


### Hth, opa, erm, Hbn, ey, slp1, Scro, D, BarHI, Tll ###
### Meis1/Meis2, Zic1-5, Fezf1, Fezf2, Arx, Pax6, Foxg1, Nkx2-1 (Ttf1), Nkx2-4, Sox12, Sox14, Sox21, Barhl1, Barhl2, Nr2e1 ###

tiff("Figure S11A_Meis.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Meis2"))
dev.off()

tiff("Figure S11A_Zic5.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Zic5"))
dev.off()

tiff("Figure S11A_Fezf2.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Fezf2"))
dev.off()

tiff("Figure S11A_Arx.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Arx"))
dev.off()

tiff("Figure S11A_Pax6.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Pax6"))
dev.off()

tiff("Figure S11A_Foxg1.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Foxg1"))
dev.off()

tiff("Figure S11A_Ttf1.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Ttf1"))
dev.off()

tiff("Figure S11A_Sox12.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Sox12"))
dev.off()

tiff("Figure S11A_Nr2e1.tiff",width = 1500, height=1500, res = 300)
VlnPlot(RG_sub, features = c("Nr2e1"))
dev.off()

### retina ###

genes<-read.table("data2/genes.tsv")
cells<-read.table("data2/barcodes.tsv")
matrix<-readMM("data2/matrix.mtx")

rownames(matrix)<-cells[,1]
colnames(matrix)<-genes[,2]
retina<-CreateSeuratObject(counts=t(matrix))

retina<-NormalizeData(retina)
retina <- FindVariableFeatures(retina, selection.method = "vst", nfeatures = 2000)
retina <- ScaleData(retina)
retina <- RunPCA(retina, features = VariableFeatures(object = retina))
retina <- FindNeighbors(retina, dims = 1:50)
retina <- FindClusters(retina, resolution = 0.5)
retina <- RunUMAP(retina, dims = 1:50)

DimPlot(retina, reduction = "umap", raster = F)
FeaturePlot(retina, features = c("Pax6","Ccnd1", "Cdk4"))

retina<-AddMetaData(retina, metadata = cells$age, col.name = "age")
retina<-AddMetaData(retina, metadata = cells$umap2_CellType, col.name = "celltype")
Idents(retina)<-retina@meta.data$celltype
#Idents(retina)<-retina@meta.data$age

retina_RPCs<-subset(retina, idents = c("Early RPCs", "Late RPCs"))
retina_RPCs <- FindVariableFeatures(retina_RPCs, selection.method = "vst", nfeatures = 2000)
retina_RPCs <- ScaleData(retina_RPCs)
retina_RPCs <- RunPCA(retina_RPCs, features = VariableFeatures(object = retina))
retina_RPCs <- FindNeighbors(retina_RPCs, dims = 1:50)
retina_RPCs <- FindClusters(retina_RPCs, resolution = 0.5)
retina_RPCs <- RunUMAP(retina_RPCs, dims = 1:50)
DimPlot(retina_RPCs, reduction = "umap", raster = F)

FeaturePlot(retina_RPCs, features = c("Ikzf1", "Pou2f2", "Casz1"))

Idents(retina_RPCs)<-retina_RPCs@meta.data$age

retina_RPCs_sub<-subset(retina_RPCs, idents = c("E11", "E12","E14", "E16","E18", "P0"))


tiff("Figure S11C_Ikzf1.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Ikzf1"))
dev.off()

tiff("Figure S11C_Pou2f2.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Pou2f2"))
dev.off()

tiff("Figure S11C_Casz1.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Casz1"))
dev.off()

tiff("Figure S11D_Meis.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Meis2"))
dev.off()

tiff("Figure S11D_Zic5.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Zic5"))
dev.off()

tiff("Figure S11D_Fezf2.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Fezf2"))
dev.off()

tiff("Figure S11D_Arx.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Arx"))
dev.off()

tiff("Figure S11D_Pax6.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Pax6"))
dev.off()

tiff("Figure S11D_Foxg1.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Foxg1"))
dev.off()

tiff("Figure S11D_Ttf1.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Ttf1"))
dev.off()

tiff("Figure S11D_Sox12.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Sox12"))
dev.off()

tiff("Figure S11D_Nr2e1.tiff",width = 1500, height=1500, res = 300)
VlnPlot(retina_RPCs_sub, features = c("Nr2e1"))
dev.off()

