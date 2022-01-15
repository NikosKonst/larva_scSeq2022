library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggthemes)

larva<-readRDS("larvaOL.integrated150.rds")

Idents(larva)<-larva@meta.data$orig.ident

split<-SplitObject(larva, split.by = "ident")

merged_object1<-merge(split$Lib1, split$Lib2)
merged_object2<-merge(merged_object1, split$Lib3)
merged_object3<-merge(merged_object2, split$Lib4)
merged_object4<-merge(merged_object3, split$Lib5)
merged_object5<-merge(merged_object4, split$Lib6)
merged_object6<-merge(merged_object5, split$Lib7)
merged_object7<-merge(merged_object6, split$Lib8)
merged_object8<-merge(merged_object7, split$Lib9)
merged_object<-merge(merged_object8, split$Lib10)


DefaultAssay(merged_object)<-"RNA"
merged_object <- NormalizeData(object = merged_object, verbose = FALSE)
merged_object<- FindVariableFeatures(merged_object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged_object <- ScaleData(object = merged_object, verbose = FALSE)
merged_object <- RunPCA(object = merged_object, npcs = 150, verbose = FALSE)
merged_object<-RunTSNE(object = merged_object, dims = 1:150)
merged_object <- RunUMAP(object = merged_object, dims = 1:150)
merged_object<-FindNeighbors(object = merged_object, dims = 1:150)
merged_object<-FindClusters(object = merged_object, resolution = 2)

Idents(merged_object)<-merged_object@meta.data$orig.ident

tiff("new Figure S1A.tiff",width = 4500, height=3000, res=300)
DimPlot(merged_object, reduction = "tsne", pt.size = 1.5)
dev.off()

Idents(larva)<-larva@meta.data$orig.ident
tiff("new Figure S1A'.tiff",width = 4500, height=3000, res=300)
DimPlot(larva, reduction = "tsne", pt.size = 1.5)
dev.off()

Idents(merged_object)<-merged_object@meta.data$RNA_snn_res.2
merged_object_clusters <- SplitObject(merged_object, split.by = "ident")

saveRDS(merged_object, "larva_merged.rds")

merged_object<-readRDS("larva_merged.rds")
larva<-readRDS("larvaOL.integrated150.rds")

Idents(merged_object)<-merged_object@meta.data$RNA_snn_res.2
tab<-matrix(nrow=length(table(Idents(merged_object))), ncol=10)

for (i in 0:(length(table(Idents(merged_object)))-1)) {
A=SubsetData(merged_object, ident.use = i)
Idents(A)<-A@meta.data$orig.ident
tab[i,]<-table(factor(Idents(A), levels = c("Lib1", "Lib2", "Lib3", "Lib4", "Lib5", "Lib6", "Lib7", "Lib8", "Lib9", "Lib10")))
print(i)
}

row.names(tab)=c(0:(length(table(Idents(merged_object)))-1))
colnames(tab) = c("Lib1", "Lib2", "Lib3", "Lib4", "Lib5", "Lib6", "Lib7", "Lib8", "Lib9", "Lib10")
Idents(larva)<-larva@meta.data$integrated_snn_res.2
tab2<-matrix(nrow=length(table(Idents(larva))), ncol=10)

for (i in 0:(length(table(Idents(larva)))-1)) {
  A=SubsetData(larva, ident.use = i)
  Idents(A)<-A@meta.data$orig.ident
  tab2[i,]<-table(factor(Idents(A), levels = c("Lib1", "Lib2", "Lib3", "Lib4", "Lib5", "Lib6", "Lib7", "Lib8", "Lib9", "Lib10")))
  print(i)
}
row.names(tab2)=c(0:(length(table(Idents(larva)))-1))
colnames(tab2) = c("Lib1", "Lib2", "Lib3", "Lib4", "Lib5", "Lib6", "Lib7", "Lib8", "Lib9", "Lib10")

head(tab)

counts=as.vector(t(tab))
libraries=rep(c("Lib1", "Lib2", "Lib3", "Lib4", "Lib5", "Lib6", "Lib7", "Lib8", "Lib9", "Lib10"),length(table(Idents(merged_object))) )
identities=rep(0:(length(table(Idents(merged_object)))-1), each = 10)

tab_df<-data.frame(counts, libraries, identities)

plot1<- ggplot(tab_df, aes(fill=libraries, y=counts, x=identities)) + 
  geom_bar(position="fill", stat="identity") + theme(panel.background = element_rect(fill = "white",
                                                                                     colour = "white",
                                                                                     size = 0.5, linetype = "solid"))+ scale_x_continuous(expand = c(0, 0))

ggsave("new Figure S1B.tiff",plot = plot1, units = "cm", width = 54, height = 12)

counts2=as.vector(t(tab2))
libraries2=rep(c("Lib1", "Lib2", "Lib3", "Lib4", "Lib5", "Lib6", "Lib7", "Lib8", "Lib9", "Lib10"),length(table(Idents(larva))) )
identities2=rep(0:(length(table(Idents(larva)))-1), each = 10)

tab2_df<-data.frame(counts2, libraries2, identities2)
colnames(tab2_df)<-c("counts", "libraries", "identities")

plot2 <- ggplot(tab2_df, aes(fill=libraries, y=counts, x=identities)) + 
  geom_bar(position="fill", stat="identity") + theme(panel.background = element_rect(fill = "white",
                                                                                     colour = "white",
                                                                                     size = 0.5, linetype = "solid"))   + scale_x_continuous(expand = c(0, 0))
ggsave("new Figure S1B'.tiff",plot = plot2, units = "cm", width = 54, height = 12)

ggarrange(plot1, plot2)

boxplot(tab_df$counts~tab_df$libraries, col ="green")

boxplot(tab2_df$counts~tab2_df$libraries)

tiff("new Figure S1C.tiff",width = 9000, height=3000, res=300)


tab_all <- rbind(tab_df, tab2_df)

ggplot(data = tab_all, aes(x=libraries, y=counts)) + 
  geom_boxplot(aes(fill=origin)) + theme(panel.background = element_rect(fill = "white",
                                                                         colour = "white",
                                                                         size = 0.5, linetype = "solid"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black")) + scale_y_continuous(expand = c(0, 0)) 
ggsave("new Figure S1C.tiff",units = "cm", width = 54, height = 18)
