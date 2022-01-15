library(Seurat)
library(stringr)
library(ggplot2)

# indicate dir where files are
home_dir = "."
setwd(home_dir)

############################################################
## Set up file structure and load data
############################################################

# read HTO
batch_name = c("16TH-8Hipp-11Th-11St","4CP-4Hipp-16CP-8Th","4St-8CP-11CP-16Hipp","8St-4-Th-11Hipp-16St")
hashtag = sapply(batch_name,function(x) NULL)
singlet = sapply(batch_name,function(x) NULL)

############################################################
## Run pipeline
############################################################

for (batch in batch_name){
  
  hto <- Read10X(paste("../data/",batch,"/hto/",sep=""), gene.column=1)
  hto <- hto[1:(nrow(hto)-1),]
  rownames(hto) = str_split_fixed(rownames(hto), "-", 2)[,1]
  
  # read expression
  gene <- Read10X(paste("../data/",batch,"/expression/",sep=""))
  colnames(gene) = str_split_fixed(colnames(gene), "-", 2)[,1]
  
  # Select cell barcodes detected by both RNA and HTO In the example datasets we have already
  # filtered the cells for you, but perform this step for clarity.
  joint.bcs <- intersect(colnames(gene), colnames(hto))
  
  # Subset RNA and HTO counts by joint cell barcodes
  gene <- gene[, joint.bcs]
  hto <- as.matrix(hto[, joint.bcs])
  
  # Confirm that the HTO have the correct names
  rownames(hto)
  
  ############################################################
  ## Setup Seurat object and add in the HTO data
  ############################################################
  
  # Setup Seurat object
  hashtag[[batch]] <- CreateSeuratObject(counts = gene, project = batch)
  
  # Add batch information
  hashtag[[batch]]$batch <- batch
  
  ############################################################
  ## Adding HTO as an independent assay
  ############################################################
  
  # Add HTO data as a new assay independent from RNA
  hashtag[[batch]][["HTO"]] <- CreateAssayObject(counts = hto)

  ############################################################
  ## Normalization
  ############################################################

  # Normalize RNA data with log normalization
  hashtag[[batch]] <- NormalizeData(hashtag[[batch]], normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Find variable features
  hashtag[[batch]] <- FindVariableFeatures(hashtag[[batch]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  # Identify the 20 most highly variable genes and plot variable features
  pdf(paste("Top20_variable_genes_",batch,".pdf",sep=""), width = 12, height = 8)
  p <- LabelPoints(plot = VariableFeaturePlot(hashtag[[batch]]), points = head(VariableFeatures(hashtag[[batch]]), 20), repel = TRUE)
  print(p)
  dev.off()
  
  # Scale all genes
  hashtag[[batch]] <- ScaleData(hashtag[[batch]], features = rownames(hashtag[[batch]]))
  
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  hashtag[[batch]] <- NormalizeData(hashtag[[batch]], assay = "HTO", normalization.method = "CLR")
  
  ############################################################
  ## Demultiplex cells based on HTO enrichment
  ############################################################
  
  # If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
  # clustering function for large applications You can also play with additional parameters (see
  # documentation for HTODemux()) to adjust the threshold for classification Here we are using the
  # default settings
  
  hashtag[[batch]] = HTODemux(hashtag[[batch]], assay = "HTO", positive.quantile = 0.99)
  
  ############################################################
  ## Visualize demultiplexing results
  ############################################################
  
  # Global classification results
  table(hashtag[[batch]]$HTO_classification.global)
  
  # HTO classification results
  table(hashtag[[batch]]$HTO_classification)
  table(hashtag[[batch]]$hash.ID)
  
  ## Visualize enrichment for selected HTOs with ridge plots
  # Group cells based on the max HTO signal
  Idents(hashtag[[batch]]) <- "HTO_maxID"
  pdf(paste("Ridgeplot_",batch,".pdf",sep=""))
  p=RidgePlot(hashtag[[batch]], assay = "HTO", features = rownames(hashtag[[batch]][["HTO"]]), ncol = 2)
  print(p)
  dev.off()
  
  ## Compare number of UMIs for singlets, doublets and negative cells
  Idents(hashtag[[batch]]) <- "HTO_classification.global"
  pdf(paste("Vlnplot_",batch,".pdf",sep=""))
  p=VlnPlot(hashtag[[batch]], features = "nCount_RNA", pt.size = 0.1, log = TRUE)
  print(p)
  dev.off()

  ############################################################
  ## Perform QC
  ############################################################
  
  # Estimate MT %
  hashtag[[batch]][["percent.mt"]] <- PercentageFeatureSet(hashtag[[batch]], pattern = "^MT-")
  
  # Print basic QC
  pdf(paste("QC_",batch,".pdf",sep=""), width = 15, height = 10)
  p=VlnPlot(hashtag[[batch]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p)
  plot1 <- FeatureScatter(hashtag[[batch]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(hashtag[[batch]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p=CombinePlots(plots = list(plot1, plot2))
  print(p)
  dev.off()
  
  # Remove samples with poor QC
  hashtag[[batch]] <- subset(hashtag[[batch]], subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 5)

  ############################################################
  ## Visualite doublet/singlet
  ############################################################
  
  ## Generate t-SNE plot
  # First, we will remove negative cells from the object
  hashtag.subset <- subset(hashtag[[batch]], idents = "Negative", invert = TRUE)
  
  # Calculate a distance matrix using HTO
  hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = hashtag.subset, assay = "HTO"))))
  
  # Calculate tSNE embeddings with a distance matrix
  hashtag.subset <- RunTSNE(hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
  pdf(paste("tSNE_singlet_doublet_",batch,".pdf",sep=""))
  p=DimPlot(hashtag.subset)
  print(p)
  dev.off()
  
  # To increase the efficiency of plotting, you can subsample cells using the num.cells argument
  png(paste("HTOHeatmap_",batch,".png",sep=""))
  p=HTOHeatmap(hashtag[[batch]], assay = "HTO")
  print(p)
  dev.off()
  
  ############################################################
  ## Visualite singlet
  ############################################################
  
  ## Cluster and visualize cells using the usual scRNA-seq workflow, and examine for the potential presence of batch effects.
  
  # Extract the singlets
  singlet[[batch]] <- subset(hashtag[[batch]], idents = "Singlet")
  
}

batch = "merged"

merged <- merge(singlet[[1]], y = c(singlet[[2]], singlet[[3]], singlet[[4]]), 
                  add.cell.ids = batch_name, project = "singlet", merge.data = TRUE)

save(merged,file="human_snSeq.Rdata")
