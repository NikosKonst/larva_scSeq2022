library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(monocle3)
library(GO.db)
source("functions.r")
library(plotly)

##___Open L3 + P15 combined dataset
larva_dataset <- readRDS("OL.combined_published.annotations.rds")

tiff("FigureS10A.tiff",width = 6000, height=2000, res=300)
FeaturePlot(OL.combined, features = "bsh")
dev.off()

png("Figure S10A'.tiff", width = 10000, height=5000, res=200)
DimPlot(larva_dataset, label = T)
dev.off()

##___remove low-quality cells
larva_dataset[['percent.mito']] <- PercentageFeatureSet(larva_dataset, pattern = "^mt:")
larva_dataset <- subset(larva_dataset, nCount_RNA > 1500 & percent.mito < 10)


##___select clusters corresponding to Mi1 cells (identified via bsh expression)
## pupal cells with a RNA count inferior to 3000 are discarded
#L3
Mi1L3Seurat = subset(larva_dataset,Azalia_dim150_res10 %in% c("L3_77","L3_91","L3_86","L3_47")) ###We do not select cluster L3_194 because even though these cells look like Mi1, they express Wnt4 and cluster appart from the others
colnames(Mi1L3Seurat@assays[["RNA"]]@data) <- substr(colnames(Mi1L3Seurat), 4,25)
colnames(Mi1L3Seurat@assays[["RNA"]]@counts) <- substr(colnames(Mi1L3Seurat), 4,25) ##harmonizing names
Mi1L3Seurat <- AddMetaData(Mi1L3Seurat, "L3", col.name = "age")
#P15
Pupa15 <- readRDS("P15_10Xv3Adjusted.rds")
Pupa15 <- UpdateSeuratObject(Pupa15)
Mi1P15Seurat = subset(Pupa15, NNPredsFixed == 142)
Mi1P15Seurat = subset(Mi1P15Seurat, nCount_RNA >3000)
Mi1P15Seurat <- AddMetaData(Mi1P15Seurat, "P15", col.name = "age")
rm(Pupa15)
#P30
Pupa30 <- readRDS("P30_10Xv3Adjusted.rds")
Pupa30 <- UpdateSeuratObject(Pupa30)
Mi1P30Seurat = subset(Pupa30, NNPredsFinal == 142)
Mi1P30Seurat = subset(Mi1P30Seurat, nCount_RNA >3000)
Mi1P30Seurat <- AddMetaData(Mi1P30Seurat, "P30", col.name = "age")
rm(Pupa30)
#P40
Pupa40 <- readRDS("P40_10Xv3Adjusted.rds")
Pupa40 <- UpdateSeuratObject(Pupa40)
Mi1P40Seurat = subset(Pupa40, NNPredsFinal == 142)
Mi1P40Seurat = subset(Mi1P40Seurat, nCount_RNA >3000)
Mi1P40Seurat <- AddMetaData(Mi1P40Seurat, "P40", col.name = "age")
rm(Pupa40)
#P50
Pupa50 <- readRDS("P50_10Xv3Adjusted.rds")
Pupa50 <- UpdateSeuratObject(Pupa50)
Mi1P50Seurat = subset(Pupa50, NNPredsFinal == 142)
Mi1P50Seurat = subset(Mi1P50Seurat, nCount_RNA >3000)
Mi1P50Seurat <- AddMetaData(Mi1P50Seurat, "P50", col.name = "age")
rm(Pupa50)
#P70
Pupa70 <- readRDS("P70_10Xv3Adjusted.rds")
Pupa70 <- UpdateSeuratObject(Pupa70)
Mi1P70Seurat = subset(Pupa70, NNPredsFixed == 142)
Mi1P70Seurat = subset(Mi1P70Seurat, nCount_RNA >3000)
Mi1P70Seurat <- AddMetaData(Mi1P70Seurat, "P70", col.name = "age")
rm(Pupa70)

mergedSeurat = merge(Mi1L3Seurat, y = c(Mi1P15Seurat, Mi1P30Seurat, Mi1P40Seurat, Mi1P50Seurat, Mi1P70Seurat), project = "merged", add.cell.id = c("L3","P15","P30","P40","P50","P70"))

saveRDS(mergedSeurat, "Mi1_all_ages_merged.rds")

mergedSeurat<-readRDS("Mi1_all_ages_merged.rds")

##___Input these Mi1 cells into monocle
data <- GetAssayData(object = mergedSeurat, slot = "counts",assay = "RNA")
pData <- data.frame(cell_id = colnames(data),
                    row.names = colnames(data),
                    orig.ident = mergedSeurat[["orig.ident"]]$orig.ident,
                    age = mergedSeurat[["age"]]$age,
                    nCount_RNA = mergedSeurat[["nCount_RNA"]]$nCount_RNA)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
cdsMi1 <- new_cell_data_set(data, pData, fData)


set.seed(12)
cdsMi1 <- preprocess_cds(cdsMi1, num_dim = 100)
cdsMi1 <- reduce_dimension(cdsMi1)


plot_cells(cdsMi1[,cdsMi1[["age"]]=="L3"], color_cells_by = "orig.ident")

##___Batch alignment
##we here separate the cells by age to batch correct them by looking only at cells that are the same age

cdsMi1L3 <- cdsMi1[,cdsMi1[["age"]]=="L3"]
cdsMi1P15 <- cdsMi1[,cdsMi1[["age"]]=="P15"]
cdsMi1P30 <- cdsMi1[,cdsMi1[["age"]]=="P30"]
cdsMi1P40 <- cdsMi1[,cdsMi1[["age"]]=="P40"]
cdsMi1P50 <- cdsMi1[,cdsMi1[["age"]]=="P50"]
cdsMi1P70 <- cdsMi1[,cdsMi1[["age"]]=="P70"]

#the next lines are necessary to remove unwanted empty factors that make the code crash
cdsMi1L3@colData@listData[["orig.ident"]] <- factor(as.character(cdsMi1L3@colData@listData[["orig.ident"]]))
cdsMi1P15@colData@listData[["orig.ident"]] <- factor(as.character(cdsMi1P15@colData@listData[["orig.ident"]]))
cdsMi1P30@colData@listData[["orig.ident"]] <- factor(as.character(cdsMi1P30@colData@listData[["orig.ident"]]))
cdsMi1P40@colData@listData[["orig.ident"]] <- factor(as.character(cdsMi1P40@colData@listData[["orig.ident"]]))
cdsMi1P50@colData@listData[["orig.ident"]] <- factor(as.character(cdsMi1P50@colData@listData[["orig.ident"]]))
cdsMi1P70@colData@listData[["orig.ident"]] <- factor(as.character(cdsMi1P70@colData@listData[["orig.ident"]]))


cdsMi1L3 <- align_cds(cdsMi1L3, num_dim = 100, alignment_group = "orig.ident")
cdsMi1L3 <- reduce_dimension(cdsMi1L3)
## Uncomment the next line to see the result
#plot_cells(cdsMi1L3, color_cells_by = "orig.ident")
cdsMi1P15 <- align_cds(cdsMi1P15, num_dim = 100, alignment_group = "orig.ident")
cdsMi1P15 <- reduce_dimension(cdsMi1P15)
cdsMi1P30 <- align_cds(cdsMi1P30, num_dim = 100, alignment_group = "orig.ident")
cdsMi1P30 <- reduce_dimension(cdsMi1P30)
cdsMi1P40 <- align_cds(cdsMi1P40, num_dim = 100, alignment_group = "orig.ident")
cdsMi1P40 <- reduce_dimension(cdsMi1P40)
cdsMi1P50 <- align_cds(cdsMi1P50, num_dim = 100, alignment_group = "orig.ident")
cdsMi1P50 <- reduce_dimension(cdsMi1P50)
cdsMi1P70 <- align_cds(cdsMi1P70, num_dim = 100, alignment_group = "orig.ident")
cdsMi1P70 <- reduce_dimension(cdsMi1P70)
#NOTE : an error can occur here or later if there is not enough cells compared to the number of dimensions you want to use

## reinject the batch corrected sata into the cds object
cdsMi1@int_colData <- rbind(cdsMi1L3@int_colData,cdsMi1P15@int_colData,cdsMi1P30@int_colData,cdsMi1P40@int_colData,cdsMi1P50@int_colData,cdsMi1P70@int_colData)
cdsMi1@colData <-  rbind(cdsMi1L3@colData, cdsMi1P15@colData, cdsMi1P30@colData,cdsMi1P40@colData,cdsMi1P50@colData,cdsMi1P70@colData)
cdsMi1@assays <- cbind(cdsMi1L3@assays, cdsMi1P15@assays, cdsMi1P30@assays,cdsMi1P40@assays,cdsMi1P50@assays,cdsMi1P70@assays)

set.seed(585)
cdsMi1 <- reduce_dimension(cdsMi1)
plot_cells(cdsMi1, color_cells_by = "age")


##___Giving pseudotime to data

#############trying to give a pseudotime to pupal data

#cluster cells
set.seed(15)
cdsMi1 <- cluster_cells(cdsMi1, resolution = 1.5e-2)
plot_cells(cdsMi1, color_cells_by = "partition")
#compute trajectories
cdsMi1 <- learn_graph(cdsMi1,learn_graph_control = list(minimal_branch_len= 20,L1.sigma= 0.2))
plot_cells(cdsMi1)
#order trajectories, the beginning of P15 and L3 trajectories must be pointed at manually. To identify them we look at Ase
plot_cells(cdsMi1, genes = "ase")
cdsMi1 <- order_cells(cdsMi1, root_cells = c("P15_1_CTGCGGAAGCAACGGT","L3_GCAGTTAGTCGCGGTT_7"))
plot_cells(cdsMi1, color_cells_by = "partition")


tiff("Figure 6A.tiff", width = 2000, height=2000, res=300)
plot_cells(cdsMi1, color_cells_by = "age", cell_size = 0.5, show_trajectory_graph = T, label_cell_groups = F, label_leaves = F, label_roots = F, label_groups_by_cluster = F, label_branch_points = F, labels_per_group = F)
dev.off()

tiff("Figure S10B_ase.tiff", width = 2000, height=2000, res=300)
plot_cells(cdsMi1, genes = "ase", cell_size = 0.5, show_trajectory_graph = F, label_cell_groups = F, label_leaves = F, label_roots = F, label_groups_by_cluster = F, label_branch_points = F, labels_per_group = F)
dev.off()

tiff("Figure S10B'_bsh.tiff", width = 2000, height=2000, res=300)
plot_cells(cdsMi1, genes = "bsh", cell_size = 0.5, show_trajectory_graph = F, label_cell_groups = F, label_leaves = F, label_roots = F, label_groups_by_cluster = F, label_branch_points = F, labels_per_group = F)
dev.off()

##we now use Ggamma30A as a reference gene to align pseudotime

# Note: make sure the partitions are correct #
partL3 <- names(partitions(cdsMi1)[partitions(cdsMi1) ==2]) #PL3
partP15 <- names(partitions(cdsMi1)[partitions(cdsMi1) ==1]) #P15
partP30 <- names(partitions(cdsMi1)[partitions(cdsMi1) ==3]) #P30
partP40 <- names(partitions(cdsMi1)[partitions(cdsMi1) ==6]) #P40
partP50 <- names(partitions(cdsMi1)[partitions(cdsMi1) ==5]) #P50
partP70 <- names(partitions(cdsMi1)[partitions(cdsMi1) ==4]) #P70
ClEndP15 <- names(clusters(cdsMi1)[clusters(cdsMi1) ==1]) ## synchronized P15cells

### ! be careful ! clusters and partitions can change over rerunning


plot_cells(cdsMi1, genes = "Ggamma30A")
reference_cells <- colnames(cdsMi1[,partL3])[pseudotime(cdsMi1[,partL3])>7]
plot_genes_in_pseudotime(cdsMi1["Ggamma30A",partL3])+ scale_y_continuous()


l3data = exprs(cdsMi1)["Ggamma30A",reference_cells]
l3timeref= pseudotime(cdsMi1[,reference_cells])
linearmodl3 = lm(formula = exp ~ time, data = data.frame(exp =l3data, time = l3timeref))

P15data = exprs(cdsMi1)["Ggamma30A",partP15]
P15timeref= pseudotime(cdsMi1[,partP15])
linearmodP15 = lm(formula = exp ~ time, data = data.frame(exp =P15data, time = P15timeref))


P15time = pseudotime(cdsMi1[,names(partitions(cdsMi1)[partitions(cdsMi1) ==1])]) # the  following partition should agree with the partition for P15 before
P15timadj = (P15time*linearmodP15$coefficients[2] + linearmodP15$coefficients[1] - linearmodl3$coefficients[1])/linearmodl3$coefficients[2]

#verification of alignment
test = rbind(data.frame(time = l3timeref, expr = l3data, group = 1), data.frame(time = P15timadj, expr = exprs(cdsMi1)["Ggamma30A",partP15], group = 2))
ggplot(data=test, aes(x=time, y=expr, color = group))  + geom_point()
test2 =  rbind(data.frame(time = l3timeref, expr = exprs(cdsMi1)["noe",reference_cells], group = 1), data.frame(time = P15timadj, expr = exprs(cdsMi1)["noe",partP15], group = 2))
ggplot(data=test2, aes(x=time, y=expr, color = group))  + geom_point()


##modification of pseudotime in data
cdsMi1@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][partP15] <- P15timadj
cdsMi1@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][partP30] <-  mean(pseudotime(cdsMi1)[ClEndP15]) + 15/3
cdsMi1@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][partP40] <-  mean(pseudotime(cdsMi1)[ClEndP15]) + 25/3
cdsMi1@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][partP50] <-  mean(pseudotime(cdsMi1)[ClEndP15]) + 35/3
cdsMi1@principal_graph_aux@listData[["UMAP"]][["pseudotime"]][partP70] <-  mean(pseudotime(cdsMi1)[ClEndP15]) + 55/3

plot_genes_in_pseudotime(cdsMi1["Cha",],trend_formula = "~ splines::ns(pseudotime, df=4)") #check that it worked


###__Differential expression analysis
#first identify variable/differentially-expressed genes along pseudotime ("principal graph") or localized in UMAP ("knn")
#we will use the combination of the two

gene_expression_L3 = as.data.frame(Matrix::rowSums(exprs(cdsMi1[,partL3])>0)/ncol(cdsMi1[,partL3])) ###counts the fraction of cells in which each gene is detected
colnames(gene_expression_L3) <- "expression"
gene_expression_L3$gene_short_name <- row.names(gene_expression_L3)

gene_expression = as.data.frame(Matrix::rowSums(exprs(cdsMi1)>0)/ncol(cdsMi1))
colnames(gene_expression) <- "expression"
gene_expression$gene_short_name <- row.names(gene_expression)

#first let's use "principal graph"
DE_genes_L3_pg <- graph_test(cdsMi1[,partL3], neighbor_graph="principal_graph")
DE_genes_allages_pg <- graph_test(cdsMi1, neighbor_graph="principal_graph")


#removing genes expressed in less than 2.5% of the cells and with a q value less than 0.05
signif_genes_L3_pg <- DE_genes_L3_pg %>% merge(gene_expression_L3) %>%
                                  filter(q_value <0.05) %>%
                                  filter(expression >0.025)
signif_genes_L3_pg <- arrange(signif_genes_L3_pg, desc(morans_I))
signif_genes_L3_pg$gene_short_name <- as.character(signif_genes_L3_pg$gene_short_name)

signif_genes_pg <- DE_genes_allages_pg %>% merge(gene_expression) %>%
  filter(q_value <0.05) %>%
  filter(expression >0.025)
signif_genes_pg <- arrange(signif_genes_pg, desc(morans_I))
signif_genes_pg$gene_short_name <- as.character(signif_genes_pg$gene_short_name)

#then same with "knn"
DE_genes_L3_knn <- graph_test(cdsMi1[,partL3], neighbor_graph="knn")
DE_genes_allages_knn <- graph_test(cdsMi1, neighbor_graph="knn")

#removing genes expressed in less than 2.5% of the cells and with a q value less than 0.05
signif_genes_L3_knn <- DE_genes_L3_knn %>% merge(gene_expression_L3) %>%
  filter(q_value <0.05) %>%
  filter(expression >0.025)
signif_genes_L3_knn <- arrange(signif_genes_L3_knn, desc(morans_I))
signif_genes_L3_knn$gene_short_name <- as.character(signif_genes_L3_knn$gene_short_name)

signif_genes_knn <- DE_genes_allages_knn %>% merge(gene_expression) %>%
  filter(q_value <0.05) %>%
  filter(expression >0.025)
signif_genes_knn <- arrange(signif_genes_knn, desc(morans_I))
signif_genes_knn$gene_short_name <- as.character(signif_genes_knn$gene_short_name)


####The genes we will look at are the one differentally expressed in one method or the other
signif_gene_names_L3 <- unique(c(signif_genes_L3_knn$gene_short_name, signif_genes_L3_pg$gene_short_name))
signif_gene_names <- unique(c(signif_genes_knn$gene_short_name, signif_genes_pg$gene_short_name))


##We now group these DE genes into modules
# we use genes and cells from every Mi1 age and not just L3 (works best that way)
set.seed(106) #105
gene_modules <- find_gene_modules(cdsMi1[signif_gene_names,], resolution=c(10^seq(-5,-4), 2*10^-4, 3*10^-4, 4*10^-4, 5*10^-4,6*10^-4, 7*10^-4,8*10^-4),umap.n_neighbors = 45L, random_seed = 42) ##High n_neighbors to capture more global behaviours
ggplot(gene_modules, aes(x = dim_1, y= dim_2, color = module)) + geom_point()

plot_cells(cdsMi1, genes = filter(gene_modules, module %in% c(1,2,3,4,5,6,7)), show_trajectory_graph = F, min_expr = 60)


scale_modules = FALSE
agg_mat = as.matrix(aggregate_gene_expression(cdsMi1,gene_modules, norm_method = "log", scale_agg_values = scale_modules))

colors = c("red","green","blue","yellow","purple","pink","darkgreen","brown","black")
modulesplot <- ggplot()
for(i in 1:nlevels(gene_modules$module)){
  curve_i <- plot_genes_in_pseudotime2(cdsMi1, min_expr = 1,module_mat = agg_mat[i,], scaled = scale_modules, returncurve = T, trend_formula = "~ splines::ns(pseudotime, df=5)")
  modulesplot <- modulesplot + geom_line(aes(x = pseudotime, y = expectation), data = curve_i, colour = colors[i]) +
    theme(axis.line = element_line(size = 0.3, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white")) +
          scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
          scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
}

#tiff("modulesplot.tiff", width = 2000, height=2000, res=300)
modulesplot
#dev.off()

scale_modules = TRUE
agg_mat = as.matrix(aggregate_gene_expression(cdsMi1,gene_modules, norm_method = "log", scale_agg_values = scale_modules))

colors = c("red","green","blue","yellow","purple","pink","darkgreen","brown","black")

moduleCurves = data.frame()
moduleFig <- plot_ly() %>%
  layout(xaxis = list(showgrid = F),
         yaxis = list(showgrid = F))

for(i in 1:nlevels(gene_modules$module)){
  curve <- plot_genes_in_pseudotime2(cdsMi1, min_expr = 1,module_mat = agg_mat[i,], scaled = scale_modules, returncurve = T, trend_formula = "~ splines::ns(pseudotime, df=5)")
  curve$module = i
  curve <- curve[order(curve$pseudotime),]
  moduleCurves = rbind(moduleCurves, curve)
  
  moduleFig <- moduleFig %>% add_trace(data = curve,type = "scatter", mode="lines",x=~pseudotime, y=~expectation, name = paste("Module ",i))
  
}

moduleFig
#orca(moduleFig, file = "moduleFig.png")



#visualize individual modules
plot_genes_in_pseudotime2(cdsMi1, min_expr = 1,module_mat = agg_mat[1,], scaled = scale_modules, returncurve = F, trend_formula = "~ splines::ns(pseudotime, df=5)")


###___Save modules to analyze them with Panther

conversiontable = read.csv("ConversionTable.csv") #names of genes need to be harmonized
#there is a twist : sometimes Panther recognizes the reference_id , and sometimes the converted_id. To maximize the number of identified genes you put both columns into Panther and it will eliminate doublets

for (i in 1:nlevels(gene_modules$module)){
  module_i <- conversiontable[conversiontable$reference_symbol %in% gene_modules$id[gene_modules$module == i], c("reference_symbol", "reference_id","converted_id")]
  module_i_panther = c(as.character(module_i$reference_id), as.character(module_i$converted_id))
  write.table(module_i_panther, paste("Mi1modules/module",i,"_panther_L3P70.csv",sep=""), row.names = F,col.names = F, quote = F, sep = ",")
}

#you need to input reference genes into panther, they are all the genes from all the modules
reference_genes = conversiontable[conversiontable$reference_symbol %in% gene_modules$id, c("reference_symbol", "reference_id","converted_id")]
reference_genes_panther = c(as.character(reference_genes$reference_id), as.character(reference_genes$converted_id))
write.table(reference_genes_panther, "Mi1modules/reference_genes_panther_L3P70.csv", row.names = F,col.names = F, quote = F, sep = ",")


##if you want to investigate a bump or something like this, the following function can plot individual genes of a module
plot_genes_of_module(cdsMi1, gene_modules$id[gene_modules$module==3], scaled = T, logscale = F)


###__Then look at individual GO terms along trajectory
goterms = read.csv("mart_export.txt", stringsAsFactors = F)
idtoterms = data.frame(id = names(Term(GOTERM)), term = Term(GOTERM), stringsAsFactors = F)

goToPlot = c("DNA replication","ribosome assembly","axon guidance","dendrite development","chemical synaptic transmission",
              "regulation of membrane potential")

goCurves = data.frame()
goFig <- plot_ly()%>%
  layout(xaxis = list(showgrid = F),
         yaxis = list(showgrid = F))
for (name in goToPlot){
  idofinterest = filter(idtoterms, term == name)$id
  allidsofinterest = c(idofinterest,GOBPOFFSPRING[[idofinterest]])
  
  genesofinterest = goterms$Gene.stable.ID[goterms$GO.term.accession %in% allidsofinterest]
  genesofinterest = conversiontable$reference_symbol[conversiontable$reference_id %in% genesofinterest]
  gene.set = as.character(genesofinterest)
  gene.set = intersect(gene.set, signif_gene_names)
  
  GOmodule = data.frame(id = gene.set, module = 1)
  agg_matGO = as.matrix(aggregate_gene_expression(cdsMi1,GOmodule, norm_method = "log", scale_agg_values =T))
  curve = plot_genes_in_pseudotime2(cdsMi1,module_mat =  agg_matGO[1,], min_expr = 1, returncurve = T, scaled = T,trend_formula = "~ splines::ns(pseudotime, df=5)")
  curve$GO_term = name
  curve <- curve[order(curve$pseudotime),]
  goCurves = rbind(goCurves, curve)
  
  goFig <- goFig %>% add_trace(data = curve,type = "scatter", mode="lines",x=~pseudotime, y=~expectation, name = name)
  
}

goFig
orca(goFig, file = "Figure 6Β.png")

###__Now we generalize by looking at all L3 + P15 cells

#gene.set can be any module or go term gene list
gene.set <- filter(gene_modules, module == 1)$id

# Get mean expression of genes of interest per cell
##we use a mean here because gene expression is already scaled by Seurat. No gene should overcome the others
mean.exp <- scale(colMeans(x = larva_dataset[["RNA"]]@data[row.names(larva_dataset) %in% gene.set, ], na.rm = TRUE))

# Add mean expression values in 'object@meta.data$gene.set.score'
if (all(names(x = mean.exp) == rownames(x = larva_dataset@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  larva_dataset@meta.data$gene.set.score <- mean.exp
}

# Plot mean expression using Seurat::FeaturePlot()
q = quantile(larva_dataset[["gene.set.score"]]$gene.set.score, c(.05, .95))
FeaturePlot(object = larva_dataset, feature= "gene.set.score", pt.size = 0.05, min.cutoff =  q[[1]], max.cutoff = q[[2]]) #+ ggtitle(termofinterest)

###save plots
for(i in 1:nlevels(gene_modules$module)){
  gene.set <- filter(gene_modules, module == i)$id
  mean.exp <- scale(colMeans(x = larva_dataset[["RNA"]]@data[row.names(larva_dataset) %in% gene.set, ], na.rm = TRUE))
  if (all(names(x = mean.exp) == rownames(x = larva_dataset@meta.data))) {
    cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
        "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
    larva_dataset@meta.data$gene.set.score <- mean.exp
  }
  q = quantile(larva_dataset[["gene.set.score"]]$gene.set.score, c(.05, .95))
  ggsave(
    paste("module_",i,"_all_L3_P15.png"),
    FeaturePlot(object = larva_dataset, feature= "gene.set.score", pt.size = 0.05, min.cutoff =  q[[1]], max.cutoff = q[[2]]) + ggtitle(paste("module",i,sep="_")),
    width = 12,
    height = 12,
    dpi = 1000
  )
}


###__Doing the same with specific GO terms instead of modules
name = "axon guidance"
idofinterest = filter(idtoterms, term == name)$id
allidsofinterest = c(idofinterest,GOBPOFFSPRING[[idofinterest]])

genesofinterest = goterms$Gene.stable.ID[goterms$GO.term.accession %in% allidsofinterest]
genesofinterest = conversiontable$reference_symbol[conversiontable$reference_id %in% genesofinterest]
gene.set = as.character(genesofinterest)
gene.set = intersect(gene.set, signif_genes$gene_short_name)

mean.exp <- scale(colMeans(x = larva_dataset[["RNA"]]@data[row.names(larva_dataset) %in% gene.set, ], na.rm = TRUE))

if (all(names(x = mean.exp) == rownames(x = larva_dataset@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  larva_dataset@meta.data$gene.set.score <- mean.exp
}

q = quantile(larva_dataset[["gene.set.score"]]$gene.set.score, c(.3, .95))
FeaturePlot(object = larva_dataset, feature= "gene.set.score", raster = F, pt.size = 0.05, cols = c("blue2", "firebrick1"),min.cutoff =  q[[1]], max.cutoff = q[[2]]) + ggtitle(name)

#automating saving

for (termofinterest in goToPlot){
  idofinterest = filter(idtoterms, term == termofinterest)$id
  allidsofinterest = c(idofinterest,GOBPOFFSPRING[[idofinterest]])
  
  genesofinterest = goterms$Gene.stable.ID[goterms$GO.term.accession %in% allidsofinterest]
  genesofinterest = conversiontable$reference_symbol[conversiontable$reference_id %in% genesofinterest]
  genesofinterest = as.character(genesofinterest)
  genesofinterest = intersect(genesofinterest, rownames(larva_dataset))

  mean.exp <- scale(colMeans(x = larva_dataset[["RNA"]]@data[row.names(larva_dataset) %in% genesofinterest, ], na.rm = TRUE))
  
  if (all(names(x = mean.exp) == rownames(x = larva_dataset@meta.data))) {
    larva_dataset@meta.data$gene.set.score <- mean.exp
  }
  
  q = quantile(larva_dataset[["gene.set.score"]]$gene.set.score, c(.3, .95))
  
  ggsave(
    paste(termofinterest,"_all_L3_P15.png",sep=''),
    FeaturePlot(object = larva_dataset, feature= "gene.set.score",cols = c("blue2", "firebrick1"),raster = F, pt.size = 0.3, min.cutoff = q[[1]], max.cutoff = q[[2]]) + ggtitle(termofinterest),
    width = 12,
    height = 12,
    dpi = 300
  )
}



### Figure 6D ###

tiff("Figure 6D.tiff",width = 6000, height=2000, res=300)
LayeredFeaturePlot(object = larva_dataset, features = c("Cha", "VGlut", "Gad1"), colors = c("green", "lightblue", "red"),reduction = "umap", pt.brightness = 10, min.cutoff = 0, max.cutoff = 2 )
dev.off()


### prep for Figure 6F-F' ###
load("human_snSeq.Rdata")

#### We will focus on Cortical plate cells

CP_cells <- subset(merged, brain_region == "CP")
CP_cells <- NormalizeData(CP_cells, normalization.method = "LogNormalize", scale.factor = 10000)
CP_cells <- FindVariableFeatures(CP_cells, selection.method = "vst", nfeatures = 2000)
CP_cells <- ScaleData(CP_cells)
CP_cells <- RunPCA(CP_cells, features = VariableFeatures(object = CP_cells), npcs = 50)
CP_cells <- RunUMAP(CP_cells, dims = 1:25)
DimPlot(CP_cells, group.by = "weeks")

rm(CP_cells) ##remove this seurat object for memory purposes

#### We continue by looking at CP week 19 as it presents the fullest trajectory
CP19_seurat <- subset(merged, Dissection == "11_CP") #matching of numbers : 16 = 18 weeks; 11 = 19 weeks; 4 = 23 weeks; 8 = 24 weeks
CP19_seurat <- NormalizeData(CP19_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
CP19_seurat <- FindVariableFeatures(CP19_seurat, selection.method = "vst", nfeatures = 2000)
CP19_seurat <- ScaleData(CP19_seurat)
CP19_seurat <- RunPCA(CP19_seurat, features = VariableFeatures(object = CP19_seurat), npcs = 30)
CP19_seurat <- RunUMAP(CP19_seurat, dims = 1:12, n.neighbors = 45L)
DimPlot(CP19_seurat)
CP19_seurat <- FindNeighbors(CP19_seurat)
CP19_seurat <- FindClusters(CP19_seurat, random.seed = 42, resolution = 0.5)
DimPlot(CP19_seurat, label = T)
FeaturePlot(CP19_seurat, features = "NEUROD2")
##we will keep only clusters from the main neuronal trajectory, they are 0,1,2,3,5 and 7 for this run
CP19_seurat_sub <- subset(CP19_seurat, seurat_clusters %in% c(0,1,2,3,5,7))
CP19_seurat_sub <- RunUMAP(CP19_seurat_sub, dims = 1:10, n.neighbors = 45L, min.dist = 0.2) #we rerun the UMAP once so that no cell is leftover far from the other because of the subset
DimPlot(CP19_seurat_sub)

tiff("Figure 6F.tiff", width = 2000, height=2000, res=300)
LayeredFeaturePlot(CP19_seurat_sub, reduction = "umap", features = c("rna_PAX6", "rna_EOMES", "rna_NEUROD2"), colors = c("green", "purple", "red"),pt.size = 1)
dev.off()

## Giving data to Monocle3
data <- GetAssayData(object = CP19_seurat_sub, slot = "counts",assay = "RNA")
pData <- data.frame(cell_id = colnames(data), row.names = colnames(data), weeks = factor(CP19_seurat_sub@meta.data[colnames(data),"weeks"]), nCount_RNA =CP19_seurat_sub@meta.data[colnames(data),"nCount_RNA"] )
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
CP19_cds <- new_cell_data_set(data, pData, fData)
CP19_cds <- preprocess_cds(CP19_cds, num_dim = 30)

## As Monocle3 doesn't give as good a trajectory for the CP neurons (probably because of the different preprocessing method),
## we use the UMAP from Seurat and inject it into Monocle3 before proceeding to trajectory detection and analysis

CP19_cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- CP19_seurat_sub@reductions[["umap"]]@cell.embeddings
rownames(CP19_cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
colnames(CP19_cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]]) <- NULL

## clustering and trajectories
set.seed(1502)
CP19_cds <- cluster_cells(CP19_cds, resolution = 0.01)

plot_cells(CP19_cds, color_cells_by = "partition") ##the neuronal trajectory should be one partition
plot_cells(CP19_cds, genes = c("PAX6","EOMES","TBR1")) ##PAX6 marks stem cells, EOMES intermediate progenitors, TBR1 neurons

CP19_cds <- learn_graph(CP19_cds, learn_graph_control = list(L1.sigma= 0.4,minimal_branch_len = 30))
tiff("Figure S10C.tiff", width = 2000, height=2000, res=300)
plot_cells(CP19_cds, cell_size = 0.5, show_trajectory_graph = T, label_roots = F, label_leaves = F, label_cell_groups = F, label_groups_by_cluster = F, labels_per_group = F)
dev.off()
CP19_cds <- order_cells(CP19_cds) ##choosing the beginning point of the curve (see PAX6 expression)


CP19_cds <- detect_genes(CP19_cds)
CP19_cds <- CP19_cds[rownames(CP19_cds)[rowData(CP19_cds)$num_cells_expressed>ncol(CP19_cds)*0.01],] ##only keep genes expressed in at least 1% of cells
gene_significance_pg <- graph_test(CP19_cds, neighbor_graph="principal_graph")
gene_significance_knn <- graph_test(CP19_cds, neighbor_graph="knn")

signif_genes <- c(filter(gene_significance_pg, q_value <0.05),filter(gene_significance_knn, q_value <0.05))


set.seed(53) #for reproducibility
#gene_modules <- find_gene_modules(CP19_cds[as.character(signif_genes$gene_short_name),], resolution=c(10^seq(-5,-4), 2*10^-4, 3*10^-4, 4*10^-4, 5*10^-4, 6*10^-4), random_seed = 42L, umap.fast_sgd = F, init = "normlaplacian")
gene_modules <- find_gene_modules(CP19_cds[as.character(signif_genes$gene_short_name),], resolution=c(10^seq(-5,-4), 2*10^-4, 3*10^-4, 4*10^-4, 5*10^-4, 6*10^-4, 7*10^-4, 8*10^-4), random_seed = 42L, umap.fast_sgd = F, init = "normlaplacian")
nlevels(gene_modules$module) #7 modules
ggplot(gene_modules, aes(x = dim_1, y= dim_2, color = module)) + geom_point()

scale_modules = F
agg_mat = as.matrix(aggregate_gene_expression(CP19_cds,gene_modules, norm_method = "log", scale_agg_values = scale_modules))
colors = c("red","green","blue","yellow","purple","pink","darkgreen","brown","black")
modulesplot <- ggplot()
for(i in 1:nlevels(gene_modules$module)){
  curve_i <- plot_genes_in_pseudotime2(CP19_cds, min_expr = 1,module_mat = agg_mat[i,], scaled = scale_modules, returncurve = T, trend_formula = "~ splines::ns(pseudotime, df=5)")
  if(i ==6){
    curve_i <- plot_genes_in_pseudotime2(CP19_cds, min_expr = 1,module_mat = agg_mat[i,], scaled = scale_modules, returncurve = T, trend_formula = "~ splines::ns(pseudotime, df=10)") ###we need an exception to see the very short change
  }
  modulesplot <- modulesplot + geom_line(aes(x = pseudotime, y = expectation), data = curve_i, colour = colors[i])
}
modulesplot


scale_modules = TRUE
agg_mat = as.matrix(aggregate_gene_expression(CP19_cds,gene_modules, norm_method = "log", scale_agg_values = scale_modules))

colors = c("red","green","blue","yellow","purple","pink","darkgreen","brown","black")

moduleCurves = data.frame()
moduleFig <- plot_ly() %>%
  layout(xaxis = list(showgrid = F),
         yaxis = list(showgrid = F))

for(i in 1:nlevels(gene_modules$module)){
  curve <- plot_genes_in_pseudotime2(CP19_cds, min_expr = 1,module_mat = agg_mat[i,], scaled = scale_modules, returncurve = T, trend_formula = "~ splines::ns(pseudotime, df=5)")
  curve$module = i
  curve <- curve[order(curve$pseudotime),]
  moduleCurves = rbind(moduleCurves, curve)
  
  moduleFig <- moduleFig %>% add_trace(data = curve,type = "scatter", mode="lines",x=~pseudotime, y=~expectation, name = paste("Module ",i))
  
}
moduleFig
orca(moduleFig, file = "Figure S5C'.png")

#save for GO analysis
write.table(gene_modules$id,"modules/reference_genes.csv", row.names=F, quote = F, col.names = F)
for(i in 1:nlevels(gene_modules$module)){
  write.table(gene_modules$id[gene_modules$module==i], paste("modules/human_module_",i,".csv",sep=''), row.names = F,col.names =F, quote = F)
}

##GO
goToPlot = c("DNA replication","ribosome assembly","axon guidance","dendrite development","chemical synaptic transmission",
             "regulation of membrane potential")
goCurves = data.frame()
goFig <- plot_ly()%>%
  layout(xaxis = list(showgrid = F),
         yaxis = list(showgrid = F))
for (termofinterest in goToPlot){
  idofinterest = filter(idtoterms, term == termofinterest)$id
  allidsofinterest = c(idofinterest,GOBPOFFSPRING[[idofinterest]])
  
  genesofinterest = goterms$Gene.name[goterms$GO.term.accession %in% allidsofinterest]
  genesofinterest = as.character(genesofinterest)
  genesofinterest = intersect(genesofinterest,signif_genes$gene_short_name)
  
  GOmodule = data.frame(id = genesofinterest, module = 1)
  agg_matGO = as.matrix(aggregate_gene_expression(CP19_cds,GOmodule, norm_method = "log", scale_agg_values =T))
  curve = plot_genes_in_pseudotime2(CP19_cds,module_mat =  agg_matGO[1,], min_expr = 1, returncurve = T, scaled = T,trend_formula = "~ splines::ns(pseudotime, df=5)")
  curve$GO_term = termofinterest
  curve <- curve[order(curve$pseudotime),]
  goCurves = rbind(goCurves, curve)
  
  goFig <- goFig %>% add_trace(data = curve,type = "scatter", mode="lines",x=~pseudotime, y=~expectation, name = termofinterest)
  
}

goFig
orca(goFig, file = "Figure 6F'.png")


tiff("Figure S10D.tiff", width = 2000, height=2000, res=300)
FeaturePlot(CP19_seurat_sub, reduction = "umap", features = c("rna_PAX6"))
dev.off()

tiff("Figure S10D'.tiff", width = 2000, height=2000, res=300)
LayeredFeaturePlot(CP19_seurat_sub, reduction = "umap", features = c("rna_FBXO32", "rna_HOPX", "rna_NEUROD2"), colors = c("green", "red"),pt.size = 1)
dev.off()



