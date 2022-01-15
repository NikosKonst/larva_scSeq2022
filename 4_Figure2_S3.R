library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(monocle)


#### prep for Figure 2B #####
larvaOL.integrated150<-readRDS("larvaOL.integrated150.rds")

NBs<-SubsetData(larvaOL.integrated150,ident.use=9)

data <- GetAssayData(object = NBs, slot = "counts",assay = "RNA")

pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
pd <- new("AnnotatedDataFrame", data = pData)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)

cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))

expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))

pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))

cds <- cds[,pData(cds)$Total_mRNAs < 1e6]

qplot(Total_mRNAs, data = pData(cds), geom =
        "density")

upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) +
                     2*sd(log10(pData(cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) -
                     2*sd(log10(pData(cds)$Total_mRNAs)))

cds <- cds[,pData(cds)$Total_mRNAs > lower_bound &
             pData(cds)$Total_mRNAs < upper_bound]
cds <- detectGenes(cds, min_expr = 0.1)

L <- log(exprs(cds[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")

cds <-
  reduceDimension(cds, method = 'DDRTree')
cds<-orderCells(cds,reverse = F)

###  Figure 2B ####


tiff("Figure 2B.tiff",width = 4000, height=3000, res=300)
plot_pseudotime_heatmap(cds[c("hth","opa","erm","ey","hbn","esg","Oaz","slp1","slp2","scro","D","B-H1","B-H2","tll"),],cluster_rows = F,
                        show_rownames = T)
dev.off()

### Figure S3 ###

flyTFs<-read.csv("FlyBase_TFs.csv", row.names = 1)
head(flyTFs)
flyTFs2<-flyTFs[,2]
flyTFs3<-as.character(flyTFs2)


for (i in 1:length(flyTFs3)) {
  A=flyTFs3[i]
  cds_subset <- cds[A,]
  tiff(sprintf(sprintf("TFs/test_%s.tiff",A)))
  plot_genes_in_pseudotime(cds_subset) + geom_smooth(size=2)
  dev.off()
  print(i)
}


