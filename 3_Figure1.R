library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(monocle)
library(reshape2)
library(ggrepel)
library(rlang)
library(ranger)

source("LayeredFeaturePlot_white_background_pt.brightness.R")
source("classifier_utils.R")
OL.combined<-readRDS(file="OL.combined_published.annotations.rds")


#### Figure 1B #####

tiff("Figure 1B.tiff",width = 6000, height=2000, res=300)
LayeredFeaturePlot(OL.combined,features = c("rna_dpn","rna_ase","rna_elav", "rna_repo"),reduction="umap", colors = c("green", "cyan", "purple","red1"), pt.brightness = 10)
dev.off()


#### Figure 1C ####
tiff("Figure 1C.tiff",width = 6000, height=2000, res=300)
DimPlot(object = OL.combined, reduction = "umap", pt.size = 0.3,label = T,label.size = 1,group.by = "L3_vs_P15",cols=c("#00BFC4","#F8766D"))
dev.off()


#### Figure 1D ####
# annotate clusters based on markers and correspondence with P15 data #
lamina<-c("L3_57", "L3_5", "NB2", "L3_65", "L1","L2","L3","L4","L5", "L3_66", "L3_17")
lobula_plate <- c("L3_37",  "Imm1", "L3_35", "C2", "C3", "T2", "T3", "T2a?", "Imm8", "L3_55", "L3_11", "L3_48", "Imm9", "L3_8", "Imm10", "T4-5c/d", "T4-5a/b")
lamina_wide_field<-  c("L3_49","Lawf1","Lawf2")
glia       <-        c("L3_44", "L3_45", "L3_67", "L3_56", "L3_63", "L3_87", "L3_91", "Glia", "Glia-LQ?", "Adult 188", "Not neuron or glia")
lobula <- c("LC4", "LC10a", "LC10b?", "LC16", "64", "111", "Adult 97 (LC6)", "L3_80")
unknown_CB<-                c("L3_15","97","91","78","30","27","153","149","106","102","101","99","98","87","86","85","74","73","68","67","65","50","49","44","5","L3_85", "L3_2", "L3_3", "L3_88", "P70 cluster39", "Adult 208", "PR", "P70 cluster4", "P70 cluster5", "LQ" ,"L3_83","L3_64","L3_89","L3_4","L3_26","L3_77","L3_18","L3_92","L3_93","L3_26","L3_39","L3_50","L3_51","L3_79","L3_90", "LLPC1", "LPC1", "LPLC1", "LPLC2")
#Medulla<-              rest

Idents(object = OL.combined) <- OL.combined@meta.data$Nikos_Neset_ID2

Neuropil<-as.character(OL.combined@active.ident)

for (i in 1:length(Neuropil)) {
  if (Neuropil[i]%in%lamina) {
    Neuropil[i]="lamina"
  }
  else { if (Neuropil[i]%in%lobula_plate) {
    Neuropil[i]="lobula_plate"
  }
    else { if (Neuropil[i]%in%lamina_wide_field) {
      Neuropil[i]="lamina"
    }
      else { if (Neuropil[i]%in%glia) {
        Neuropil[i]="glia"
      }
        else { if (Neuropil[i]%in%unknown_CB) {
          Neuropil[i]="CB"
        }
          else { if (Neuropil[i]%in%lobula) {
            Neuropil[i]="lobula"
          }
          else {
            Neuropil[i]="medulla"
          }
        }
      }
    }
  }}}
Neuropil1<-as.factor(Neuropil)

OL.combined$Neuropil<-Neuropil1

tiff("Figure 1D.tiff",width = 6000, height=2000, res=300)
DimPlot(object = OL.combined, reduction = "umap", pt.size = 1,label = T,label.size = 0,group.by = "Neuropil")
dev.off()


### Figure 1E ###

mne <- readRDS("medulla_ne_reanalysis.rds")

DimPlot(mne, reduction = "umap", pt.size = 1.5) + ggplot2::scale_color_manual(values = c(hsv(90/360, 1, 0.8), hsv(210/360, 1, 0.8), hsv(330/360, 1, 0.8)))
ggsave("Figure1E.png", width = 6.5, height = 6)

LayeredFeaturePlot(mne,
                   c("rna_Optix", "rna_Rx", "rna_Vsx1"),
                   reduction = "umap", background.color = "white", pt.size = 1.5)
ggsave("Figure1E_2.png", width = 5, height = 5)


#### Figure 1F ### 

tiff("Figure 1F.tiff",width = 6000, height=2000, res=300)
LayeredFeaturePlot(OL.combined, features = c("hth","ey","tll", "D", "slp1"),reduction="umap", colors = c("green", "purple", "blue","cyan", "red1"),pt.brightness = 10)
dev.off()
