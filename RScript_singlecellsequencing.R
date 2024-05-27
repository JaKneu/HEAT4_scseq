#R Script Cohort 1, HF and Ctrl
################################################################################################
#R version 4.1.3 (2022-03-10)
#RStudio version 2022.02.0+443


#Load libraries
################################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scales)
library(cowplot)
library(sctransform)
library(glmGamPoi)


#Load data and create seurat objects
###############################################################################################
#ctrl1, exemplary for ctrl1-ctrl3
ctrl1.data <- Read10X(data.dir = ".../filtered_feature_bc_matrix")
ctrl1 <- CreateSeuratObject(counts = ctrl1.data, project = "ctrl1", min.cells=3, min.features=200)
ctrl1$sampleid <- "Ctrl1"
ctrl1$condition <- "Ctrl"
ctrl1[["percent.mt"]] <- PercentageFeatureSet(ctrl1, pattern = "^MT-")
VlnPlot(ctrl1, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)
ctrl1 <- subset(ctrl1, subset = nFeature_RNA > 200  & percent.mt < 10)
ctrl1 <- NormalizeData(ctrl3, normalization.method = "LogNormalize", scale.factor = 10000)
ctrl1 <- SCTransform(ctrl1, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE)
ctrl1 <- RenameCells(ctrl1, add.cell.id = "ctrl1")
 
#pat1, exemplary for pat1-pat8
pat1.data <- Read10X(data.dir = ".../filtered_feature_bc_matrix")
pat1 <- CreateSeuratObject(counts = pat1.data, project = "pat1", min.cells=3, min.features=200)
pat1$sampleid <- "Pat1"
pat1$condition <- "Pat"
pat1[["percent.mt"]] <- PercentageFeatureSet(pat1, pattern = "^MT-")
VlnPlot(pat1, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)
pat1 <- subset(pat1, subset = nFeature_RNA > 200  & percent.mt < 10)
pat1 <- NormalizeData(pat1, normalization.method = "LogNormalize", scale.factor = 10000)
pat1 <- SCTransform(pat1, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE)
pat1 <- RenameCells(pat1, add.cell.id = "pat1")


#Integration and dimensional reduction
###############################################################################################
seuratobject.list <- c(ctrl1, ctrl2, ctrl3, pat1, pat2, pat3, pat4, pat5, pat6, pat7, pat8)
features <- SelectIntegrationFeatures(object.list = seuratobject.list)
seuratobject.list <- PrepSCTIntegration(object.list = seuratobject.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seuratobject.list, reference = c(1, 4, 6), normalization.method = "SCT", dims = 1:50, anchor.features = features)
immune.integrated <- IntegrateData(anchorset = anchors, dims = 1:50,normalization.method = "SCT", verbose=TRUE)

DefaultAssay(immune.integrated)
immune.integrated <- RunPCA(immune.integrated, verbose = FALSE)


#Test cell cycle bias
###############################################################################################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
immune.integrated <- CellCycleScoring(immune.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(immune.integrated)
#no regression needed


#Clustering
###############################################################################################
ElbowPlot(immune.integrated, ndims=40)
immune.integrated <- RunUMAP(immune.integrated, reduction = "pca", dims = 1:12)
immune.integrated <- FindNeighbors(immune.integrated, reduction = "pca", dims = 1:12)
immune.integrated <- FindClusters(immune.integrated, resolution = c(0.4, 0.6, 0.8))

#finding best resolution
Idents(object=immune.integrated) <- "integrated_snn_res.0.4"
dimplot1 <- DimPlot(immune.integrated, reduction="umap", label = TRUE, label.size = 6)
Idents(object=immune.integrated) <- "integrated_snn_res.0.6"
dimplot2 <- DimPlot(immune.integrated, reduction="umap", label = TRUE, label.size = 6)
Idents(object=immune.integrated) <- "integrated_snn_res.0.8"
dimplot3 <- DimPlot(immune.integrated, reduction="umap", label = TRUE, label.size = 6)
dimplot1 + dimplot2 + dimplot3
#choose 0.4
Idents(object=immune.integrated) <- "integrated_snn_res.0.4"

#FIGURE S1A
FS1A <- DimPlot(immune.integrated, label=TRUE)


#Number of sequenced cells
###############################################################################################
#FIGURE 1A
table(immune.integrated$sampleid)


#Assigning cell type identity to clusters
###############################################################################################
#FIGURE S1B-E
DefaultAssay(immune.integrated) <- "RNA"
FS1B_E <- FeaturePlot(object = immune.integrated, features = c("CD14", "FCGR3A", "CD79A", "CD3E"), cols=c("grey", "#2882BD"))

#FIGURE S1F
markers.to.plot <- c("CD14", "LGALS1", "S100A12", "S100A9", "S100A8", "FCN1", "MS4A7", "CDKN1C", "HMOX1", "CD3G", "CD2", "CD3E", "CD3D", "LDHB", "CD8B", "CD8A", "CD22", "CD40", "CD72", "CD19", "MS4A1", "CD79B", "CD79A", "ITGB3", "ITGA2B", "KLRD1", "KIR3DL1", "KLRC1", "IL2RB", "FCGR3A", "SPON2", "XCL2", "FGFBP2", "IL3RA", "SERPINF1", "ITM2C", "CD1C", "FCER1A", "ITGAX", "CD33", "CD34")
FS1F <- DotPlot(immune.integrated, features = markers.to.plot, cols = c("grey", "#2882BD"), dot.scale = 10) + 
  RotatedAxis()+ scale_y_discrete(limits=c("8", "9", "4", "0", "6", "14", "3", "2", "5", "7", "12", "1", "11", "10", "13"))

#SINGLE_R
library(SingleR)
library(celldex)
exp<-GetAssayData(immune.integrated, slot = "data", assay = "RNA")
#Load reference data
refs<-list(ref1 = HumanPrimaryCellAtlasData(), ref2 = BlueprintEncodeData(), ref3 = DatabaseImmuneCellExpressionData(), ref4 = NovershternHematopoieticData(), ref5 = MonacoImmuneData())
ref1.data <- HumanPrimaryCellAtlasData()
ref2.data <- BlueprintEncodeData()
ref3.data <- DatabaseImmuneCellExpressionData()
ref4.data <- NovershternHematopoieticData()
ref5.data <- MonacoImmuneData()
#Perform predictions, exemplary for HumanPrimaryCellAtlasData
HumanPrimaryCellAtlasData.predictions <- SingleR(test=exp, assay.type.test=1, ref=ref1.data, labels=ref1.data$label.main)
immune.integrated[["HumanPrimaryCellAtlasData"]] <- HumanPrimaryCellAtlasData.predictions$labels

#Renaming of clusters
new.cluster.ids <- c("Monocytes","NK Cells","T Cells", "T Cells", "Monocytes","T Cells",  "Monocytes", "B Cells", "Monocytes","Monocytes","Myeloid Dendritic Cells", "Plasmacytoid Dendritic Cells", "Mixed", "Progenitor Cells", "T Cells")
names(new.cluster.ids) <- levels(immune.integrated)
immune.integrated <- RenameIdents(immune.integrated, new.cluster.ids)
immune.integrated <- StashIdent(object = immune.integrated, save.name = 'cluster.ident')


#Divide monocytes into subtypes
###############################################################################################
Idents(object = immune.integrated) <- immune.integrated@meta.data$cluster.ident
Idents(immune.integrated, WhichCells(object = immune.integrated, ident = "Monocytes", expression = FCGR3A > 0, slot = 'data')) <- 'CD16Pos Monocytes'
immune.integrated[["cluster.ident.0"]] <- Idents(object = immune.integrated)
DimPlot(immune.integrated, label=TRUE)

#Downsample the number of cells per condition for comparison plots
table(immune.integrated$condition)
#Ctrl 24,573, Pat 59,016
Idents(object = immune.integrated) <- immune.integrated@meta.data$condition
set.seed(12)
immune.integrated.20000 <- subset(x = immune.integrated, downsample = 20000)

#FIGURE 1B
Idents(object = immune.integrated.20000) <- immune.integrated.20000@meta.data$cluster.ident.0
F1B <- DimPlot(immune.integrated.20000, split.by="condition")


#Cell type percentages and stacked bar plots
###############################################################################################
#FIGURE S1G
immune.integrated.dataS1G <- immune.integrated@meta.data %>% group_by(sampleid, cluster.ident.0) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)
v_factor_levels <- c("T Cells","NK Cells", "Mixed",  "B Cells",   "Progenitor Cells","Plasmacytoid Dendritic Cells","Myeloid Dendritic Cells","CD16Pos Monocytes",  "Monocytes")
FS1G <- ggplot(immune.integrated.dataS1G, aes(x = sampleid, y = percent, fill = factor(cluster.ident.0, levels = v_factor_levels), order = cluster.ident.0))+
  geom_bar(stat = "identity")+scale_x_discrete(limits=c("Pat8", "Pat7", "Pat6", "Pat5", "Pat4", "Pat3", "Pat2", "Pat1", "Ctrl3", "Ctrl2", "Ctrl1"))
FS1G + coord_flip() +
  geom_col(position = position_stack(reverse = TRUE))+ scale_fill_manual(values=c("grey", "light grey", "grey55", "grey75", "grey40", "grey80", "grey45", "#2882BD", "#6CB3E1"))+theme(panel.background = element_rect(fill='transparent'), 
                                                                                                                                                                                       panel.border = element_blank(),
                                                                                                                                                                                       panel.grid.major = element_blank(),
                                                                                                                                                                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))                                                                                                                                                                                             

#FIGURE 1D + S1H
write.csv(immune.integrated.dataS1G, file = ".../scCohort1_PercentCells_sampleid.csv")

#FIGURE 1C
immune.integrated.data1C <- immune.integrated@meta.data %>% group_by(condition, cluster.ident.0) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)
v_factor_levels <- c("T Cells","NK Cells", "Mixed",  "B Cells",   "Progenitor Cells","Plasmacytoid Dendritic Cells","Myeloid Dendritic Cells","CD16Pos Monocytes",  "Monocytes")
F1C <- ggplot(immune.integrated.data1C, aes(x = condition, y = percent, fill = factor(cluster.ident.0, levels = v_factor_levels), order = cluster.ident.0))+
  geom_bar(stat = "identity")+scale_x_discrete(limits=c("Pat", "Ctrl"))
F1C + coord_flip() +geom_col(position = position_stack(reverse = TRUE)) + scale_fill_manual(values=c("grey", "light grey", "grey55", "grey75", "grey40", "grey80", "grey45", "#2882BD", "#6CB3E1"))+theme(panel.background = element_rect(fill='transparent'), 
                                                                                                                                                                                                         panel.border = element_blank(),
                                                                                                                                                                                                         panel.grid.major = element_blank(),
                                                                                                                                                                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))   
#Regulation of monocyte-specific lncRNAs in CEBPB+ vs. CEBPB- monocytes
###############################################################################################
#FIGURE 1E
FeaturePlot(object = immune.integrated, features = c("FCGR3A", "CD14", "CEBPB"), cols=c("grey", "#2882BD"),keep.scale ="all")& theme(legend.position = "right")

#FIGURE 1F + TABLE S2
##Find lncRNAs expressed in monocytes only
#Divide data set into "monocyte" and "nonmonocyte" (all Clusters besides Monocytes)
Idents(object = immune.integrated) <- immune.integrated@meta.data$cluster.ident.0
new.cluster.ids <- c("monocyte","monocyte", "nonmonocyte","nonmonocyte", "nonmonocyte","nonmonocyte", "nonmonocyte", "nonmonocyte", "nonmonocyte")
names(new.cluster.ids) <- levels(immune.integrated)
immune.integrated <- RenameIdents(immune.integrated, new.cluster.ids)
DimPlot(immune.integrated, reduction = "umap", cols=c( "#2882BD","grey"))
immune.integrated <- StashIdent(object = immune.integrated, save.name = 'monoornot')
#Load lncRNA names of 56,846 lncRNAs from Excel ("lncipedia_Gennamesonly_56846"), sheet: Gene.names
library(readxl)
lncipedia_Gennamesonly <- read_excel(".../lncipedia_Gennamesonly_56846.xlsx", sheet = "Gene.names")
lncipedia <- lncipedia_Gennamesonly[[1]]
#Find lncRNAs expressed in the data set
lncipedia.used=intersect(rownames(immune.integrated), lncipedia)
#found 18688
#Get percentage of monocytes and nonmonocytes expressing the found 18688 lncRNAs
library(scCustomize)
#exemplary:
percent.lncRNAs0 <- Percent_Expressing(immune.integrated, features = c(lncipedia.used[1:1000]), assay = "RNA", threshold = 0)
percent.lncRNAs0$Gene_ID <- row.names(percent.lncRNAs0)
write.csv(percent.lncRNAs0, file = ".../percent.lncRNAs0.csv")
#...
percent.lncRNAs18 <- Percent_Expressing(immune.integrated, features = c(lncipedia.used[18001:18688]), assay = "RNA", threshold = 0)
percent.lncRNAs18$Gene_ID <- row.names(percent.lncRNAs18)
write.csv(percent.lncRNAs18, file = ".../percent.lncRNAs18.csv")
#in Excel: ("Monocyte_specific_lncRNAs")
#Load all 19 .csv-files and combine them in one file ??? "Monocyte_specific_lncRNAs"
#filter: keep only lncRNAs which are expressed in >2% monocyte and <0.5% nonmonocyte
#Load filtered Excel file, sheet: monospecificlncnrnas 
library(readxl)
Monocyte_specific_lncRNAs <- read_excel(".../Monocyte_specific_lncRNAs.xlsx", sheet = "monospecificlncrnas")
monospecificlncRNAs <- Monocyte_specific_lncRNAs[[1]]
#list with 137 monocyte-specific lncRNAs
##Divide monocytes into CEBPB+ and CEBPB- monocytes
immune.integrated[["CEBPB_RNA"]] <- immune.integrated@assays[["RNA"]]@data["CEBPB",]
immune.integrated[["CEBPB_PN"]] <- {ifelse(immune.integrated@assays[["RNA"]]@data["CEBPB",] > 0, "CEBPBPos", "CEBPBNeg")}
immune.integrated[["ident.CEBPB"]] <-paste(immune.integrated@meta.data$monoornot, immune.integrated@meta.data$CEBPB_PN, sep = "_")
Idents(object = immune.integrated) <- immune.integrated@meta.data$ident.CEBPB
DimPlot(immune.integrated)
monocytes <- subset(immune.integrated, idents = c("monocyte_CEBPBPos", "monocyte_CEBPBNeg"))
##Calculate average expression of each monocyte-specific lncRNA (#137 in total) in the two clusters: "monocyte_CEBPBPos", "monocyte_CEBPBNeg"
averagelncRNAexpressionCEBPposnegmonos <-  as.data.frame(AverageExpression(monocytes, assays="RNA", slot="data", features=c(monospecificlncRNAs), verbose = TRUE)$RNA)
averagelncRNAexpressionCEBPposnegmonos$Gene_ID <- row.names(averagelncRNAexpressionCEBPposnegmonos)
write.csv(averagelncRNAexpressionCEBPposnegmonos, file = ".../averagelncRNAexpressionCEBPposnegmonos.csv")
#in Excel: ("lncRNAs_in_CEBPposvsneg_Monos")
#calculateFold average_expression_monocyte_CEBPBPos/ average_expression_monocyte_CEBPBPNeg
#sort after Fold (descending)
#check all lncRNAs with Fold enrichment >1.1 for LNCipedia Class
#result: top4: lnc-CEBPB-13, lnc-GHRH-4, lnc-KLF6-20 & lnc-CD300C-1 (Heat4)

#FIGURE 1G
immune.integrated[["lnc.CEBPB.13"]] <- immune.integrated@assays[["RNA"]]@data["lnc-CEBPB-13",]
immune.integrated[["lnc.GHRH.4"]] <- immune.integrated@assays[["RNA"]]@data["lnc-GHRH-4",]
immune.integrated[["lnc.KLF6.20"]] <- immune.integrated@assays[["RNA"]]@data["lnc-KLF6-20",]
immune.integrated[["lnc.CD300C.1"]] <- immune.integrated@assays[["RNA"]]@data["lnc-CD300C-1",]
write.csv(immune.integrated@meta.data, file = ".../CEBPBPosNeg_top4lncRNAExpression.csv")
#in Excel: ("lncRNAs_in_CEBPposvsneg_Monos")
#filter for Monocytes Only, calculate mean and plot data

#FIGURE 1H
FeaturePlot(immune.integrated, feature=c("lnc-CEBPB-13", "lnc-GHRH-4", "lnc-KLF6-20", "lnc-CD300C-1"),cols=c( "grey","#2882BD"), pt.size=2, order=TRUE, keep.scale ="all")& theme(legend.position = "right")


#Regulation of HEAT4 expresssion
###############################################################################################
#FIGURE 4B
#values already in Excel: ("lncRNAs_in_CEBPposvsneg_Monos")
#command: immune.integrated[["lnc.CD300C.1"]] <- immune.integrated@assays[["RNA"]]@data["lnc-CD300C-1",]

#FIGURE 4C
Idents(object = immune.integrated.20000) <- immune.integrated.20000@meta.data$cluster.ident
FeaturePlot(object = immune.integrated.20000, features = c("lnc-CD300C-1"),  order=TRUE,split.by="condition", cols=c("grey", "#2882BD"), pt.size=3, keep.scale ="all")& theme(legend.position = "right")

#FIGURE 4D
Idents(object = immune.integrated) <- immune.integrated@meta.data$cluster.ident.0
F4Dl_ <- DotPlot(immune.integrated, features=c("CD14", "FCGR3A"), cols = c("grey", "#2882BD"), idents=c("Monocytes", "CD16Pos Monocytes")) 
F4Dl_data <- F4Dl_[["data"]]
F4Dl_data.gg <- ggplot(F4Dl_data, aes(x = features.plot, y = id, color = avg.exp, size = pct.exp)) +  geom_point() +
  scale_size_continuous(limits=c(0, 100), range = c(0,15)) +  scale_colour_gradient(high = "#2882BD", low = "grey") +
  labs( size = "Percent expressed", color = "Average expression") + theme(panel.background = element_rect(fill='transparent'))
F4Dl_data.gg

F4Dr_ <- DotPlot(immune.integrated, features=c("lnc-CD300C-1"),   cols = c("grey", "#2882BD"), idents=c("Monocytes", "CD16Pos Monocytes"))
F4Dr_data <- F4Dr_[["data"]]
F4Dr_data.gg <- ggplot(F4Dr_data, aes(x = features.plot, y = id, color = avg.exp, size = pct.exp)) +  geom_point() +
  scale_size_continuous(limits=c(0, 3), range = c(0,15)) +    scale_colour_gradient(high = "#2882BD", low = "grey") +
  labs( size = "Percent expressed", color = "Average expression") +  theme(panel.background = element_rect(fill='transparent'))
F4Dr_data.gg

#FIGURE 4E
immune.integrated[["Heat4"]] <- {ifelse(immune.integrated@assays[["RNA"]]@data["lnc-CD300C-1",] > 0, "Heat4Pos", "Heat4Neg")}
immune.integrated[["ident.Heat4"]] <-paste(immune.integrated@meta.data$cluster.ident, immune.integrated@meta.data$Heat4, sep = "_")
immune.integrated[["ident.Heat4.condition"]] <-paste(immune.integrated@meta.data$ident.Heat4, immune.integrated@meta.data$condition, sep = "_")
Idents(object = immune.integrated) <- immune.integrated@meta.data$ident.Heat4.condition
DimPlot(immune.integrated, split.by="Heat4")
library(gridExtra)
Cd16all <- VlnPlot(object = immune.integrated, features = c("FCGR3A"), y.max=5,idents=c("Monocytes_Heat4Pos_Ctrl","Monocytes_Heat4Pos_Pat", "Monocytes_Heat4Neg_Ctrl","Monocytes_Heat4Neg_Pat"), group.by = "Heat4", cols=c("grey", "#2882BD"), pt.size = FALSE)
Cd16ctrl <- VlnPlot(object = immune.integrated, features = c("FCGR3A"), y.max=5,idents=c("Monocytes_Heat4Pos_Ctrl", "Monocytes_Heat4Neg_Ctrl"), group.by = "Heat4", cols=c("grey", "#2882BD"), pt.size = FALSE)
Cd16pat <- VlnPlot(object = immune.integrated, features = c("FCGR3A"),y.max=5, idents=c("Monocytes_Heat4Pos_Pat", "Monocytes_Heat4Neg_Pat"), group.by = "Heat4", cols=c("grey", "#2882BD"), pt.size = FALSE)
Cd16all + Cd16ctrl + Cd16pat
F4E <- grid.arrange(Cd16all, Cd16ctrl, Cd16pat, ncol = 3)


#Interaction of the lncRNA Heat4 with S100A9 in the cytoplasm of human monocytes
###############################################################################################
#FIGURE 6H + TABLE S5 
immune.integrated[["ident.condition"]] <-paste(immune.integrated@meta.data$cluster.ident, immune.integrated@meta.data$condition, sep = "_")
Idents(object = immune.integrated) <- immune.integrated@meta.data$ident.condition
control.monocytes <- subset(immune.integrated, idents = c("Monocytes_Ctrl"))
control.monocytes.averages <- AverageExpression(control.monocytes,  slot="counts", assays="RNA", group.by="condition",  verbose=TRUE)
write.csv(control.monocytes.averages[["RNA"]], file = ".../MS_ctrl.monocytes.averages.csv")
#in Excel: ("MS_control.monocytes.averages")
#detect top 50 expressed genes in monocytes of controls
#Venn analysis with MS found proteins bound to Heat4 
#??? S100A8 + S100A9

#FIGURE 6I + 6J as well as FIGURE S5C + S5D
Idents(object = immune.integrated.20000) <- immune.integrated.20000@meta.data$cluster.ident
DimPlot(immune.integrated.20000, split.by="condition")
F6I_FS5C <- FeaturePlot(object = immune.integrated.20000, features = c("S100A8"),split.by="condition",  cols=c("grey", "#2882BD"), pt.size=3, keep.scale ="all")& theme(legend.position = "right")
F6J_FS5D <- FeaturePlot(object = immune.integrated.20000, features = c("S100A9"),split.by="condition", cols=c("grey", "#2882BD"), pt.size=3, keep.scale ="all")& theme(legend.position = "right")
#onlyCtrlshown in 6I + 6J

#QC ribosomal genes
###############################################################################################
#FIGURE S7
#start Code Line 1
#after Line 26:
ctrl1[["percent.rp"]] <- PercentageFeatureSet(ctrl1, pattern = "^RP[SL]")


saveRDS(immune.integrated, file = ".../immune.integrated.rds")
saveRDS(immune.integrated.20000, file = ".../immune.integrated.20000.rds")
saveRDS(control.monocytes, file = ".../control.monocytes.rds")


################################################################################################
################################################################################################
################################################################################################
################################################################################################


#R Script pcDNA and OE Heat4
################################################################################################
#R version 4.1.3 (2022-03-10)
#RStudio version 2022.02.0+443

#Load data and create seurat objects
###############################################################################################
pcDNA.data <- Read10X(data.dir = ".../filtered_feature_bc_matrix")
pcDNA <- CreateSeuratObject(counts = pcDNA.data, project = "pcDNA", min.cells=3, min.features=200)
pcDNA$sampleid <- "pcDNA"
pcDNA[["percent.mt"]] <- PercentageFeatureSet(pcDNA, pattern = "^MT-")
VlnPlot(pcDNA, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)
pcDNA <- subset(pcDNA, subset = nFeature_RNA > 200  & percent.mt < 10)
pcDNA <- NormalizeData(pcDNA, normalization.method = "LogNormalize", scale.factor = 10000)
pcDNA <- SCTransform(pcDNA, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE)
pcDNA <- RenameCells(pcDNA, add.cell.id = "pcDNA")

OEH4.data <- Read10X(data.dir = ".../filtered_feature_bc_matrix")
OEH4 <- CreateSeuratObject(counts = OEH4.data, project = "OEH4", min.cells=3, min.features=200)
OEH4$sampleid <- "OEH4"
OEH4[["percent.mt"]] <- PercentageFeatureSet(OEH4, pattern = "^MT-")
VlnPlot(OEH4, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3)
OEH4 <- subset(OEH4, subset = nFeature_RNA > 200  & percent.mt < 10)
OEH4 <- NormalizeData(OEH4, normalization.method = "LogNormalize", scale.factor = 10000)
OEH4 <- SCTransform(OEH4, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = FALSE)
OEH4 <- RenameCells(OEH4, add.cell.id = "OEH4")


#Integration and dimensional reduction
###############################################################################################
seuratobject.list <- c(pcDNA, OEH4)
features <- SelectIntegrationFeatures(object.list = seuratobject.list)
seuratobject.list <- PrepSCTIntegration(object.list = seuratobject.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seuratobject.list, reference = c(1), normalization.method = "SCT", dims = 1:50, anchor.features = features)
immune.integrated <- IntegrateData(anchorset = anchors, dims = 1:50,normalization.method = "SCT", verbose=TRUE)
DefaultAssay(immune.integrated)
immune.integrated <- RunPCA(immune.integrated, verbose = FALSE)


#Test cell cycle bias
###############################################################################################
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
immune.integrated <- CellCycleScoring(immune.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(immune.integrated)
#no regression needed


#Clustering
###############################################################################################
ElbowPlot(immune.integrated, ndims=40)
immune.integrated <- RunUMAP(immune.integrated, reduction = "pca", dims = 1:14)
immune.integrated <- FindNeighbors(immune.integrated, reduction = "pca", dims = 1:14)
immune.integrated <- FindClusters(immune.integrated, resolution = c(0.4, 0.6, 0.8))

#finding best resolution
Idents(object=immune.integrated) <- "integrated_snn_res.0.4"
dimplot1 <- DimPlot(immune.integrated, reduction="umap", label = TRUE, label.size = 6)
Idents(object=immune.integrated) <- "integrated_snn_res.0.6"
dimplot2 <- DimPlot(immune.integrated, reduction="umap", label = TRUE, label.size = 6)
Idents(object=immune.integrated) <- "integrated_snn_res.0.8"
dimplot3 <- DimPlot(immune.integrated, reduction="umap", label = TRUE, label.size = 6)
dimplot1 + dimplot2 + dimplot3
#choose 0.8
Idents(object=immune.integrated) <- "integrated_snn_res.0.8"

#FIGURE S3B
FS3B <- DimPlot(immune.integrated, label=TRUE)


#Assigning cell type identity to clusters - like PBMC data set cohort 1
###############################################################################################
#FIGURE S3D-G
DefaultAssay(immune.integrated) <- "RNA"
FS3D_G <- FeaturePlot(object = immune.integrated, features = c("CD14", "FCGR3A", "CD79A", "CD3E"), cols=c("grey", "#2882BD"))

#SINGLE_R using PBMC data set of cohort 1
library(SingleR)
library(celldex)
Idents(object=immune.integrated) <- "integrated_snn_res.0.8"
exp<-GetAssayData(immune.integrated, slot = "data", assay = "RNA")
#Load reference data: PBMC data set of cohort 1 = immune.integrated.pbmc
Idents(object = immune.integrated.pbmc) <- immune.integrated.pbmc@meta.data$cluster.ident
reff <- as.SingleCellExperiment(immune.integrated.pbmc)
reff$cluster.ident
predictions <- SingleR(test=exp, ref = reff, labels=reff$cluster.ident, de.method="wilcox")
immune.integrated[["immune.integrated.likepbmc"]] <- predictions$labels
Idents(object = immune.integrated) <- immune.integrated@meta.data$immune.integrated.likepbmc
DimPlot(immune.integrated, label=TRUE)


#Divide monocytes into subtypes
###############################################################################################
Idents(immune.integrated, WhichCells(object = immune.integrated, ident = "Monocytes", expression = FCGR3A > 0, slot = 'data')) <- 'CD16Pos Monocytes'
immune.integrated[["cluster.ident.0"]] <- Idents(object = immune.integrated)
DimPlot(immune.integrated, label=TRUE)

#FIGURES3C
Idents(object = immune.integrated) <- immune.integrated@meta.data$cluster.ident.0
FS3C <- DimPlot(immune.integrated, label=TRUE)

#FIGURES3H
markers.to.plot <- c("CD14", "LGALS1", "S100A12", "S100A9", "S100A8", "FCN1", "MS4A7", "CDKN1C", "HMOX1", "CD3G", "CD2", "CD3E", "CD3D", "LDHB", "CD8B", "CD8A", "CD22", "CD40", "CD72", "CD19", "MS4A1", "CD79B", "CD79A", "ITGB3", "ITGA2B", "KLRD1", "KIR3DL1", "KLRC1", "IL2RB", "FCGR3A", "SPON2", "XCL2", "FGFBP2", "IL3RA", "SERPINF1", "ITM2C", "CD1C", "FCER1A", "ITGAX", "CD33", "CD34")
DotPlot(immune.integrated, features = markers.to.plot, cols = c("grey", "#2882BD"), dot.scale = 10) + 
  RotatedAxis()+ scale_y_discrete(limits=c("CD16Pos Monocytes", "Monocytes", "T Cells", "B Cells", "NK Cells", "Plasmacytoid Dendritic Cells", "Myeloid Dendritic Cells", "Progenitor Cells"))


#RNA expression in pcDNA vs. OE Heat4 monocytes
###############################################################################################
#FIGURES3I
#Downsample the number of cells per condition for comparison plots
table(immune.integrated$sampleid)
#pcDNA 11,081, OE Heat4 11,850
Idents(object = immune.integrated) <- immune.integrated@meta.data$sampleid
set.seed(12)
immune.integrated.10000 <- subset(x = immune.integrated, downsample = 10000)
FS3I <- FeaturePlot(object = immune.integrated.10000, features = c("lnc-CD300C-1"), cols=c("grey", "#2882BD"), pt.size=3, order=TRUE, split.by="sampleid", keep.scale ="all")& theme(legend.position = "right")

#FIGURES3K
Idents(object = immune.integrated) <- immune.integrated@meta.data$sampleid
features.to.plot.upanddown <- c( "S100A4", "CLEC10A", "TGFBR2", "S100A6","FCGR3A",  "TGFBI","CD44", "IDO1", "CD83", "CCL2",  "CCL3", "CCL4")
FS3K <- DotPlot(immune.integrated, features = features.to.plot.upanddown, cols = c("grey", "#2882BD"), dot.scale = 10) + 
  RotatedAxis() + scale_x_discrete(limits=c("S100A4", "CLEC10A", "TGFBR2", "S100A6","FCGR3A",  "TGFBI", "CD44", "IDO1", "CD83", "CCL2",  "CCL3", "CCL4"))+ scale_y_discrete(limits=c("OEH4", "pcDNA"))

#FIGURES6F
FS6F <- FeaturePlot(immune.integrated.10000, features = "S100A8", split.by="sampleid", cols=c("grey", "#2882BD"), pt.size = 1.5, keep.scale ="all")& theme(legend.position = "right")
FS6G <- FeaturePlot(immune.integrated.10000, features = "S100A9", split.by="sampleid", cols=c("grey", "#2882BD"), pt.size = 1.5, keep.scale ="all")& theme(legend.position = "right")


#Overexpression of the lncRNA Heat4 in human monocytes induces anti-inflammatory and vascular healing functions and promotes CD16+ monocyte populations
###############################################################################################
#FIGURE 5I
#Compare Heat4+ and Heat4- monocytes
immune.integrated[["Heat4"]] <- {ifelse(immune.integrated@assays[["RNA"]]@data["lnc-CD300C-1",] > 0, "Heat4Pos", "Heat4Neg")}
immune.integrated[["ident.mono.Heat4"]] <-paste(immune.integrated@meta.data$cluster.ident, immune.integrated@meta.data$Heat4, sep = "_")
Idents(object = immune.integrated) <- immune.integrated@meta.data$ident.mono.Heat4
monocytes <- subset(immune.integrated, idents = c("Monocytes_Heat4Neg", "Monocytes_Heat4Pos"))
monocytes[["ident.mono.Heat4.sampleid"]] <-paste(monocytes@meta.data$ident.mono.Heat4, monocytes@meta.data$sampleid, sep = "_")
Idents(object = monocytes) <- monocytes@meta.data$ident.mono.Heat4.sampleid
DimPlot(monocytes)
F5Iup <- DotPlot(monocytes, feature="FCGR3A", ident=c("Monocytes_Heat4Pos_pcDNA", "Monocytes_Heat4Pos_OEH4"), dot.scale = 5)
F5Iup$data
# ??? pcDNA: 10.71% of Heat4+ monocytes express FCGR3A; OE Heat4: 22.37% of Heat4+ monocytes express FCGR3A

immune.integrated[["Heat4_RNA"]] <- immune.integrated@assays[["RNA"]]@data["lnc-CD300C-1",]
F5Idown <- VlnPlot(monocytes, feature="lnc-CD300C-1", ident=c("Monocytes_Heat4Pos_pcDNA", "Monocytes_Heat4Pos_OEH4"), slot="data")
F5Idown$data
write.csv(F5Idown$data, file = ".../VlnPlotData_Heat4.csv")

saveRDS(immune.integrated, file = ".../immune.integrated.rds")
saveRDS(immune.integrated.10000, file = ".../immune.integrated.20000.rds")
saveRDS(monocytes, file = ".../monocytes.rds")


#QC ribosomal genes
###############################################################################################
#FIGURE S8
#start Code Line 332
#after Line 335:
pcDNA[["percent.rp"]] <- PercentageFeatureSet(pcDNA, pattern = "^RP[SL]")
OEH4[["percent.rp"]] <- PercentageFeatureSet(OEH4, pattern = "^RP[SL]")

