library(Seurat)
library(DoubletFinder)
library(ggplot2)


DR = 0.01;

Genes<-read.table("/Users/christopherpenfold/Desktop/Aracely/Pairwise_PGCLC_vs_Amnioids/Plots/GenesMapping.csv",sep=",",header = T)
metaD<-read.table("/Users/christopherpenfold/Downloads/All_meta.tsv",sep="\t",header = T)
D1<-readRDS('/Users/christopherpenfold/Dropbox/PGCs/Processed.dir/ALL_SCE.RDS')
rownames(D1) <- Genes$Gene
humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGCLC_18h"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_18h",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_81", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_81
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_18h_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_18h.rds")


#[1] "PGCLC_18h"        "PGCLC_32h"        "ME_r2"            "PGCLC_40h"        "PGCLC_24h"        "PGCLC_12h_r2"    
#[7] "PGCLC_48h"        "PGCLC_D4_r2"      "DE"               "PGCLC_D4"         "hESC_W5_E8"       "PreME_10h"       
#[13] "PGC_wk6.5_invivo" "4i"               "TC_25"           


Genes<-read.table("/Users/christopherpenfold/Desktop/Aracely/Pairwise_PGCLC_vs_Amnioids/Plots/GenesMapping.csv",sep=",",header = T)
metaD<-read.table("/Users/christopherpenfold/Downloads/All_meta.tsv",sep="\t",header = T)
D1<-readRDS('/Users/christopherpenfold/Dropbox/PGCs/Processed.dir/ALL_SCE.RDS')
rownames(D1) <- Genes$Gene
humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGCLC_32h"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_32h",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_95", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_95
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_32h_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_32h.rds")



humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGCLC_40h"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_40h",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_47", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_47
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_40h_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_40h.rds")





humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGCLC_24h"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_24h",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_153", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_153
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_24h_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_24h.rds")






humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGCLC_12h_r2"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_12h_r2",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_56", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_56
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_24h_doub_r2",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_12h_r2.rds")


humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGCLC_48h"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_48h",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_83", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_83
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_48h_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_48h.rds")









humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGCLC_D4_r2"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_D4_r2",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_50", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_50
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_D4_r2_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_D4_r2.rds")





humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGCLC_D4"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_D4",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_92", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_92
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_D4_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_D4.rds")






humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PreME_10h"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_PreME",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_125", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_125
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_PreME_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_PreME.rds")




humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("hESC_W5_E8"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_hESC_W5_E8",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_36", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_36
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_hESC_W5_E8_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_hESC_W5_E8.rds")





humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("4i"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_4i",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_128", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_128
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_4i_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_4i.rds")






humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("PGC_wk6.5_invivo"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_PGC_wk6.5_invivo",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_171", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_171
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_PGC_wk6.5_invivo_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_PGC_wk6.5_invivo.rds")








humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("DE"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_DE",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_88", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_88
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_DE_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_DE.rds")




humaninvit_data <- as.Seurat(D1, counts = "counts", data = "logcounts")
Idents(humaninvit_data) <- as.character(metaD$CellOrigin) 
humaninvit_data <- subset(humaninvit_data, idents = c("ME_r2"), invert = FALSE)
humaninvit_data <- NormalizeData(humaninvit_data, verbose = FALSE)
humaninvit_data <- FindVariableFeatures(humaninvit_data, selection.method = "vst", nfeatures = 10000)
humaninvit_data <- ScaleData(humaninvit_data)
humaninvit_data <- RunPCA(humaninvit_data, npcs = 30, verbose = FALSE)
humaninvit_data <- RunUMAP(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindNeighbors(humaninvit_data, reduction = "pca", dims = 1:20)
humaninvit_data <- FindClusters(humaninvit_data, resolution = 0.5)

DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_ME_r2",".pdf",sep=""),width = 10, height = 10)


sweep.res.list_18h <- paramSweep_v3(humaninvit_data, PCs = 1:20, sct = FALSE)
sweep.stats_18h <- summarizeSweep(sweep.res.list_18h, GT = FALSE)
bcmvn_18h <- find.pK(sweep.stats_18h)

annotations <- humaninvit_data@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(humaninvit_data@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
humaninvit_data <- doubletFinder_v3(humaninvit_data, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_100", sct = FALSE)
Idents(humaninvit_data) <- humaninvit_data$DF.classifications_0.25_0.09_100
DimPlot(humaninvit_data, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste("~/Desktop/UMAP_ME_r2_doub",".pdf",sep=""),width = 10, height = 10)
D18 <- WhichCells(humaninvit_data,idents="Doublet")
saveRDS(D18,"~/Desktop/Douoblet_ME_r2.rds")

### Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
#seu_kidney <- CreateSeuratObject(kidney.data)
#seu_kidney <- NormalizeData(seu_kidney)
#seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
#seu_kidney <- ScaleData(seu_kidney)
#seu_kidney <- RunPCA(seu_kidney)
#seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)

### Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
#seu_kidney <- CreateSeuratObject(kidney.data)
#seu_kidney <- SCTransform(seu_kidney)
#seu_kidney <- RunPCA(seu_kidney)
#seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)

### pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
#sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
#bcmvn_kidney <- find.pK(sweep.stats_kidney)

### pK Identification (ground-truth) ------------------------------------------------------------------------------------------
#sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
#gt.calls <- seu_kidney@meta.data[rownames(sweep.res.list_kidney[[1]]), "GT"].   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
#sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = TRUE, GT.calls = gt.calls)
#bcmvn_kidney <- find.pK(sweep.stats_kidney)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
#omotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
#nExp_poi <- round(0.075*nrow(seu_kidney@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
#seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)



t1 <- c(0.83373494,	0.161445783,	0,	0.004819277,	0)
t2 <- c(0.904761905,	0.095238095,	0,	0,	0)
t3 <- c(0.997222222,	0.002777778,	0,	0,	0)
t4 <- c(0.005208333,	0.633680556,	0.147569444,	0.213541667,	0)
t5 <- c(0.109170306,	0.598253275,	0.096069869,	0.021834061,	0.174672489)
t6 <- c(0,	0.574608409,	0.069249794,	0.356141797,	0)

ids <- c("PS","Mes","Am","PGC","ES") 
tim1 <- c(18,18,18,18,18)
tim2 <- c(96,96,96,96,96)

C <- data.frame(x=c(ids,ids,ids,ids,ids,ids),y=c(t1,t2,t3,t4,t5,t6),z=c(tim1,tim1,tim1,tim2,tim2,tim2),z2 = c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6) )

ggplot(C, aes(fill=x, y=y, x=z2)) + 
  geom_bar(position="stack", stat="identity")

