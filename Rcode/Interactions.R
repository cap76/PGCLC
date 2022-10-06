library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)

set.seed(1)

fileID2 <- "/Volumes/Overflow/Uber18G/"
mammal.combined <- readRDS(paste(fileID2,'mammal.combined.KO5.rds',sep=""))

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "~/Desktop/Thorsten/FINAL/AraInteraction/"
fileID <- saveext
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

labs <- as.character( mammal.combined$Cells )
labs[which(labs=="MeLC_R1")] <- "PreME"
labs[which(labs=="MeLC_R2")] <- "PreME"
labs[which(labs=="MeLC_R3")] <- "PreME"
labs[which(labs=="MeLC_R4")] <- "PreME"
labs[which(labs=="MELC_KO")] <- "PreME"
labs[which(labs=="ESC_R1")] <- "ESC"
labs[which(labs=="ESC_R2")] <- "ESC"
labs[which(labs=="ESC_R3")] <- "ESC"
labs[which(labs=="ESC_R4")] <- "ESC"
labs[which(labs=="ESC_KO")] <- "ESC"
labs[which(labs=="PGCLC_D1_R1")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D1_R2")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D1_R3")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D1_R4")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D1_KO")] <- "PGCLC_24h"
labs[which(labs=="PGCLC_D2_R1")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D2_R2")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D2_R3")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D2_R4")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D2_KO")] <- "PGCLC_48h"
labs[which(labs=="PGCLC_D3_R1")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D3_R2")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D3_R3")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D3_R4")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D3_KO")] <- "PGCLC_72h"
labs[which(labs=="PGCLC_D4_R1")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4_R2")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4_R3")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4_R4")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4_KO")] <- "PGCLC_96h"
labs[which(labs=="hESC")] <- "ESC"
labs[which(labs=="hPGCLC")] <- "PGCLC"
labs[which(labs=="twAMLC")] <- "AMLC"
labs[which(labs=="MELC1")] <- "MELC1"
labs[which(labs=="MELC2")] <- "MELC2"
labs[which(labs=="PreME_10h")] <- "PreME"
labs[which(labs=="hESC_W5_E8")] <- "ESC"
labs[which(labs=="PGCLC_D4_r2")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_D4")] <- "PGCLC_96h"
labs[which(labs=="PGCLC_12h_r2")] <- "PGCLC_12h"
labs[which(labs=="TE")] <- "Tb"
labs[which(labs=="PGC_wk6.5_invivo")] <- "PGC"
labs[which(labs=="ME_r2")] <- "ME"

mammal.combined$ID2 <- as.factor(labs)
Idents(mammal.combined) <- as.factor(labs)

#Rejig IDs and order
Species2 <- as.character(mammal.combined$species)
Species2[which(Species2=="Human (Ara EBs)")] <- "1) Human EBs (in vitro)"
Species2[which(Species2=="Human (Amnioid)")] <- "2) Human Amnioids (in vitro)"
Species2[which(Species2=="Human (CS7)")] <- "5) Human CS7 (in vivo)"
Species2[which(Species2=="4) Human (in vitro)")] <- "4) Human (in vitro 1)"
Species2[which(Species2=="human_invivo")] <- "5) Human (in vitro 2)"
Species2[which(Species2=="3) Cynomologous (in vitro)")] <- "6) Cynomolgus (in vitro)"
Species2[which(Species2=="3) Marmoset")] <- "7) Marmoset (in vivo)"
Species2[which(Species2=="10X_1")] <- "4) Human EBs (Clark1)"
Species2[which(Species2=="10X_2")] <- "4) Human EBs (Clark1)"
Species2[which(Species2=="10X_3")] <- "4) Human EBs (Clark2)"
Species2[which(Species2=="10X_4")] <- "4) Human EBs (Clark2)"
Species2[which(Species2=="10X_5")] <- "8) Human EBs (Clark)"
Species2[which(Species2=="Gastruloid2D")] <- "3) Gastruloid"

mammal.combined$species2 <- Species2
Idents(mammal.combined) <- mammal.combined$species2

subsetdata <- subset(mammal.combined,idents=c("1) Human EBs (in vitro)"))

#DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl5_gast2",".pdf",sep=""),width = 32, height = 10)



#id0 <- Id3
#Data2$Mapping <- id0
#Idents(Data2)<- Data2$ID2

Type <- as.character(subsetdata$Cl9)
Type[which(Type%in%c(6))] <- "eMELC"
Type[which(Type%in%c(0,30,27))] <- "aMELC"
Type[which(Type%in%c(2))] <- "PSLC"
Type[which(Type%in%c(10))] <- "pAm1"
Type[which(Type%in%c(13))] <- "eAm"
Type[which(Type%in%c(5,7,8,24))] <- "lAm"
Type[which(Type%in%c(17,20))] <- "DE"
Type[which(Type%in%c(15,16))] <- "lDE"
Type[which(Type%in%c(29,26))] <- "Other"
Type[which(Type%in%c(22))] <- "ME"
Type[which(Type%in%c(9))] <- "ePGCLC"
Type[which(Type%in%c(4))] <- "lPGCLC"

Type[which(Type%in%c(11,12))] <- "Pl_4i"
Type[which(Type%in%c(1))] <- "Pl_1"
Type[which(Type%in%c(3))] <- "Pl_3"
Type[which(Type%in%c(19))] <- "Pl_19"
Type[which(Type%in%c(18))] <- "Pl_18"
Type[which(Type%in%c(23))] <- "Pl_23"
Type[which(Type%in%c(14))] <- "Pl_14"
Type[which(Type%in%c(21))] <- "Pl_21"
Type[which(Type%in%c(25))] <- "Pl_25"
Type[which(Type%in%c(28,11,12))] <- "Pl_11"


Idents(subsetdata)<- subsetdata$Cells
DimPlot(subsetdata, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_check",".pdf",sep=""),width = 10, height = 10)
Idents(subsetdata)<- subsetdata$Cl9
DimPlot(subsetdata, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_clcheck",".pdf",sep=""),width = 10, height = 10)

Data2 <- subsetdata
Data2$uID <- Type
Idents(Data2) <- Data2$Cells


require(tidyverse)

Dat <- data.frame(y=Data2$uID,x=Data2$Cells)
Dat1 <- Dat %>%  count(y, x)
ggplot(Dat1, aes(fill=y, x=x, y=n)) + geom_bar(position="fill", stat="identity") + theme_bw()
ggsave(filename=paste(saveext,"AraTable","_fin.pdf",sep=""),width = 10, height = 10, limitsize = FALSE)



Data3 <- subset(Data2,idents="PGCLC_12h_r2")
Idents(Data3) <- Data3$uID 
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/12h/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/12h/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/12h/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/12h/barcodes.tsv",sep="") )
Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'d12_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)



Data3 <- subset(Data2,idents="PGCLC_18h")
Idents(Data3) <- Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/18h/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/18h/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/18h/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/18h/barcodes.tsv",sep="") )

Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'18h_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)



Data3 <- subset(Data2,idents="PGCLC_24h")
Idents(Data3) <-  Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/24h/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/24h/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/24h/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/24h/barcodes.tsv",sep="") )

Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'24h_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)
#DEGs <- FindAllMarkers(MARM,test.use = 'MAST',verbose = F,only.pos = T,random.seed = 1, logfc.threshold = 0.2, min.pct = 0.1,return.thresh = 0.05)
#fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_logFC > 0.1)
#fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_logFC', 'pct.1', 'pct.2')] 
#write.table(fDEGs, file ='E87_example_DEGs.tsv', sep = '\t', quote = F, row.names = F)




Data3 <- subset(Data2,idents="PGCLC_32h")
Idents(Data3) <-  Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/32h/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/32h/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/32h/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/32h/barcodes.tsv",sep="") )

Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'32h_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)
#DEGs <- FindAllMarkers(MARM,test.use = 'MAST',verbose = F,only.pos = T,random.seed = 1, logfc.threshold = 0.2, min.pct = 0.1,return.thresh = 0.05)
#fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_logFC > 0.1)
#fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_logFC', 'pct.1', 'pct.2')] 
#write.table(fDEGs, file ='E87_example_DEGs.tsv', sep = '\t', quote = F, row.names = F)






Data3 <- subset(Data2,idents="PGCLC_40h")
Idents(Data3) <- Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/40h/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/40h/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/40h/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/40h/barcodes.tsv",sep="") )

Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'40h_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)
#DEGs <- FindAllMarkers(MARM,test.use = 'MAST',verbose = F,only.pos = T,random.seed = 1, logfc.threshold = 0.2, min.pct = 0.1,return.thresh = 0.05)
#fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_logFC > 0.1)
#fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_logFC', 'pct.1', 'pct.2')] 
#write.table(fDEGs, file ='E87_example_DEGs.tsv', sep = '\t', quote = F, row.names = F)



Data3 <- subset(Data2,idents="PGCLC_48h")
Idents(Data3) <- Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/48h/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/48h/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/48h/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/48h/barcodes.tsv",sep="") )

Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'48h_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)



Data3 <- subset(Data2,idents="PGCLC_D4")
Idents(Data3) <- Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/96h/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/96h/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/96h/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/96h/barcodes.tsv",sep="") )

Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'96h_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)



Data3 <- subset(Data2,idents="hESC_W5_E8")
Idents(Data3) <- Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/ESC/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/ESC/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/ESC/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/ESC/barcodes.tsv",sep="") )

Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'ESC_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)


Data3 <- subset(Data2,idents="PreME_10h")
Idents(Data3) <- Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/PreME/",sep=""))
writeMM(Data3@assays$RNA@data, file = paste(saveext,'/PreME/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/PreME/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/PreME/barcodes.tsv",sep="") )
Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'PreME_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)


Data3 <- subset(Data2,idents="4i")
Idents(Data3) <- Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/4i/",sep=""))
writeMM(Data3@assays$RNA@data, file = paste(saveext,'/4i/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/4i/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/4i/barcodes.tsv",sep="") )
Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'4i_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)



Data3 <- subset(Data2,idents="PGCLC_D4_r2")
Idents(Data3) <- Data3$uID
#Now we can think about CellPhonedb
dir.create(paste(saveext,"/96h4i/",sep=""))

writeMM(Data3@assays$RNA@data, file = paste(saveext,'/96h4i/matrix.mtx',sep="") )
# save gene and cell names
write(x = rownames(Data3@assays$RNA@data), file = paste(saveext,"/96h4i/features.tsv",sep="") )
write(x = colnames(Data3@assays$RNA@data), file = paste(saveext,"/96h4i/barcodes.tsv",sep="") )

Data3$cell_type <-  Idents(Data3)
Data3@meta.data$Cell = rownames(Data3@meta.data)
df = Data3@meta.data[, c('Cell', 'cell_type')]
write.table(df, file =paste(saveext,'96h4i_example_meta.tsv',sep=""), sep = '\t', quote = F, row.names = F)

gnames <- rownames(GetAssayData(subset(Data2),assay="integrated"))


gnames <- rownames(GetAssayData(subset(humaninvit_data),assay="RNA"))

Idents(Data2) <- Data2$Cells
Data2 <- subset(Data2,idents=c("PGC_wk6.5_invivo"),invert=TRUE)
Idents(Data2) <- Data2$uID
Data3 <- GetAssayData(subset(Data2,idents="eMELC"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'eMLC_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="aMELC"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'aMLC_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="eAm"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'eAm_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="lAm"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'lAm_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="pAm1"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'pAm_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="ePGCLC"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'ePGCLC_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="lPGCLC"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'lPGCLC_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="PSLC"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'PSLC_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="DE"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'DE_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="Pl_1"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'Pl_1_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="Pl_3"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'Pl_3_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="Pl_18"),assay="RNA")[gnames,][gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'Pl_18_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="Other"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'Other_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="Pl_19"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'Pl_19_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="Pl_14"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'Pl_14_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)
Data3 <- GetAssayData(subset(Data2,idents="lDE"),assay="RNA")[gnames,]
write.table(as.data.frame(t(Data3) ), file =paste(saveext,'lDE_exp.csv',sep=""), sep = '\t', quote = F, row.names = F)

#PGCLC_18h         PGCLC_32h         PGCLC_40h         PGCLC_24h         PGCLC_12h         PGCLC_48h         PGCLC_96h   
#ESC               PreME 