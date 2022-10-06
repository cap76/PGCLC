
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)

set.seed(1)

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "~/Desktop/Thorsten/FINAL/AlignAraData_to_Marmoset/"
fileID <- saveext
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

fileID2 <- "/Volumes/Overflow/Uber18G/"
mammal.combined <- readRDS(paste(fileID2,'mammal.combined.KO5.rds',sep=""))

list1 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_4i.rds") 
list2 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_12h_r2.rds") 
list3 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_18h.rds") 
list4 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_24h.rds") 
list5 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_32h.rds") 
list6 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_40h.rds") 
list7 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_48h.rds") 
list8 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_D4_r2.rds") 
list9 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_D4.rds") 
list10 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_DE.rds") 
list11 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_hESC_W5_E8.rds") 
list12 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_ME_r2.rds") 
list13 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_PGC_wk6.5_invivo.rds")
list14 <- readRDS("/Volumes/GoogleDrive/My\ Drive/Chris\ and\ Ara\'s\ shared\ folder/2022-01-21\ Paper/Doublets/Run2/Douoblet_PreME.rds")


Idents(mammal.combined) <- mammal.combined$species
subuset1 <- subset(mammal.combined,idents = "Human (Ara EBs)")

allcells <- c(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,list13,list14)

newid <- as.character(Idents(subuset1))
newid[1:length(newid)] <- "Singlet"
Idents(subuset1) <- newid
Idents(subuset1, cells = str_replace(allcells,'-1','-7') ) <- "Doublet"

DimPlot(subuset1, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste("~/Desktop/Doublet",".pdf",sep=""),width = 15, height = 15,limitsize = FALSE)


