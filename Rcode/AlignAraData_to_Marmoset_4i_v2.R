library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(enrichR)

set.seed(1)

#Tediously add all the folders can shortcut this later when we have finalised everything.
saveext = "~/Desktop/Thorsten/FINAL/AlignAraData_to_Marmoset_w4i/"
fileID <- saveext
dir.create(saveext)
dir.create(paste(saveext,"/Markers/",sep=""))
dir.create(paste(saveext,"/DimRed/",sep=""))

fileID2 <- "/Volumes/Overflow/Uber18G_FINAL2/"
mammal.combined <- readRDS(paste(fileID2,'mammal.combined.KO3.rds',sep=""))

mammal.combined$Cells <- Idents(mammal.combined)

mammal.combined <- FindClusters(mammal.combined, resolution = 0.9)
mammal.combined$Cl9 <- Idents(mammal.combined)

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
Species2[which(is.na(mammal.combined$species)==1)] <- "3) Marmoset"
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

Data2 <- subset(mammal.combined, idents=c("1) Human EBs (in vitro)","7) Marmoset (in vivo)"))

Idents(Data2) <- Data2$ID2
Data2 <- subset(Data2, idents=c("DE","PGC","ME"),invert=TRUE)

Data2 <- FindNeighbors(Data2, reduction = "pca", dims = 1:20, k.param = 100)



#uID <- Data2$Cl9
#uID[which()]
Idents(Data2) <- Data2$Cl9
Cl8 <- WhichCells(Data2,idents=8)
SNN_K100 <- Data2@graphs$integrated_snn
othercells <- sort(SNN_K100[which(colnames(Data2)==Cl8[1]),which(Data2$species2=="7) Marmoset (in vivo)")],decreasing = TRUE)[1:5]

Idents(Data2) <- Data2$species2
Idents(Data2,cells=Cl8[1]) <- "B"
Idents(Data2,cells=names(othercells)) <- "A"
DimPlot(Data2, reduction = "umap", label = TRUE, repel = TRUE)# + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/test",".pdf",sep=""),width = 15, height = 15, useDingbats=FALSE)

#othercells <- sort(SNN_K100[which(colnames(Data2)%in%Cl8),which(Data2$species2=="7) Marmoset (in vivo)")],decreasing = TRUE)[1:5]




KNN_K100 <- Data2@graphs$integrated_nn
SNN_K100 <- Data2@graphs$integrated_snn

Data2 <- FindNeighbors(Data2, reduction = "pca", dims = 1:20, k.param = 50)

KNN_K50 <- Data2@graphs$integrated_nn
SNN_K50 <- Data2@graphs$integrated_snn

Data2 <- FindNeighbors(Data2, reduction = "pca", dims = 1:20, k.param = 150)

KNN_K150 <- Data2@graphs$integrated_nn
SNN_K150 <- Data2@graphs$integrated_snn

saveRDS(KNN_K100,file=paste(saveext,"KNN_K100.rds",sep=""))
saveRDS(SNN_K100,file=paste(saveext,"SNN_K100.rds",sep=""))

saveRDS(KNN_K150,file=paste(saveext,"KNN_K150.rds",sep=""))
saveRDS(SNN_K150,file=paste(saveext,"SNN_K150.rds",sep=""))

saveRDS(KNN_K50,file=paste(saveext,"KNN_K50.rds",sep=""))
saveRDS(SNN_K50,file=paste(saveext,"SNN_K50.rds",sep=""))

#uID <- as.character(Data2$Cl9s)
#uID[which(Data2$species2=="7) Marmoset (in vivo)")] <- as.character(Data2$ID2[which(Data2$species2=="7) Marmoset (in vivo)")])

uID <- as.character(Data2$Cl9)
uID[which(Data2$species2=="7) Marmoset (in vivo)")] <- as.character(Data2$ID2[which(Data2$species2=="7) Marmoset (in vivo)")])

#Cluster IDs
#2 = base
#10 = pre Am
#9=early pgc
#4 = pgc
#13=mid amnion
#5/8/7=late amnion
#6=early mes
#0,30,27=late mesoderm
#29/17/20 endo

uID[which(Data2$ID2=="ESC" & Data2$species2=="1) Human EBs (in vitro)")] <- "ESC"
uID[which(Data2$ID2=="PreME" & Data2$species2=="1) Human EBs (in vitro)")] <- "PreME"
uID[which(Data2$Cl9=="5" & Data2$species2=="1) Human EBs (in vitro)")] <- "Basal"

uID[which(Data2$Cl9=="10" & Data2$species2=="1) Human EBs (in vitro)")] <- "PreAMLC"

uID[which(Data2$Cl9=="8" & Data2$species2=="1) Human EBs (in vitro)")] <- "EarlyPGCLC"

uID[which(Data2$Cl9=="7" & Data2$species2=="1) Human EBs (in vitro)")] <- "PGCLC"

uID[which(Data2$Cl9=="11" & Data2$species2=="1) Human EBs (in vitro)")] <- "EarlyAMLC"

uID[which(Data2$Cl9=="9" & Data2$species2=="1) Human EBs (in vitro)")] <- "LateAMLC"
uID[which(Data2$Cl9=="3" & Data2$species2=="1) Human EBs (in vitro)")] <- "LateAMLC"
uID[which(Data2$Cl9=="6" & Data2$species2=="1) Human EBs (in vitro)")] <- "LateAMLC"


uID[which(Data2$Cl9=="4" & Data2$species2=="1) Human EBs (in vitro)")] <- "EarlyMELC"
uID[which(Data2$Cl9=="0" & Data2$species2=="1) Human EBs (in vitro)")] <- "LateMELC"

uID[which(Data2$Cl9=="14" & Data2$species2=="1) Human EBs (in vitro)")] <- "DELC"



Idents(Data2) <- uID
DimPlot(Data2, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_split",".pdf",sep=""),width = 30, height = 15)


Idents(Data2) <- Data2$Cl9
DimPlot(Data2, reduction = "umap", split.by = "species2", label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl9_split",".pdf",sep=""),width = 30, height = 15)
#Data2 <- FindClusters(Data2, resolution = 0.9)

Idents(Data2) <- uID

L1 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L1b <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
P1 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L2 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L2b <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
P2 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L3 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
L3b <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))
P3 <- matrix(, ncol = length(table(uID)), nrow = length(Idents(Data2)))

Lb <- table(Idents(Data2))
nn <- length(uID)

for (i in 1:length(Idents(Data2)) ) {
  L1[i,] <- table(Idents(Data2)[which(SNN_K50[i,]>0)])
  L1b[i,] <- L1[i,] / Lb
  P1[i,] <- phyper(L1[i,]-1, Lb, nn-Lb, 50, log = TRUE)
  
  L2[i,] <- table(Idents(Data2)[which(SNN_K100[i,]>0)])
  L2b[i,] <- 0.5*(L2b[i,] / Lb)
  P2[i,] <- phyper(L2[i,]-1, Lb, nn-Lb, 100, log = TRUE)

  L3[i,] <- table(Idents(Data2)[which(SNN_K150[i,]>0)])
  L3b[i,] <- (L3b[i,] / Lb)/3
  P3[i,] <- phyper(L3[i,]-1, Lb, nn-Lb, 150, log = TRUE)
}

colnames(L1)  <- rownames(Lb)
colnames(L1b) <- rownames(Lb)
colnames(P1)  <- rownames(Lb)
colnames(L2)  <- rownames(Lb)
colnames(L2b) <- rownames(Lb)
colnames(P2)  <- rownames(Lb)
colnames(L3)  <- rownames(Lb)
colnames(L3b) <- rownames(Lb)
colnames(P3)  <- rownames(Lb)
rownames(L1) <- colnames(Data2)
rownames(L2) <- colnames(Data2)
rownames(L3) <- colnames(Data2)
rownames(L1b) <- colnames(Data2)
rownames(L2b) <- colnames(Data2)
rownames(L3b) <- colnames(Data2)
rownames(P1) <- colnames(Data2)
rownames(P2) <- colnames(Data2)
rownames(P3) <- colnames(Data2)

saveRDS(L1,paste(fileID,'L1.rds',sep=""))
saveRDS(L1b,paste(fileID,'L1b.rds',sep=""))
saveRDS(P1,paste(fileID,'P1.rds',sep=""))
saveRDS(L2,paste(fileID,'L2.rds',sep=""))
saveRDS(L2b,paste(fileID,'L2b.rds',sep=""))
saveRDS(P2,paste(fileID,'P2.rds',sep=""))
saveRDS(L3,paste(fileID,'L3.rds',sep=""))
saveRDS(L3b,paste(fileID,'L3b.rds',sep=""))
saveRDS(P3,paste(fileID,'P3.rds',sep=""))

write.table(L1,file=paste(fileID,'L1.csv',sep=""),sep=",")
write.table(L1b,file=paste(fileID,'L1b.csv',sep=""),sep=",")
write.table(P1,file=paste(fileID,'P1.csv',sep=""),sep=",")
write.table(L2,file=paste(fileID,'L2.csv',sep=""),sep=",")
write.table(L2b,file=paste(fileID,'L2b.csv',sep=""),sep=",")
write.table(P2,file=paste(fileID,'P2.csv',sep=""),sep=",")
write.table(L3,file=paste(fileID,'L3.csv',sep=""),sep=",")
write.table(L3b,file=paste(fileID,'L3b.csv',sep=""),sep=",")
write.table(P3,file=paste(fileID,'P3.csv',sep=""),sep=",")

inds <- which(Data2$species2=="7) Marmoset (in vivo)")

write.table(Lb,file=paste(fileID,'Lbase.csv',sep=""),sep=",")


write.table(L1[inds,],file=paste(fileID,'L1_.csv',sep=""),sep=",")
write.table(L1b[inds,],file=paste(fileID,'L1b_.csv',sep=""),sep=",")
write.table(P1[inds,],file=paste(fileID,'P1_.csv',sep=""),sep=",")
write.table(L2[inds,],file=paste(fileID,'L2_.csv',sep=""),sep=",")
write.table(L2b[inds,],file=paste(fileID,'L2b_.csv',sep=""),sep=",")
write.table(P2[inds,],file=paste(fileID,'P2_.csv',sep=""),sep=",")
write.table(L3[inds,],file=paste(fileID,'L3_.csv',sep=""),sep=",")
write.table(L3b[inds,],file=paste(fileID,'L3b_.csv',sep=""),sep=",")
write.table(P3[inds,],file=paste(fileID,'P3_.csv',sep=""),sep=",")




inds1 <- which(Data2$species2=="1) Human EBs (in vitro)")
inds2 <- which(Data2$species2=="7) Marmoset (in vivo)")
k = 5
C1 <- matrix(, ncol = k, nrow = length(Idents(Data2)))
S1 <- matrix(, ncol = k, nrow = length(Idents(Data2)))
for (i in 1:length(Idents(Data2)) ) {
  a1 <- sort(SNN_K50[i,inds2],decreasing = TRUE)[1:k]
  C1[i,] <- names(a1)
  S1[i,] <- a1
}
write.table(Idents(Data2)[inds1],file=paste(fileID,'AllMapIdents.csv',sep=""),sep=",")
write.table(C1[inds1,],file=paste(fileID,'C1_.csv',sep=""),sep=",")
write.table(S1[inds1,],file=paste(fileID,'S1.csv',sep=""),sep=",")



k = 10
C1 <- matrix(, ncol = k, nrow = length(Idents(Data2)))
S1 <- matrix(, ncol = k, nrow = length(Idents(Data2)))
for (i in 1:length(Idents(Data2)) ) {
  a1 <- sort(SNN_K100[i,inds2],decreasing = TRUE)[1:k]
  C1[i,] <- names(a1)
  S1[i,] <- a1
}
write.table(C1[inds1,],file=paste(fileID,'C1_100.csv',sep=""),sep=",")
write.table(S1[inds1,],file=paste(fileID,'S1_100.csv',sep=""),sep=",")



#Load in the datasets and get the intersecting genes that are in the gene models
raw_counts1<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/cynomolgous/featurecountsNakamura3.csv",sep=",",header = T, row.names=1)
raw_countsA1<-read.table("/Users/christopherpenfold/Desktop/Aracely/featurecountsLi.csv",sep=",",header = T, row.names=1)
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/featurecountsAll_CAPProcessed.csv",sep=",",header = T, row.names=1)
#raw_counts4 <- readRDS("/Volumes/GoogleDrive/Shared\ drives/MALKOWSKA-update/5\ species\ metabolism\ paper\ /Datasets/Human/Human_CS7.rds")
raw_counts4 <- readRDS("/Users/christopherpenfold/Downloads/expression_values_1.rds")
#raw_counts4 <- readRDS("/Users/christopherpenfold/Downloads/expression_values_1.rds")
colnames(raw_counts4)[which(colnames(raw_counts4)=="TBXT")] <- "T"
rownames(raw_counts1)[which(rownames(raw_counts1)=="ENSMFAG00000033399")] = "NANOG"

#ENSCJAG00000040521==GATA1

#Look at these that map in all experiments
common_genes1 <- intersect(rownames(raw_counts1),rownames(raw_counts3))
common_genes2 <- intersect(rownames(raw_countsA1),colnames(raw_counts4))
common_genes <- intersect(common_genes1,common_genes2)

X = GetAssayData(Data2, assay = "RNA")[common_genes,inds1]
write.table(as.data.frame(X),file=paste(fileID,'Expression.csv',sep=""),sep=",")




#Other ideas. KNNs for each cell in a set, then minimise distance (e.g., ordering)
sadasdasdasdas


Inds1 <- which(mammal.combined$species=="3) Marmoset"  )
Den2 <- matrix(, nrow = length(Inds1), ncol = length((as.numeric(unique(Idents(mammal.combined))))))
Den3 <- matrix(, nrow = length(Inds1), ncol = length((as.numeric(unique(Idents(mammal.combined))))))
#Den4 <- matrix(, nrow = length(Inds1), ncol = length((as.numeric(unique(Idents(mammal.combined))))))

L1 <- matrix(, nrow = 1, ncol = length((as.numeric(unique(Idents(mammal.combined))))))

for (i in seq(1,length(levels(Idents(mammal.combined))),by = 1) ) {
    Inds2 <- which(mammal.combined$species=="human_synth" & Idents(mammal.combined)==levels(Idents(mammal.combined))[i] )    
    #Den1[,i] <- rowSums(as.data.frame(KNN50_K20[Inds1,Inds2]))
    Den2[,i] <- rowSums(as.data.frame(KNN50_K100[Inds1,Inds2]))
    Den3[,i] <- rowSums(as.data.frame(KNN50_K500[Inds1,Inds2]))
    #Den4[,i] <- rowSums(as.data.frame(KNN50_K1000[Inds1,Inds2]))
    L1[i] <- length(Inds2)
}

#write.csv(as.data.frame(Den1), file=paste(saveext,"/Den1.csv",sep=""))
write.csv(as.data.frame(Den2), file=paste(saveext,"/Den2.csv",sep=""))
write.csv(as.data.frame(Den3), file=paste(saveext,"/Den3.csv",sep=""))
#write.csv(as.data.frame(Den4), file=paste(saveext,"/Den4.csv",sep=""))
write.csv(as.data.frame(L1), file=paste(saveext,"/ClLen.csv",sep=""))
write.csv(as.data.frame(colnames(mammal.combined)[Inds1]), file=paste(saveext,"/EBIDs.csv",sep=""))



#Cluster summaries

#PreME 40,30,10,31,32,1

#18 DE plpo
#11,21,34 DE
#37,2 ESc

#8,15 12h
 #8 MELC bias

#5,7,28 18h
  #7,pgc

#3,6,16,14,19
  #16 sox17
  #3,6 are MELC
  #14,9?

#4,27,29,35 MELC2
#13,26,24,22,23 EPGCLC
#0 LPGCLC

#17,20,12 AMLC

#33 bridge AMLC PGCLC

#7,36,25,9 MELC1

library("pheatmap")
Idents(mammal.combined) <- mammal.combined$species
subs <- subset(mammal.combined,idents=c("3) Marmoset","human_synth"))

Idents(subs) <- colnames(subs)
avexp  <- AverageExpression(object = subs, slot = "data")


DN1 <- GetAssayData(subs, assay = "integrated")
C1 <- cor(as.matrix(DN1[,which(subs$species=="3) Marmoset")]),as.matrix(DN1[,which(subs$species=="Human (Ara EBs)")]), method = "pearson")
write.csv(as.data.frame( C1 ), file=paste(saveext,"/CorrInt.csv",sep=""))

#d1 <- as.matrix(D[which(subs$species=="3) Marmoset"),])
#d2 <- as.matrix(D[which(subs$species=="Human (Ara EBs)"),])
#as.matrix(D[,which(subs$species=="3) Marmoset")]),as.matrix(D[,which(subs$species=="Human (Ara EBs)")])
D <- as.data.frame(subs[["umap"]]@cell.embeddings)
DM <- as.matrix( dist(D) )
DM2 <- DM[which(subs$species=="3) Marmoset"), which(subs$species=="Human (Ara EBs)")]
write.csv(as.data.frame( DM2 ), file=paste(saveext,"/DistInt.csv",sep=""))

cell1 <- subs$Cells[which(subs$species=="3) Marmoset")]
cell2 <- subs$Cells[which(subs$species=="Human (Ara EBs)")]

cl1 <- subs$seurat_clusters[which(subs$species=="3) Marmoset")]
cl2 <- subs$seurat_clusters[which(subs$species=="Human (Ara EBs)")]

write.csv(as.data.frame( cl1 ), file=paste(saveext,"/Cl1.csv",sep=""))
write.csv(as.data.frame( cl2 ), file=paste(saveext,"/Cl2.csv",sep=""))


write.csv(as.data.frame( cell1 ), file=paste(saveext,"/Cell1.csv",sep=""))
write.csv(as.data.frame( cell2 ), file=paste(saveext,"/Cell2.csv",sep=""))

newVector1 <- vector(mode = "numeric", length = dim(C1)[1])
newVector2 <- vector(mode = "numeric", length = dim(C1)[2])

s1 <- seq(from = 1, to = dim(C1)[2], by = 1)
s2 <- seq(from = 1, to = dim(C1)[1], by = 1)

for (i in s2) {
  #s2[i] <- which(C1[,i]==max(C1[,i]))
  s2[i] <- which(DM2[,i]==min(DM2[,i]))  
}

#Now scatter plots
ID <- Idents(mammal.combined)
species2 <- mammal.combined$species
species2[which(species2=="human_synth")] <- "Amnioid"
species2[which(species2=="Human (Ara EBs)")] <- "EB"
species2[which(species2=="3) Marmoset")] <- "Marmoset"


#Data2 <- FindClusters(Data2, resolution = 0.5)
Idents(Data2) <- Data2$Cl5
DimPlot(Data2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl5_split",".pdf",sep=""),width = 30, height = 15)
#Data2 <- FindClusters(Data2, resolution = 0.9)

Idents(Data2) <- Data2$Cl9
DimPlot(Data2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) 
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl9_split",".pdf",sep=""),width = 30, height = 15)
#Data2 <- FindClusters(Data2, resolution = 1.5)
#DimPlot(Data2, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + NoLegend()
#ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl15_split",".pdf",sep=""),width = 30, height = 15)







mammal.combined <- FindClusters(mammal.combined, resolution = 0.5)
mammal.combined$Cl5 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl5_gast2",".pdf",sep=""),width = 32, height = 10)

mammal.combined <- FindClusters(mammal.combined, resolution = 1)
mammal.combined$Cl10 <- Idents(mammal.combined)
DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl10_gast2",".pdf",sep=""),width = 32, height = 10)


NID <- paste(species2,"_",mammal.combined$Cl10,sep="")
Idents(mammal.combined) <- NID



#ID2 <- as.character(ID)
#ID2[which(species2=="EB")] <- as.character(mammal.combined$Cl5[which(species2=="EB")])

#ID2[which(species2=="CS7" & mammal.combined$Cl5==11)] <- "PGC"

#14,2,13,12,7

#Pull global markers from Ara's data
DefaultAssay(mammal.combined) <- "RNA"

cluster1 <- FindMarkers(mammal.combined, ident.1 = "EB_14", ident.2 = c("EB_2","EB_13","EB_12","EB_7","EB_19","EB_9","EB_0","EB_5"), min.pct = 0.1, only.pos = TRUE, test.use = "MAST")
cluster2 <- FindMarkers(mammal.combined, ident.1 = "EB_2", ident.2 = c("EB_14","EB_13","EB_12","EB_7","EB_19","EB_9","EB_0","EB_5"), min.pct = 0.1, only.pos = TRUE, test.use = "MAST")
cluster3 <- FindMarkers(mammal.combined, ident.1 = c("EB_12","EB_13"), ident.2 = c("EB_2","EB_14","EB_7","EB_19","EB_9","EB_0","EB_5"), min.pct = 0.1, only.pos = TRUE, test.use = "MAST")
cluster5 <- FindMarkers(mammal.combined, ident.1 = "EB_7", ident.2 = c("EB_2","EB_13","EB_12","EB_14","EB_19","EB_9","EB_0","EB_5"), min.pct = 0.1, only.pos = TRUE, test.use = "MAST")
cluster6 <- FindMarkers(mammal.combined, ident.1 = "EB_19", ident.2 = c("EB_2","EB_13","EB_12","EB_14","EB_7","EB_9","EB_0","EB_5"), min.pct = 0.1, only.pos = TRUE, test.use = "MAST")
cluster7 <- FindMarkers(mammal.combined, ident.1 = "EB_9", ident.2 = c("EB_2","EB_13","EB_12","EB_14","EB_7","EB_19","EB_0","EB_5"), min.pct = 0.1, only.pos = TRUE, test.use = "MAST")
cluster8 <- FindMarkers(mammal.combined, ident.1 = "EB_0", ident.2 = c("EB_2","EB_13","EB_12","EB_14","EB_7","EB_19","EB_9","EB_5"), min.pct = 0.1, only.pos = TRUE, test.use = "MAST")
cluster9 <- FindMarkers(mammal.combined, ident.1 = "EB_5", ident.2 = c("EB_2","EB_13","EB_12","EB_14","EB_7","EB_19","EB_9","EB_0"), min.pct = 0.1, only.pos = TRUE, test.use = "MAST")

MarkerList <- c(rownames(cluster1),rownames(cluster2),rownames(cluster3),rownames(cluster4),rownames(cluster5),rownames(cluster6),rownames(cluster7),rownames(cluster8),rownames(cluster9))


#Idents(mammal.combined) <- NID # paste(species2,"_",ID2)

DimPlot(mammal.combined, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE)
ggsave(filename=paste(saveext,"/DimRed/UMAP_split_Cl10_gast2wtf",".pdf",sep=""),width = 32, height = 10)

NID2 <- NID
NID2[which(species2=="Marmoset")] <- colnames(mammal.combined)[which(species2=="Marmoset")]

Idents(mammal.combined) <- NID2 # paste(species2,"_",ID2)

DefaultAssay(mammal.combined) <- "integrated"
avg.Cl <- AverageExpression(mammal.combined, verbose = FALSE)
avg.Cl$gene <- rownames(avg.Cl)

Ind1 <- c(grep("ERR", colnames(avg.Cl$integrated)),
grep("SLX", colnames(avg.Cl$integrated)),
grep("CME", colnames(avg.Cl$integrated)),
grep("X35", colnames(avg.Cl$integrated)))

Ind2 <- grep("EB", colnames(avg.Cl$integrated))

b <- avg.Cl$integrated[,Ind2]
c <- avg.Cl$integrated[,Ind1]
C1 <- cor(log2(b+1),log2(c+1), method = "pearson")

library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.3, 0.7, length.out = 20)
pheatmap(C1,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sampler_corr_gast",".pdf",sep=""),width=75,height=10)



b <- avg.Cl$RNA[rownames(avg.Cl$integrated),Ind2]
c <- avg.Cl$RNA[rownames(avg.Cl$integrated),Ind1]
C2 <- cor(log2(b+1),log2(c+1), method = "pearson")

library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.3, 0.7, length.out = 20)
pheatmap(C2,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sampler_corr_gast2",".pdf",sep=""),width=75,height=10)


C3 <- cor(log2(b+1),log2(c+1), method = "spearman")

library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0.3, 0.7, length.out = 20)
pheatmap(C3,color =  redblue1(20),breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sampler_corr_gast3",".pdf",sep=""),width=75,height=10)





b <- avg.Cl$RNA[MarkerList,Ind2]
c <- avg.Cl$RNA[MarkerList,Ind1]



C4 <- cor(log2(b+1),log2(c+1), method = "pearson")

library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0, 0.5, length.out = 20)
pheatmap(C4,color =  redblue1(20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sampler_corr_gast4",".pdf",sep=""),width=75,height=10)


C5 <- cor(log2(b+1),log2(c+1), method = "spearman")

library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0, 0.5, length.out = 20)
pheatmap(C5,color =  redblue1(20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sampler_corr_gast5",".pdf",sep=""),width=75,height=10)



c2 <- c
c2[c2!=0] <- 1

MarkerList2 <- MarkerList[which(rowSums(c2) > 5)]


b <- avg.Cl$RNA[MarkerList2,Ind2]
c <- avg.Cl$RNA[MarkerList2,Ind1]

C6 <- cor(log2(b+1),log2(c+1), method = "pearson")

library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0, 0.6, length.out = 20)
pheatmap(C6,color =  redblue1(20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sampler_corr_gast6",".pdf",sep=""),width=75,height=10)


C7 <- cor(log2(b+1),log2(c+1), method = "spearman")

library(pheatmap)
redblue1<-colorRampPalette(c("#245199","#FFFFFF","#D20000"))
mat_breaks <- seq(0, 0.6, length.out = 20)
pheatmap(C7,color =  redblue1(20), breaks = mat_breaks, border_color = NA, cluster_rows=FALSE,cluster_cols=FALSE, filename = paste(saveext,"/DimRed/Sampler_corr_gast7",".pdf",sep=""),width=75,height=10)




write.csv(as.data.frame( C1 ), file=paste(saveext,"/NormData.csv",sep=""))
write.csv(as.data.frame( C2 ), file=paste(saveext,"/RNAData.csv",sep=""))
write.csv(as.data.frame( C3 ), file=paste(saveext,"/RNAData_spear.csv",sep=""))

write.csv(as.data.frame( C4 ), file=paste(saveext,"/RNAData_marker.csv",sep=""))
write.csv(as.data.frame( C5 ), file=paste(saveext,"/RNAData_spear_marker.csv",sep=""))

write.csv(as.data.frame( C6 ), file=paste(saveext,"/RNAData_marker.csv",sep=""))
write.csv(as.data.frame( C7 ), file=paste(saveext,"/RNAData_spear_marker.csv",sep=""))




BS<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/Keycorrect_CPAll2.csv",sep=",",header = T, row.names=1)
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/featurecountsAll_CAPProcessed.csv",sep=",",header = T, row.names=1)
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
isinvit <- BS$All * BS$QC
labs <- BS$Type
labs2 <- BS$Annotation2
labs <- labs[which(isinvit>0)]
labs2 <- labs2[which(isinvit>0)]
raw_counts3<-read.table("/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/featurecountsAll_CAPProcessed.csv",sep=",",header = T, row.names=1)
raw_counts3<- raw_counts3[,which(isinvit>0)]
marmoset_data <- CreateSeuratObject(counts = raw_counts3, assay = "RNA",min.cells = 0, min.features = 0)
marmoset_data$species <- "3) Marmoset"
marmoset_data$ID1 <- labs2 #"newworld"
marmoset_data$ID2 <- labs2 #"primate"
marmoset_data <- subset(marmoset_data, subset = nFeature_RNA > 0)
Idents(marmoset_data) <- labs2
marmoset_data <- subset(marmoset_data, idents = c("EmDisc_CS5","EmDisc_CS6","EmDisc_CS7","PGC_CS5","PGC_CS6"), invert = FALSE)
#marmoset_data <- subset(marmoset_data, idents = c("Myo_CS7","Tb_CS3","Tb_CS5","Tb_CS6","Tb_CS7","ReGland_CS5","ReGland_CS7","Gland_CS5","Gland_CS6","Gland_CS7","4-cell_CS2","8-cell_CS2","Zy_CS1","cMor_CS3"), invert = TRUE)
marmoset_data <- NormalizeData(marmoset_data, verbose = FALSE)
marmoset_data <- FindVariableFeatures(marmoset_data, selection.method = "vst", nfeatures = 20000)

marmoset_data <- ScaleData(marmoset_data, verbose = FALSE)
marmoset_data <- RunPCA(marmoset_data, npcs = 30, verbose = FALSE)
marmoset_data <- RunUMAP(marmoset_data, reduction = "pca", dims = 1:20)
marmoset_data <- RunTSNE(marmoset_data, reduction = "pca", dims = 1:20)
marmoset_data <- FindNeighbors(marmoset_data, reduction = "pca", dims = 1:20)
#mammal.combined$Cells <- Idents(mammal.combined)

DimPlot(marmoset_data, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/MARMPCA_split_k100",".pdf",sep=""),width = 32, height = 10)
DimPlot(marmoset_data, reduction = "umap", split.by = "species", label = TRUE, repel = TRUE,pt.size = 6) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/MARMUMAP_split_k100",".pdf",sep=""),width = 32, height = 10)

FeaturePlot(marmoset_data,  reduction = "pca",features = "SOX17", split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 6, min.cutoff = 0)
ggsave(filename=paste(saveext,"/Markers/MARMMarker_SOX17_k100.pdf",sep=""),width = 32, height = 10,limitsize = FALSE)
FeaturePlot(marmoset_data,  reduction = "pca",features = "NANOS3", split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 6, min.cutoff = 0)
ggsave(filename=paste(saveext,"/Markers/MARMMarker_NANOS3_k100.pdf",sep=""),width = 32, height = 10,limitsize = FALSE)
FeaturePlot(marmoset_data,  reduction = "pca",features = "VTCN1", split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 6, min.cutoff = 0)
ggsave(filename=paste(saveext,"/Markers/MARMMarker_VTCN1_k100.pdf",sep=""),width = 32, height = 10,limitsize = FALSE)
FeaturePlot(marmoset_data,  reduction = "pca",features = "TFAP2C", split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 6, min.cutoff = 0)
ggsave(filename=paste(saveext,"/Markers/MARMMarker_TFAP2C_k100.pdf",sep=""),width = 32, height = 10,limitsize = FALSE)
FeaturePlot(marmoset_data,  reduction = "pca",features = "TFAP2A", split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 6, min.cutoff = 0)
ggsave(filename=paste(saveext,"/Markers/MARMMarker_TFAP2A_k100.pdf",sep=""),width = 32, height = 10,limitsize = FALSE)
FeaturePlot(marmoset_data,  reduction = "pca",features = "MIXL1", split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 6, min.cutoff = 0)
ggsave(filename=paste(saveext,"/Markers/MARMMarker_MIXL1_k100.pdf",sep=""),width = 32, height = 10,limitsize = FALSE)
FeaturePlot(marmoset_data,  reduction = "pca",features = "T", split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 6, min.cutoff = 0)
ggsave(filename=paste(saveext,"/Markers/MARMMarker_T_k100.pdf",sep=""),width = 32, height = 10,limitsize = FALSE)
FeaturePlot(marmoset_data,  reduction = "pca",features = "POU5F1", split.by = "species", cols =  c("lightgrey", "darkblue"), pt.size = 6, min.cutoff = 0)
ggsave(filename=paste(saveext,"/Markers/MARMMarker_OCT4_k100.pdf",sep=""),width = 32, height = 10,limitsize = FALSE)





marmoset_data <- FindClusters(marmoset_data, resolution = 0.5)
DimPlot(mammal.combined, reduction = "pca", split.by = "species", label = TRUE, repel = TRUE,pt.size = 2) #+ NoLegend()
ggsave(filename=paste(saveext,"/DimRed/MarmPCA_0p5_Cl1",".pdf",sep=""),width = 32, height = 10)