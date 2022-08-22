######################## alternative splicing analysis ##################################

# Install and activate the required packages.
library("reshape")
library(DESeq2)
library(grid)
library(genefilter)
library(pheatmap)
library(VennDiagram)
library(ggplot2)
library("RColorBrewer")
library(EnhancedVolcano)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(yarrr)
library(plyr)
library(goseq)
library(clusterProfiler)
library(tximport)
library(GenomicFeatures)
library(ggfortify)
library(rgl)
library(car)
library(plot3D)

# Before importing the rMTAS outpot files (JC.txt), make sure the first row (header) has been removed, the files are renamed as fixed.txt
# Here we detect each type of alternative splicing events/genes across all clean samples (all samples excluding the outlier DS2_P15) as an example.

#A3SS
A3SS=read.table("A3SS_all_fixed.txt",header=FALSE)
colnames(A3SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
A3SS=transform(A3SS,IJC=colsplit(IJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
A3SS=transform(A3SS,SJC=colsplit(SJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
A3SS=transform(A3SS,IncLevel=colsplit(IncLevel,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))

A3SS$IJC_sum=rowSums(A3SS$IJC)
A3SS$SJC_sum=rowSums(A3SS$SJC)
A3SS$IncLevel_mean=rowMeans(A3SS$IncLevel,na.rm = TRUE)
write.csv(A3SS,file="A3SS_all.csv")
A3SS_count_over_20=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]
write.csv(A3SS_count_over_20,file="A3SS_all_count_over_20.csv")
A3SS_PSI_10=A3SS_count_over_20[A3SS_count_over_20$IncLevel_mean>0.1&0.9>A3SS_count_over_20$IncLevel_mean,]
write.csv(A3SS_PSI_10,file="A3SS_all_PSI_10.csv")

#A5SS
A5SS=read.table("A5SS_all_fixed.txt",header=FALSE)
colnames(A5SS)=c("ID","GeneID","geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES","flankingEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
A5SS=transform(A5SS,IJC=colsplit(IJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
A5SS=transform(A5SS,SJC=colsplit(SJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
A5SS=transform(A5SS,IncLevel=colsplit(IncLevel,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))

A5SS$IJC_sum=rowSums(A5SS$IJC)
A5SS$SJC_sum=rowSums(A5SS$SJC)
A5SS$IncLevel_mean=rowMeans(A5SS$IncLevel,na.rm = TRUE)
write.csv(A5SS,file="A5SS_all.csv")
A5SS_count_over_20=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]
write.csv(A5SS_count_over_20,file="A5SS_all_count_over_20.csv")
A5SS_PSI_10=A5SS_count_over_20[A5SS_count_over_20$IncLevel_mean>0.1&0.9>A5SS_count_over_20$IncLevel_mean,]
write.csv(A5SS_PSI_10,file="A5SS_all_PSI_10.csv")

#SE
SE=read.table("SE_all_fixed.txt",header=FALSE)
colnames(SE)=c("ID","GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
SE=transform(SE,IJC=colsplit(IJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
SE=transform(SE,SJC=colsplit(SJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
SE=transform(SE,IncLevel=colsplit(IncLevel,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))

SE$IJC_sum=rowSums(SE$IJC)
SE$SJC_sum=rowSums(SE$SJC)
SE$IncLevel_mean=rowMeans(SE$IncLevel,na.rm = TRUE)
write.csv(SE,file="SE_all.csv")
SE_count_over_20=SE[SE$IJC_sum>20&SE$SJC_sum>20,]
write.csv(SE_count_over_20,file="SE_all_count_over_20.csv")
SE_PSI_10=SE_count_over_20[SE_count_over_20$IncLevel_mean>0.1&0.9>SE_count_over_20$IncLevel_mean,]
write.csv(SE_PSI_10,file="SE_all_PSI_10.csv")

#MXE
MXE=read.table("MXE_all_fixed.txt",header=FALSE)
colnames(MXE)=c("ID","GeneID","geneSymbol","chr","strand","1stExonStart_0base","1stExonEnd","2stExonStart_0base","2stExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
MXE=transform(MXE,IJC=colsplit(IJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
MXE=transform(MXE,SJC=colsplit(SJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
MXE=transform(MXE,IncLevel=colsplit(IncLevel,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))

MXE$IJC_sum=rowSums(MXE$IJC)
MXE$SJC_sum=rowSums(MXE$SJC)
MXE$IncLevel_mean=rowMeans(MXE$IncLevel,na.rm = TRUE)
write.csv(MXE,file="MXE_all.csv")
MXE_count_over_20=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]
write.csv(MXE_count_over_20,file="MXE_all_count_over_20.csv")
MXE_PSI_10=MXE_count_over_20[MXE_count_over_20$IncLevel_mean>0.1&0.9>MXE_count_over_20$IncLevel_mean,]
write.csv(MXE_PSI_10,file="MXE_all_PSI_10.csv")

#RI
RI=read.table("RI_all_fixed.txt",header=FALSE)
colnames(RI)=c("ID","GeneID","geneSymbol","chr","strand","riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","ID","IJC","SJC","lncFormLen","SkipFormLen","PValue","FDR","IncLevel","IncLevelDifference")
RI=transform(RI,IJC=colsplit(IJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
RI=transform(RI,SJC=colsplit(SJC,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))
RI=transform(RI,IncLevel=colsplit(IncLevel,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60','WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60','DS1_PP50','DS2_PP50','DS3_PP50','DS4_PP50','WS1_PP50','WS2_PP50','WS3_PP50','WS4_PP50','DS1_P15','DS3_P15','DS4_P15','WS1_P15','WS2_P15','WS3_P15','WS4_P15','DS1_P50','DS2_P50','DS3_P50','DS4_P50','WS1_P50','WS2_P50','WS3_P50','WS4_P50')))

RI$IJC_sum=rowSums(RI$IJC)
RI$SJC_sum=rowSums(RI$SJC)
RI$IncLevel_mean=rowMeans(RI$IncLevel,na.rm = TRUE)
write.csv(RI,file="RI_all.csv")
RI_count_over_20=RI[RI$IJC_sum>20&RI$SJC_sum>20,]
write.csv(RI_count_over_20,file="RI_all_count_over_20.csv")
RI_PSI_10=RI_count_over_20[RI_count_over_20$IncLevel_mean>0.1&0.9>RI_count_over_20$IncLevel_mean,]
write.csv(RI_PSI_10,file="RI_all_PSI_10.csv")

# Number of genes associated with the splicing events
genes_A3SS=unique(as.data.frame(A3SS_PSI_10[,2]))
colnames(genes_A3SS)="gene"
genes_A5SS=unique(as.data.frame(A5SS_PSI_10[,2]))
colnames(genes_A5SS)="gene"
genes_MXE=unique(as.data.frame(MXE_PSI_10[,2]))
colnames(genes_MXE)="gene"
genes_SE=unique(as.data.frame(SE_PSI_10[,2]))
colnames(genes_SE)="gene"
genes_RI=unique(as.data.frame(RI_PSI_10[,2]))
colnames(genes_RI)="gene"
genes_all_clean=rbind(genes_A3SS,genes_A5SS)
genes_all_clean=rbind(genes_all_clean,genes_MXE)
genes_all_clean=rbind(genes_all_clean,genes_SE)                 
genes_all_clean=rbind(genes_all_clean,genes_RI) 
genes_all_clean=unique(genes_all_clean)

# Create colomn data to plot all filtered splicing events across all clean samples.
A3SS_PSI10=read.csv(file="A3SS_all_PSI_10.csv")
A5SS_PSI10=read.csv(file="A5SS_all_PSI_10.csv")
MXE_PSI10=read.csv(file="MXE_all_PSI_10.csv")
RI_PSI10=read.csv(file="RI_all_PSI_10.csv")
SE_PSI10=read.csv(file="SE_all_PSI_10.csv")
A3SS_PSI10$GeneID <- paste("A3SS", substring(A3SS_PSI10$GeneID,6), sep="_")
A5SS_PSI10$GeneID <- paste("A5SS", substring(A5SS_PSI10$GeneID,6), sep="_")
MXE_PSI10$GeneID <- paste("MXE", substring(MXE_PSI10$GeneID,6), sep="_")
RI_PSI10$GeneID <- paste("RI", substring(RI_PSI10$GeneID,6), sep="_")
SE_PSI10$GeneID <- paste("SE", substring(SE_PSI10$GeneID,6), sep="_")
all_PSI10=rbind(A3SS_PSI10[80:110],A5SS_PSI10[80:110],MXE_PSI10[82:112],SE_PSI10[80:110],RI_PSI10[80:110])
colnames(all_PSI10)=substring(colnames(all_PSI10),10)
all_PSI10_no_na=all_PSI10[complete.cases(all_PSI10),]
Stage <- factor(rep(c("Wr60","PP50","P15","P50"),c(8,8,7,8)))
Form=c(rep(rep(c("DS","WS"),c(4,4)),2),rep(c("DS","WS"),c(3,4)),rep(c("DS","WS"),c(4,4)))
col_data <- data.frame(row.names = substring(colnames(A3SS_PSI10[80:110]),10),Stage,Form)
col_data

# Plot heatmap
ann_colors = list(Form = c(DS="#E69F00", WS="chartreuse3"),Stage=c(Wr60="lightskyblue",PP50="navy",P15="burlywood",P50="violetred3"))
temp_hm=pheatmap(all_PSI10_no_na,border=NA,scale = "row",annotation_colors = ann_colors,color = colorRampPalette(c("purple","purple", "black", "yellow", "yellow"))(100), annotation=col_data, show_rownames=FALSE,show_colnames = TRUE,annotation_names_col=FALSE,annotation_legend = TRUE)
ggsave(temp_hm,width=6,height=6,filename='all_clean.tiff')
ggsave(temp_hm,width=6,height=6,filename='all_clean.eps')

# Plot 3D PCA
mm_splicing=prcomp(t(all_PSI10_no_na))
summary(mm_splicing) #Record the percentage of the first three major PCs
mm_splicing=as.data.frame(mm_splicing$x)
mm_splicing$Form=c(rep("DS",4),rep("WS",4),rep("DS",4),rep("WS",4),rep("DS",3),rep("WS",4),rep("DS",4),rep("WS",4))
mm_splicing$Stage=c(rep("Wr60",8),rep("PP50",8),rep("P15",7),rep("P50",8))
mm_splicing$color=c(rep("lightskyblue",8),rep("navy",8),rep("burlywood",7),rep("violetred3",8))
mm_splicing$shape=c(rep(4,4),rep(15,4),rep(4,4),rep(15,4),rep(4,3),rep(15,4),rep(4,4),rep(15,4))
scatter3D(x=mm_splicing$PC1,y=mm_splicing$PC2,z=mm_splicing$PC3,par(mar=c(5,0,5,0)),col=mm_splicing$color,
          colvar = NULL,phi=0,bty = "g",type="h",cex=2,cex.lab=1.5, cex.main=2.3,
          main="Alternative splicing",xlab="PC1: 16%",ylab="PC2: 11%",zlab="PC3: 8%",ticktype="detailed",pch=mm_splicing$shape) # Type in the PC values as calculated previously.