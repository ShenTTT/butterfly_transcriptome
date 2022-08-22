#before read the data make sure the first row of the rmats ouput has been removed
library("reshape")

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

###number of genes
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
