######################## DSG analysis ##################################

# Install and activate the required packages.
library("reshape")

# Here we detect each type of differential splicing events/genes between DS and WS forms of Wr60 samples as an example.

#A3SS
A3SS=read.table("A3SS.MATS.JC.txt",header=TRUE)
A3SS=transform(A3SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A3SS=transform(A3SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A3SS=transform(A3SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A3SS=transform(A3SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A3SS=transform(A3SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A3SS=transform(A3SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))

nrow(A3SS)
A3SS=A3SS[A3SS$FDR<0.05,]
write.csv(A3SS,file="A3SS_compare_Wr60_FDR.csv")
nrow(A3SS)
A3SS$IJC_sum=rowSums(A3SS$IJC_SAMPLE_1)+rowSums(A3SS$IJC_SAMPLE_2)
A3SS$SJC_sum=rowSums(A3SS$SJC_SAMPLE_1)+rowSums(A3SS$SJC_SAMPLE_2)
A3SS_count_over_20=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]
write.csv(A3SS_count_over_20,file="A3SS_compare_Wr60_FDR_count20.csv")

#A5SS
A5SS=read.table("A5SS.MATS.JC.txt",header=TRUE)
A5SS=transform(A5SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A5SS=transform(A5SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A5SS=transform(A5SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A5SS=transform(A5SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A5SS=transform(A5SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A5SS=transform(A5SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))

nrow(A5SS)
A5SS=A5SS[A5SS$FDR<0.05,]
write.csv(A5SS,file="A5SS_compare_Wr60_FDR.csv")
nrow(A5SS)
A5SS$IJC_sum=rowSums(A5SS$IJC_SAMPLE_1)+rowSums(A5SS$IJC_SAMPLE_2)
A5SS$SJC_sum=rowSums(A5SS$SJC_SAMPLE_1)+rowSums(A5SS$SJC_SAMPLE_2)
A5SS_count_over_20=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]
write.csv(A5SS_count_over_20,file="A5SS_compare_Wr60_FDR_count20.csv")

#SE
SE=read.table("SE.MATS.JC.txt",header=TRUE)
SE=transform(SE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
SE=transform(SE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
SE=transform(SE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
SE=transform(SE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
SE=transform(SE,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
SE=transform(SE,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))

nrow(SE)
SE=SE[SE$FDR<0.05,]
write.csv(SE,file="SE_compare_Wr60_FDR.csv")
nrow(SE)
SE$IJC_sum=rowSums(SE$IJC_SAMPLE_1)+rowSums(SE$IJC_SAMPLE_2)
SE$SJC_sum=rowSums(SE$SJC_SAMPLE_1)+rowSums(SE$SJC_SAMPLE_2)
SE_count_over_20=SE[SE$IJC_sum>20&SE$SJC_sum>20,]
write.csv(SE_count_over_20,file="SE_compare_Wr60_FDR_count20.csv")

#MXE
MXE=read.table("MXE.MATS.JC.txt",header=TRUE)
MXE=transform(MXE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
MXE=transform(MXE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
MXE=transform(MXE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
MXE=transform(MXE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
MXE=transform(MXE,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
MXE=transform(MXE,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))

nrow(MXE)
MXE=MXE[MXE$FDR<0.05,]
write.csv(MXE,file="MXE_compare_Wr60_FDR.csv")
nrow(MXE)
MXE$IJC_sum=rowSums(MXE$IJC_SAMPLE_1)+rowSums(MXE$IJC_SAMPLE_2)
MXE$SJC_sum=rowSums(MXE$SJC_SAMPLE_1)+rowSums(MXE$SJC_SAMPLE_2)
MXE_count_over_20=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]
write.csv(MXE_count_over_20,file="MXE_compare_Wr60_FDR_count20.csv")
#RI
RI=read.table("RI.MATS.JC.txt",header=TRUE)
RI=transform(RI,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
RI=transform(RI,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
RI=transform(RI,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
RI=transform(RI,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
RI=transform(RI,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
RI=transform(RI,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))

nrow(RI)
RI=RI[RI$FDR<0.05,]
write.csv(RI,file="RI_compare_Wr60_FDR.csv")
nrow(RI)
RI$IJC_sum=rowSums(RI$IJC_SAMPLE_1)+rowSums(RI$IJC_SAMPLE_2)
RI$SJC_sum=rowSums(RI$SJC_SAMPLE_1)+rowSums(RI$SJC_SAMPLE_2)
RI_count_over_20=RI[RI$IJC_sum>20&RI$SJC_sum>20,]
write.csv(RI_count_over_20,file="RI_compare_Wr60_FDR_count20.csv")

# Generate all DSG data
A3SS_count_over_20$Type=rep("A3SS",nrow(A3SS_count_over_20))
A5SS_count_over_20$Type=rep("A5SS",nrow(A5SS_count_over_20))
MXE_count_over_20$Type=rep("MXE",nrow(MXE_count_over_20))
SE_count_over_20$Type=rep("SE",nrow(SE_count_over_20))
RI_count_over_20$Type=rep("RI",nrow(RI_count_over_20))

DSG=rbind(rbind(rbind(rbind(A3SS_count_over_20[,c(2,20,23,26)]
          ,A5SS_count_over_20[,c(2,20,23,26)])
          ,MXE_count_over_20[,c(2,22,25,28)])
          ,SE_count_over_20[,c(2,20,23,26)])
          ,RI_count_over_20[,c(2,20,23,26)])
DSG$absIncleveldifference=abs(DSG$IncLevelDifference)
write.csv(DSG,file="DS_events_Wr60.csv") # Total DS events between seasonal forms during Wr60

# The DS event with the largest absIncleveldifference of a DSG represent the differential splicing level of the gene.
DSG_unique_gene=DSG[order(DSG[,'GeneID'],-DSG[,'absIncleveldifference']),]
DSG_unique_gene=DSG_unique_gene[!duplicated(DSG_unique_gene$GeneID),]
write.csv(DSG_unique_gene,file="DSG_Wr60_unique.csv") # Total unique DS genes between seasonal forms during Wr60


# Generate all gene splicing data, including those not differentially spliced. This is used for the scatter plot for DEG vs DSG comparison.
A3SS=read.table("A3SS.MATS.JC.txt",header=TRUE)
A3SS=transform(A3SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A3SS=transform(A3SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A3SS=transform(A3SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A3SS=transform(A3SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A3SS=transform(A3SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A3SS=transform(A3SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A3SS$IJC_sum=rowSums(A3SS$IJC_SAMPLE_1)+rowSums(A3SS$IJC_SAMPLE_2)
A3SS$SJC_sum=rowSums(A3SS$SJC_SAMPLE_1)+rowSums(A3SS$SJC_SAMPLE_2)
A3SS_all=A3SS[A3SS$IJC_sum>20&A3SS$SJC_sum>20,]
A5SS=read.table("A5SS.MATS.JC.txt",header=TRUE)
A5SS=transform(A5SS,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A5SS=transform(A5SS,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A5SS=transform(A5SS,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A5SS=transform(A5SS,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A5SS=transform(A5SS,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
A5SS=transform(A5SS,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
A5SS$IJC_sum=rowSums(A5SS$IJC_SAMPLE_1)+rowSums(A5SS$IJC_SAMPLE_2)
A5SS$SJC_sum=rowSums(A5SS$SJC_SAMPLE_1)+rowSums(A5SS$SJC_SAMPLE_2)
A5SS_all=A5SS[A5SS$IJC_sum>20&A5SS$SJC_sum>20,]
SE=read.table("SE.MATS.JC.txt",header=TRUE)
SE=transform(SE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
SE=transform(SE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
SE=transform(SE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
SE=transform(SE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
SE=transform(SE,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
SE=transform(SE,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
SE$IJC_sum=rowSums(SE$IJC_SAMPLE_1)+rowSums(SE$IJC_SAMPLE_2)
SE$SJC_sum=rowSums(SE$SJC_SAMPLE_1)+rowSums(SE$SJC_SAMPLE_2)
SE_all=SE[SE$IJC_sum>20&SE$SJC_sum>20,]
MXE=read.table("MXE.MATS.JC.txt",header=TRUE)
MXE=transform(MXE,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
MXE=transform(MXE,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
MXE=transform(MXE,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
MXE=transform(MXE,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
MXE=transform(MXE,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
MXE=transform(MXE,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
MXE$IJC_sum=rowSums(MXE$IJC_SAMPLE_1)+rowSums(MXE$IJC_SAMPLE_2)
MXE$SJC_sum=rowSums(MXE$SJC_SAMPLE_1)+rowSums(MXE$SJC_SAMPLE_2)
MXE_all=MXE[MXE$IJC_sum>20&MXE$SJC_sum>20,]
RI=read.table("RI.MATS.JC.txt",header=TRUE)
RI=transform(RI,IJC_SAMPLE_1=colsplit(IJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
RI=transform(RI,SJC_SAMPLE_1=colsplit(SJC_SAMPLE_1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
RI=transform(RI,IJC_SAMPLE_2=colsplit(IJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
RI=transform(RI,SJC_SAMPLE_2=colsplit(SJC_SAMPLE_2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
RI=transform(RI,IncLevel1=colsplit(IncLevel1,split=",",names=c('DS1_Wr60','DS2_Wr60','DS3_Wr60','DS4_Wr60')))
RI=transform(RI,IncLevel2=colsplit(IncLevel2,split=",",names=c('WS1_Wr60','WS2_Wr60','WS3_Wr60','WS4_Wr60')))
RI$IJC_sum=rowSums(RI$IJC_SAMPLE_1)+rowSums(RI$IJC_SAMPLE_2)
RI$SJC_sum=rowSums(RI$SJC_SAMPLE_1)+rowSums(RI$SJC_SAMPLE_2)
RI_all=RI[RI$IJC_sum>20&RI$SJC_sum>20,]

# Combine
A3SS_all$Type=rep("A3SS",nrow(A3SS_all))
A5SS_all$Type=rep("A5SS",nrow(A5SS_all))
MXE_all$Type=rep("MXE",nrow(MXE_all))
SE_all$Type=rep("SE",nrow(SE_all))
RI_all$Type=rep("RI",nrow(RI_all))

SG_all=rbind(rbind(rbind(rbind(A3SS_all[,c(2,20,23,26)]
                               ,A5SS_all[,c(2,20,23,26)])
                         ,MXE_all[,c(2,22,25,28)])
                   ,SE_all[,c(2,20,23,26)])
             ,RI_all[,c(2,20,23,26)])
SG_all$absIncleveldifference=abs(SG_all$IncLevelDifference)

# Only keep the splicing event with the max asbIncleveldiff for a unique gene list.
SG_unique_gene=SG_all[order(SG_all[,'GeneID'],-SG_all[,'absIncleveldifference']),]
SG_unique_gene=SG_unique_gene[!duplicated(SG_unique_gene$GeneID),]
write.csv(SG_unique_gene,file="Splicing_Wr60_unique.csv") # This is used for the DEG vs DSG scatter plot.




