#######################comparison between DEGs and DSGs###################

# Install and activate the required packages.
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
library(dbplyr)
library(goseq)
library(clusterProfiler)
library(tximport)
library(GenomicFeatures)

# Read the DEG and DSG output files for comparison, the files used here are generated in the DEG and DSG analysis. Here we used the Wr60 results as an example.
DEG_Wr60=read_csv("DEG_Wr60_WSvsDS.csv",col_names = TRUE)
DEG_Wr60=as.data.frame(DEG_Wr60)
DEG_all_Wr60=read_csv("GE_Wr60_WSvsDS.csv",col_names = TRUE)
DEG_all_Wr60=as.data.frame(DEG_all_Wr60)
DSG_Wr60=read_csv("DSG_Wr60_unique.csv",col_names = TRUE)
DSG_Wr60=as.data.frame(DSG_Wr60)
DSG_Wr60$GeneID=substring(DSG_Wr60$GeneID,6)
DSG_all_Wr60=read_csv("SG_Wr60_unique.csv",col_names = TRUE)
DSG_all_Wr60=as.data.frame(DSG_all_Wr60)
DSG_all_Wr60$GeneID=substring(DSG_all_Wr60$GeneID,6)

# Create Venn plot showing how DEGs and DSGs overlap.
myCol <- c("orange","navy")
venn_DE=list(DSG_Wr60 = DSG_Wr60$GeneID,DEG_Wr60 = DEG_Wr60$Row.names)
venn.diagram(cat.cex = 2,cex = 2,cat.pos = c(0,0),compression = "lzw",x=venn_DE ,"Wr60.png",col = myCol,lwd=10,fontfamily = "Helvetica",cat.fontfamily ="Helvetica" )

# Create scatter plot to show whether DEGs are also DSGs, and vise versa.
merge_Wr60=data.frame(gene=c(union(DEG_Wr60$Row.names,DSG_Wr60$GeneID)))
merge_Wr60$FC=DEG_Wr60$log2FoldChange[match(merge_Wr60$gene,DEG_Wr60$Row.names)]
merge_Wr60$difflevel=DSG_Wr60$IncLevelDifference[match(merge_Wr60$gene,DSG_Wr60$GeneID)]
# Write type of gene.
for (n in 1:nrow(merge_Wr60)){
  if(is.na(merge_Wr60[n,2])==FALSE&is.na(merge_Wr60[n,3])==FALSE) {
      merge_Wr60[n,4]="DEG+DSG"
  }
  if(is.na(merge_Wr60[n,2])==TRUE&is.na(merge_Wr60[n,3])==FALSE) {
    merge_Wr60[n,4]="DSG"}
    if(is.na(merge_Wr60[n,2])==FALSE&is.na(merge_Wr60[n,3])==TRUE) {
      merge_Wr60[n,4]="DEG"}
}
# Assign non-significant values.
for (n in 1:nrow(merge_Wr60)){
  if(is.na(merge_Wr60[n,2])==TRUE){
    merge_Wr60[n,2]=DEG_all_Wr60$log2FoldChange[match(merge_Wr60[n,1],DEG_all_Wr60$X1)]
  }
  if(is.na(merge_Wr60[n,3])==TRUE){
    merge_Wr60[n,3]=DSG_all_Wr60$IncLevelDifference[match(merge_Wr60[n,1],DSG_all_Wr60$GeneID)]
  }
}
# Asign na values as 0.
for (n in 1:nrow(merge_Wr60)){
  if(is.na(merge_Wr60[n,2])==TRUE){
    merge_Wr60[n,2]=0
  }
  if(is.na(merge_Wr60[n,3])==TRUE){
    merge_Wr60[n,3]=0
  }
}
merge_Wr60$V4=factor(merge_Wr60$V4,levels = c("DSG","DEG","DEG+DSG"))
names(merge_Wr60)[4]=paste("Type")
summary(merge_Wr60)

#scatter plot
plot=ggplot(merge_Wr60%>%arrange(Type),aes(x=FC,y=difflevel,color=Type))+geom_point(size=0.1)+theme_bw()+scale_color_manual(values=c('orange','navy', 'chartreuse3'))
plot=plot+theme(legend.position = c(0.2,0.85),legend.title = element_blank())+ guides(colour = guide_legend(override.aes = list(size=1)))+theme(panel.grid = element_blank())
plot=plot+xlab("Gene expression (log2FC)")+ylab(expression("Alternative splicing ("~Delta~"PSI)"))+theme(legend.key.size = unit(0.1,"line"))
plot=plot+xlim(-10,10)+ylim(-1,1)
plot=plot+coord_fixed(ratio=10)
plot
ggsave(plot,width=3.5,height=3.5,filename='Wr60.pdf')
