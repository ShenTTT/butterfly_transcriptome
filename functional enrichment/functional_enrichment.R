#################### functional enrichment#####################################################
# Install and activate required packages.
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

# Import the GO/KEGG annotations.
go2name_data=read.csv("GO2name.csv",header = F,sep=",")
kegg2name_data=read.csv("KEGG2name.csv",header = F,sep=",")
GO2gene_data=read.table("go2gene.txt",header = FALSE)
kegg2gene_data=read.table("kegg2gene.txt",header = FALSE)

# Import gene lists of interest. Here we use shortlisted DSGs between seasonal forms from each stage as an example.
mylist_all <- read.csv("mylist.csv",header = T,sep=",")
names(mylist_all)=c("Wr60","PP50","P15","P50")

# GO enrichment analysis
background=GO2gene_data
clusterenrich <- compareCluster(mylist_all, fun='enricher',TERM2GENE=background,TERM2NAME=go2name_data,pvalueCutoff = 0.1, pAdjustMethod = "BH") # pvalueCutoff can be lifted (changed to 1) if very few terms can be enriched.
# Export all GO enrichment results.
write.csv(as.data.frame(clusterenrich),file="GO_DSG.csv")
# Draw dot plot
dot=dotplot(clusterenrich,showCategory=8,includeAll=TRUE,color="pvalue")
dot=dot+scale_color_continuous(low='orange', high='navy')
dot=dot+ggpubr::rotate_x_text()
dot=dot+scale_y_discrete(label = function(x) stringr::str_trunc(x, 50))
dot=dot+scale_x_discrete(label=c("Wr60","PP50","P15","P50"))
dot
# Save image
ggsave(dot,width=6.5,height=4.5,filename='GO_DSG.tiff')
ggsave(dot,width=6.5,height=4.5,filename='GO_DSG.eps')

# KEGG enrichment analysis
background_kegg=kegg2gene_data
clusterenrich_kegg <- compareCluster(mylist_all, fun='enricher',TERM2GENE=background_kegg,TERM2NAME=kegg2name_data,pvalueCutoff = 1, pAdjustMethod = "BH")
# Export all KEGG enrichment results
write.csv(as.data.frame(clusterenrich_kegg),file="KEGG_DSG.csv")
# Draw dot plot
dot=dotplot(clusterenrich_kegg, showCategory=100,includeAll=TRUE,color="pvalue")
dot=dot+scale_color_continuous(low='orange', high='navy')
dot=dot+ggpubr::rotate_x_text()
dot=dot+scale_yiscrete(label = function(x) stringr::str_trunc(x, 50))
dot=dot+scale_x_discrete(label=c("Wr60","PP50","P15")) #In this case one of the comparisons, P50, has no enriched KEGG terms
dot
# Save image
ggsave(dot,width=6.5,height=4.5,filename='KEGG_DSG.tiff')
ggsave(dot,width=6.5,height=4.5,filename='KEGG_DSG.eps')
