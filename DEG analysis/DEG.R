######################## DEG analysis ##################################

### Install and activate the required packages.
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

### Generate gene level counts from transcript abundance.
tx2gene=read.table("tx2gene_v1.2.txt",header = TRUE) # tx2gene.txt is a two-column table of all transcript ID (column 1) vs gene ID (column 2)
txi=tximport(c(DS1_Wr60="DS1_Wr60_quant.sf",DS2_Wr60="DS2_Wr60_quant.sf",DS3_Wr60="DS3_Wr60_quant.sf",DS4_Wr60="DS4_Wr60_quant.sf",DS1_PP50="DS1_PP50_quant.sf",DS2_PP50="DS2_PP50_quant.sf",DS3_PP50="DS3_PP50_quant.sf",DS4_PP50="DS4_PP50_quant.sf",DS1_P15="DS1_P15_quant.sf",DS2_P15="DS2_P15_quant.sf",DS3_P15="DS3_P15_quant.sf",DS4_P15="DS4_P15_quant.sf",DS1_P50="DS1_P50_quant.sf",DS2_P50="DS2_P50_quant.sf",DS3_P50="DS3_P50_quant.sf",DS4_P50="DS4_P50_quant.sf",WS1_Wr60="WS1_Wr60_quant.sf",WS2_Wr60="WS2_Wr60_quant.sf",WS3_Wr60="WS3_Wr60_quant.sf",WS4_Wr60="WS4_Wr60_quant.sf",WS1_PP50="WS1_PP50_quant.sf",WS2_PP50="WS2_PP50_quant.sf",WS3_PP50="WS3_PP50_quant.sf",WS4_PP50="WS4_PP50_quant.sf",WS1_P15="WS1_P15_quant.sf",WS2_P15="WS2_P15_quant.sf",WS3_P15="WS3_P15_quant.sf",WS4_P15="WS4_P15_quant.sf",WS1_P50="WS1_P50_quant.sf",WS2_P50="WS2_P50_quant.sf",WS3_P50="WS3_P50_quant.sf",WS4_P50="WS4_P50_quant.sf"),type="salmon",tx2gene=tx2gene)

### Make DEseq2 object
# Make column data
Stage <- factor(rep(rep(c("Wr60","PP50","P15","P50"),2),c(4,4,4,4,4,4,4,4)))
Form=c(rep(c("DS","WS"),c(16,16)))
col_data <- data.frame(row.names = colnames(txi$counts),Stage,Form)

# Make deseq object for:
# all samples
dds_all <- DESeqDataSetFromTximport(txi = txi,colData = col_data,design = ~ Form)
dds_all$Form <- relevel(dds_all$Form, ref = "DS") # Set DS form as the reference.
# clean samples (All samples excluding the outlier DS2_P15, see supplementary materials for more detail)
dds_clean=dds_all[,c(1:9,11:32)]
dds_clean$Form <- relevel(dds_clean$Form, ref = "DS")
# Wr60 samples (as an example to find DEGs between seasonal forms: Wr60_WS vs Wr60_DS)
dds_Wr60=dds_all[,c(1:4,17:20)]
dds_Wr60$Form <- relevel(dds_Wr60$Form, ref = "DS")

### DEG analysis 
# Filter 0 count data
dds=dds_Wr60 # Use Wr60 samples as an example to find DEGs between seasonal forms: Wr60_WS vs Wr60_DS, the same method applies to all other stages.
dds_filter <- dds[ rowSums(counts(dds))>1, ] 
# Perform DE analysis
dds_out <- DESeq(dds_filter)
resultsNames(dds_out)
res <- results(dds_out)
# Rank the list from the highest padj to the lowest.
res_desq <- res[order(res$padj),]
# Select DE genes 
diff_gene_deseq2 <- subset(res_desq, padj<0.05) # Note that both padj<0.05 and abs(log2FoldChange)>1 were used to select shortlisted DEGs for functional enrichment analysis.
# Check the results
summary(diff_gene_deseq2)
# Export the DEG data
res_DE_data <- merge(as.data.frame(diff_gene_deseq2),as.data.frame(counts(dds_out,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(res_DE_data,file = "DEG_Wr60_WSvsDS.csv",row.names = F)

# Generate all non-zero normalized gene expression data (used for DEG vs DSG comparisons)
nor_counts_data <- as.data.frame(counts(dds_out,normalize=TRUE))
write.csv(nor_counts_data,file = "all_confi_nor_gene_countsm.csv",row.names = T)
resdata <- as.data.frame(res[complete.cases(res),])
write.csv(resdata,file = "GE_Wr60_WSvsDS.csv",row.names = T)

### Generate heatmap and PCA plot for all clean samples.
# Heatmap
dds=dds_clean
dds_filter <- dds[ rowSums(counts(dds))>1, ] 
dds_out <- DESeq(dds_filter)
rld <- rlog(dds_out,blind=T) # Perform rlog transformation
mat <- assay(rld)
# Annotate color
ann_colors = list(Form = c(DS="#E69F00", WS="chartreuse3"),Stage=c(Wr60="lightskyblue",PP50="navy",P15="burlywood",P50="violetred3"))
# Plot heatmap
temp_hm=pheatmap(mat,border=NA,scale = "row",annotation_colors = ann_colors,color = colorRampPalette(c("purple","purple", "black","yellow","yellow"))(100), annotation=col_data, show_rownames=FALSE,show_colnames = TRUE,annotation_names_col=FALSE,annotation_legend = TRUE)
# Save images
ggsave(temp_hm,width=6,height=6,filename='all_heat.tiff')
ggsave(temp_hm,width=6,height=6,filename='all_heat.eps')

# 3D PCA
# Calculate PCs
mm_gene=prcomp(t(assay(rld)))
summary(mm_gene) #Record the percentage of the first three major PCs
mm_gene=as.data.frame(mm_gene$x)
# Label the datapoints
mm_gene$Form=c(rep("DS",15),rep("WS",16))
mm_gene$Stage=c(rep("Wr60",4),rep("PP50",4),rep("P15",3),rep("P50",4),rep("Wr60",4),rep("PP50",4),rep("P15",4),rep("P50",4))
mm_gene$color=c(rep("lightskyblue",4),rep("navy",4),rep("burlywood",3),rep("violetred3",4),rep("lightskyblue",4),rep("navy",4),rep("burlywood",4),rep("violetred3",4))
mm_gene$shape=c(rep(4,15),rep(15,16))
# Plot PCA
scatter3D(x=mm_gene$PC1,y=mm_gene$PC2,z=mm_gene$PC3,par(mar=c(5,0,5,0)),col=mm_gene$color,
          colvar = NULL,phi=0,bty = "g",type="h",cex=2, cex.lab=1.5, cex.main=2.3,
          main="Gene expression",xlab="PC1: 41%",ylab="PC2: 24%",zlab="PC3: 10%",ticktype="detailed",pch=mm_gene$shape) # Type in the PC values as calculated previously.
