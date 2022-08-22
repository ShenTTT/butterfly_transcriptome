# Purpose
These scripts are used to clean raw small RNA (sRNA)-seq data for the purpose of miRNA annotation, and to annotate the composition of sRNAs.

# Perform sRNA-seq data curation to annotate miRNAs

1. Trim adaptors from the raw sRNA-seq data as described.

2. Remove rRNAs, tRNAs, snRNAs, and snoRNAs using sortmeRNA (sortmerna.sh) to generate clean sRNA-seq data for miRNA annotation. The complete list of known sRNA databases were listed in the original publication. Please refer to the sortmeRNA manual for the proper format of the reference databases. Different catagories of sRNAs are removed one after another. For instance, clean reads after the removal of rRNAs are further processed to remove tRNAs, and then snRNAs, and so on. 

# Annotate sRNA composition 

1. Map the adaptor trimmed sRNA-seq reads onto the genome using STAR (STAR.sh)

2. Annotate the composition of sRNAs using unitas (unitas.sh)




