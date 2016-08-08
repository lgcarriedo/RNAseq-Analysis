# RNAseq-Analysis
# The RISE stress experiments were repeated in two different years: 2014 and 2015. 

# This analysis was done in four steps:
# 1) I have combined all the data from two years together to perform one differential gene expression analysis. I used the edgeR #‘exactTest’ function on the data.

#2) I then performed a gene ontology (GO) enrichment on all the DE genes that were sig. expressed by an FDR cut-off of <0.05. 

#3) Once I had the GO terms, I found all the genes in the dataset that were classified under particular ‘GO’ terms and calculated an #average logFC of each term. 

#4) Then I wanted to know if the same genes were represented across all ‘GO’ terms within each species.

#To replicate the analysis, use: “RISE_EdgeR_CombinedYears.R”.
