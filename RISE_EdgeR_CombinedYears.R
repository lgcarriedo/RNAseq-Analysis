# Analyzing *BOTH* years of gene expression together
# Start with edgeR exactTest
# Then Do a 'GO' analysis
# Then plot avg. logFC of the genes per term that comes out as enriched. 
# August 1, 2016
# RISE DISSERTATION CHAPTER
#setwd("/Volumes/LC_KINGSTON/RISE_DISSERTATION_CHAPTER_ANALYSIS/RNAseq_Analysis_CombinedYears")

#load the libraries#
library(limma)
require(splines)
require(edgeR)
library(reshape2)
library(ggplot2)

#All the Phosphorous RNAseq count data for both years:
#expr <-read.csv("/Users/leonelacarriedo/Desktop/Reciprocal_Mapped_Reads/RISE_CountData/M82_2014&2015_Phosphorus_CountData.csv",row.names=1)
#expr <-read.csv("/Users/leonelacarriedo/Desktop/Reciprocal_Mapped_Reads/RISE_CountData/PIM_2014&2015_Phosphorus_CountData.csv",row.names=1)
#expr <-read.csv("/Users/leonelacarriedo/Desktop/Reciprocal_Mapped_Reads/RISE_CountData/PEN_2014&2015_Phosphorus_CountData.csv",row.names=1)
#expr <-read.csv("/Users/leonelacarriedo/Desktop/Reciprocal_Mapped_Reads/RISE_CountData/NEO_2014&2015_Phosphorus_CountData.csv",row.names=1)

#Nitrogen RNAseq count data for both years
#Nitrogen
#expr <-read.csv("/Users/leonelacarriedo/Desktop/Reciprocal_Mapped_Reads/RISE_CountData/M82_2014&2015_Nitrogen_CountData_2015_switched.csv",row.names=1)
#expr <- read.csv("/Users/leonelacarriedo/Desktop/Reciprocal_Mapped_Reads/RISE_CountData/PIM_2014&2015_Nitrogen_CountData_switched.csv",row.names=1)
#expr <-read.csv("/Users/leonelacarriedo/Desktop/Reciprocal_Mapped_Reads/RISE_CountData/PEN_2014&2015_Nitrogen_CountData.csv",row.names=1)
#expr<-read.csv("/Users/leonelacarriedo/Desktop/Reciprocal_Mapped_Reads/RISE_CountData/Neo_2014&2015_Nitrogen_CountData.csv",row.names=1)

head(expr)[,1:5]
names(expr) #names look good.
dim(expr)

#Remove the S. penn specific genes
expr<-expr[!expr$INFO=="NA",]
expr<-expr[complete.cases(expr),] #remove all introduced NAs
expr<-expr[,!colnames(expr)=="INFO"] #remove the "Sopen_unique" column
rownames(expr)<-expr$ITAG
expr<-expr[,!colnames(expr)=="ITAG"] #remove the "X" column

##### set up the model #####
rep  <- factor(substring(colnames(expr), 13,13)) #pull out the "rep" info
trt  <- relevel(factor(substring(colnames(expr),7,10)),ref="cont")

data1 <- expr

group <- as.factor(trt)

#### make the model####
expr.dge <- DGEList(counts = data1,group=group)

expr.dge <- expr.dge[rowSums(cpm(expr.dge$counts) > 5) >= 5, ]
dim(expr.dge)

expr.dge$samples$lib.size <- colSums(expr.dge$counts)  

#### get Normalization factors for library size ####
expr.norm <- calcNormFactors(expr.dge, refColumn = 1, method = "TMM")

#estimate common and then tagwise Dispersion of each gene
expr.c <- estimateCommonDisp(expr.norm)

expr.tgw <- estimateTagwiseDisp(expr.c)

#plotBCV(expr.tgw)

plotMDS(expr.tgw)

#perform the exactTest to perform DE between treatments
et <- exactTest(expr.tgw)

topTags(et)

#how many genes pass the applied FDR cut-off?
sum(topTags(et, n = Inf)$table$FDR < 0.05)

#Extract all the DE genes
gen_DE <-topTags(et, n=Inf)

#Add annotations to the DE tables
annotation <- read.csv("~/Google Drive/OLD_RNAseq_Analyses/RISE- Joint gene expression Analysis 2015/2016_RNAseq_RISE_Results/ITAG2.3_HRD.csv")
gen_DEannotated <- merge(gen_DE, annotation, by.x=0, by.y=1, sort = F)

# save the output as a .csv
#write.csv(gen_DEannotated, "NEO_ExactTest_Nitrogen_DE_YrsCombined_DE_08012016.csv")

##########################
#DO GO ANALYSIS#
##########################
#source("https://bioconductor.org/biocLite.R")
#biocLite("goseq")

library(goseq)
library(GO.db)

universe<-read.csv("NEO_ExactTest_PhosphorusDE_YrsCombined_DE_08012016.csv",row.names = 2)
#rownames(universe)<-gen_DEannotated$Row.names
universe<-universe[,-1] #remove "X" column
head(universe)[,1:4]

sol_trans_len <- read.csv("~/Downloads/gofun/R scripts for DE GO/ITAG2.3_cds.mod.length.tsv", sep="\t", h=T)
colnames(sol_trans_len)<-c("Gene","length")
sol_trans_lenf <- sol_trans_len[,1]
sol_trans_len2 <- sol_trans_len[,2]
genes.in.annot<-sol_trans_lenf

#Go file
go <- read.table("~/Downloads/gofun/R scripts for DE GO/GO.table.tsv", h=T)
head(go)[,1:2]
colnames(go) <- c("Gene", "GO")
go.list <- strsplit(as.character(go[,2]),split=",",fixed=T)
names(go.list) <- as.character(go[,1])
length(go.list)

go.list2<-go.list[names(go.list)%in%row.names(universe)] #this reduces the go.list to reflect the tot. num genes
length(go.list2) 
head(go.list2)
go.list<-go.list2

#table of DE genes, not FDR filtered 
genes <- universe
names(genes)
head(genes)[,1:4]

#genes.DE <- as.numeric(genes$FDR < 0.01) # this assigns a 1 or a 0 to the genes that meet the criteria
genes.DE <- as.numeric(genes$FDR <0.05)
length(genes.DE)
names(genes.DE) <- rownames(genes)
bias<-sol_trans_len2
names(bias) <- genes.in.annot
new_bias <- bias[names(bias)%in%names(genes.DE)]
length(new_bias)

#something weird happened, where the new_bias and genes.DE did not match...so let's fix it.
new_bias<-bias[which(names(bias) %in% names(genes.DE))]
length(new_bias) 

genes.DE<-genes.DE[which (names(new_bias) %in% names(genes.DE))]
length(genes.DE)

#Calculates a Probability Weighting Function for a set of genes based on a given set of biased data 
#(usually gene length) and each genes status as differentially expressed or not.
pwf<- nullp(genes.DE,bias.data=new_bias)
head(pwf)

#Does selection-unbiased testing for category enrichment amongst differentially expressed (DE) genes 
#for RNA-seq data. By default, tests gene ontology (GO) categories, but any categories may be tested
go.analysis <- goseq(pwf,gene2cat=go.list)
print(length(go.analysis$category[p.adjust(go.analysis$over_represented_pvalue,method="fdr")<0.05]))

go.analysis$p_val_FDR <-p.adjust(go.analysis$over_represented_pvalue,method="fdr")

enriched.GO_category <-go.analysis$category[go.analysis$p_val_FDR < 0.05]
enriched.GO_p_value_FDR <-go.analysis$p_val_FDR[go.analysis$p_val_FDR < 0.05]
enriched.GO_numDEInCat <-go.analysis$numDEInCat[go.analysis$p_val_FDR < 0.05]
enriched.GO_numInCat <-go.analysis$numInCat[go.analysis$p_val_FDR < 0.05]
enriched.GO<-as.data.frame(cbind(enriched.GO_category, enriched.GO_p_value_FDR, enriched.GO_numDEInCat, enriched.GO_numInCat))
enriched.GO$term <- Term(as.character(enriched.GO$enriched.GO_category))
enriched.GO$ont <- Ontology(as.character(enriched.GO$enriched.GO_category))
#enriched.GO.BP<-enriched.GO[enriched.GO$ont=="BP",]

head(enriched.GO)
dim(enriched.GO)

unlist.go.genes<-do.call(rbind, lapply(go.list2, data.frame, stringsAsFactors=FALSE))
head(unlist.go.genes)
colnames(unlist.go.genes)[1]<- "GO"
unlist.go.genes$ITAG<-substr(row.names(unlist.go.genes),1,18)

#Here I am going to look at the differentially expressed genes within the GO categories
expsn<-merge(enriched.GO, unlist.go.genes, by.x="enriched.GO_category", by.y="GO")
#look at the DE gene list & apply an FDR cut-off
genes2<-genes[genes$FDR<0.05,]
genes2$ITAG<-rownames(genes2)
expsn<- merge(expsn, genes2, by.x="ITAG", by.y="ITAG")
expsn<- expsn[!duplicated(expsn$ITAG),]
head(expsn)
dim(expsn)

#write.csv(expsn, "NEO_Pdef_YrsCombined_GO_EnrichedGenes_08012016.csv")

############################################
#Calculate Avg. logFC of each sig GO term #
############################################

GOgenes<-read.csv("M82_Ndef_YrsCombined_GO_EnrichedGenes_08012016.csv")
GOgenes<- GOgenes[!duplicated(GOgenes$ITAG),]

GOgenes<-GOgenes[,-1]

#take a look at specific columns
head(GOgenes)[,c(1,6,8)]

means <- as.data.frame(tapply(GOgenes$logFC,list(GOgenes$term),mean))

#create an sem function
sem <- function(x) { sd(x) /sqrt(length(x))}
sems <- as.data.frame(tapply(GOgenes$logFC,list(GOgenes$term),sem))

avgs<-cbind(means, sems)
names(avgs)<-c("Avg_logFC","Avg_logFC_sems")

#save the 'avgs' df
#m82P<-avgs
m82P$species<-"M82"
m82P$Stress_type<-"Phosphorus"

#pimP<-avgs
pimP$species<-"PIM"
pimP$Stress_type<-"Phosphorus"

#penP<-avgs
penP$species<-"PEN"
penP$Stress_type<-"Phosphorus"

#neoP<-avgs
neoP$species<-"NEO"
neoP$Stress_type<-"Phosphorus"

##########
#m82N<-avgs
m82N$species<-"M82"
m82N$Stress_type<-"Nitrogen"

#pimN<-avgs
pimN$species<-"PIM"
pimN$Stress_type<-"Nitrogen"

#penN<-avgs
penN$species<-"PEN"
penN$Stress_type<-"Nitrogen"

#neoN<-avgs
neoN$species<-"NEO"
neoN$Stress_type<-"Nitrogen"

#avgP_GOExp<-rbind(m82P, pimP, penP, neoP)
#avgN_GOExp<-rbind(m82N, pimN, penN, neoN)

all_avgs<-rbind(avgP_GOExp, avgN_GOExp)

#write.csv(all_avgs, "AverageLogFC_GOterms_NutrientDep_08012016.csv")
all_avgs<-read.csv("AverageLogFC_GOterms_NutrientDep_08012016.csv")

#all_avgs$term<-rownames(all_avgs)
#avgP_GOExp$term<-rownames(avgP_GOExp)
#avgN_GOExp$term<-rownames(avgN_GOExp)

avgP_GOExp<-all_avgs[all_avgs$Stress_type=="Phosphorus",]

avgP_GOExp$species<- ordered(avgP_GOExp$species, levels = c("M82", "PIM","NEO", "PEN"))
#avgN_GOExp$species<- ordered(avgN_GOExp$species, levels = c("M82", "PIM","NEO", "PEN"))

pl<- ggplot(droplevels(avgP_GOExp), aes(x=term, y=Avg_logFC, fill = species, ymin = Avg_logFC - Avg_logFC_sems,
                                   ymax = Avg_logFC + Avg_logFC_sems))  
pl <- pl + geom_bar(stat="identity") + facet_wrap( ~ species, scales="free", drop = TRUE, ncol=4)  
pl <- pl + geom_errorbar(position=position_dodge(width=.9),width=.25)
pl <- pl + theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ylab("Average logFC") + xlab("")
pl <- pl + theme(strip.text.x = element_text(size = 18,hjust = 0.5, vjust = 0.5)) + theme(strip.text.y = element_text(size = 18,hjust = 0.5, vjust = 0.5))
pl <- pl + theme(axis.text.x  = element_text(vjust=0.5, size=10))
pl <- pl + theme(axis.text.y  = element_text(vjust=0.5, size=18))
pl <- pl + theme(axis.title.y = element_text(size = 18))
pl


M82<-read.csv("M82_ExactTest_Phosphorus_DE_YrsCombined_DE_08012016.csv")
PIM<-read.csv("PIM_ExactTest_PhosphorusDE_YrsCombined_DE_08012016.csv")
PEN<-read.csv("PEN_ExactTest_PhosphorusDE_YrsCombined_DE_08012016.csv")
NEO<-read.csv("NEO_ExactTest_PhosphorusDE_YrsCombined_DE_08012016.csv")

M82<-M82[M82$FDR<0.05,]
PIM<-PIM[PIM$FDR<0.05,]
PEN<-PEN[PEN$FDR<0.05,]
NEO<-NEO[NEO$FDR<0.05,]

dim(M82)
dim(PIM)
dim(PEN)
dim(NEO)

M82<-as.data.frame(M82[2])
PIM<-as.data.frame(PIM[2])
PEN<-as.data.frame(PEN[2])
NEO<-as.data.frame(NEO[2])

write.csv(NEO, "NEO_ITAGS_postFDR0.05_08022016.csv")
