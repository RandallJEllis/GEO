setwd("")
library(RankProd)
library(simpleaffy)
library(mouse4302cdf)
library(AnnotationDbi)
library(mouse4302.db)
library(metaArray)
library(MergeMaid)
library(moe430a.db)
library(moe430acdf)
library(annotate)
library(mouse430a2.db)

#RankProduct analysis in R. 1, 2, and 3 are placeholders for GSE accession numbers.

# load files
GSE <- read.affy(covdesc="phenodata.txt", path="GSE1/GSE1_RAW")
GSE1_exp <- exprs(rma(GSE1))

GSE2 <- read.affy(covdesc="phenodata.txt", path="GSE2/GSE2_RAW")
GSE2_exp <- exprs(rma(GSE2))

GSE3 <- read.affy(covdesc="phenodata.txt", path="GSE3/GSE3_RAW")
GSE3_exp <- exprs(rma(GSE3))

## We merge all of them

GSEmerged<-mergeExprs(GSE1_exp,GSE2_exp,GSE3_exp)

# intersection gives us a expression set with only common genes to all the studies in the merged object

GSEcommon<-intersection(GSEmerged)

write.exprs(GSEcommon,file="metadata.txt")

metasamples<-read.delim("metasamples.txt")
cl<-c(metasamples$Class)
origin<-c(metasamples$Origin)

data<-read.table("metadata.txt", check.names=FALSE)
data <- data.frame(data)

RP.advance.out<-RPadvance(data,cl,origin,gene.names=row.names(data), num.perm=100,logged=T)

plotRP(RP.advance.out, cutoff=0.05)

gene.table <- topGene(RP.advance.out, cutoff=0.05, method="pfp",logged=TRUE,logbase=2)

probe_upreg=row.names(gene.table$Table1)
#gene.symbols <- getSYMBOL(probe_upreg, "mouse430a2")
gene.symbols_upreg <- getSYMBOL(probe_upreg, "mouse430a2")
results_upreg <- cbind(gene.symbols_upreg, gene.table$Table1)
head(results_upreg)
symbols_upreg = unlist(mget(probe_upreg, mouse430a2SYMBOL, ifnotfound=NA))
gene.table$Table1 = cbind(symbols_upreg, gene.table$Table1)
write.table(gene.table$Table1, file="all_upreg_genes2.txt", quote=F, sep="\t", col.names=T)

probe_downreg=row.names(gene.table$Table2)
#gene.symbols <- getSYMBOL(probe_downreg, "mouse430a2")
gene.symbols_downreg <- getSYMBOL(probe_downreg, "mouse430a2")
results_downreg <- cbind(gene.symbols_downreg, gene.table$Table2)
head(results_downreg)
symbols_downreg = unlist(mget(probe_downreg, mouse430a2SYMBOL, ifnotfound=NA))
gene.table$Table2 = cbind(symbols_downreg, gene.table$Table2)
write.table(gene.table$Table2, file="all_downreg_genes2.txt", quote=F, sep="\t", col.names=T)