library(GEOquery)
library(simpleaffy)
library(limma)
library(annotate)
library(mouse430a2.db)
library(RColorBrewer)

#retrieve GEO data
getGEOSuppFiles("GSE72517")
untar("GSE72517_RAW.tar", exdir="GSE72517_RAW")
cels <- list.files("GSE72517_RAW/", pattern = "[gz]")
sapply(paste("GSE72517_RAW", cels, sep="/"), gunzip)

# Create metadata in terminal/cygwin, go to folder with .CEL files
# ls *.CEL > phenodata.txt
# Open this file in a text editor. Currently this is a single column listing the file names. The final file needs 3 tab-delimited columns named 'Name' 'FileName' and 'Target'. The FileName and Name columns will be identical in this case, although the Name column can be used to provide a more human readable sample name if desired. The Target column is information about the samples that define our groups for differential expression analysis.

# load files
GSE72517 <- read.affy(covdesc="phenodata.txt", path="GSE72517/GSE72517_RAW")
GSE72517_exp <- exprs(rma(GSE72517))

#Retain only top 25% expressing genes
cutoff = quantile(GSE72517_exp, 0.75) #define cutoff for top 25% of expressing genes
n=5 #number of replicates in each group (number of expression values for each gene you want to subset as being over the cutoff); its recommended that if one group has more samples than another, use the higher number.
isexpr = rowSums(GSE72517_exp > cutoff) >= n #the rows where expression values for n samples exceed the cutoff
GSE72517_exp2 <- GSE72517_exp[isexpr, ] #subset rows

# find differentially expressed probesets
samples72517 <- GSE72517$Target
samples72517 # to verify Targets
samples_factors72517 <- as.factor(samples72517) # turn Targets into factors
samples_factors72517 # to verify factors

# create experimental design
design72517 <- model.matrix(~0 + samples_factors72517)
colnames(design72517) <- c("control_0hrs", "control_8hrs", "control_72hrs", "control_7days", "alcohol_0hrs", "alcohol_8hrs", "alcohol_72hrs", "alcohol_7days")
#colnames(design) <- c("aj_ova_24", "aj_ova_6", "aj_pbs_24", "aj_pbs_6", "c3h_ova_24", "c3h_ova_6", "c3h_pbs_24", "c3h_pbs_6")
design72517 # to verify design

# load the data into limma
# fit the linear model to the filtered expression set
fit72517 <- lmFit(GSE72517_exp, design72517)
# set up a contrast matrix to compare tissues v cell line
contrast.matrix72517 <- makeContrasts(alcohol7days_control7days = alcohol_7days - control_7days, levels=design72517)
#contrast.matrix <- makeContrasts(aj_24 = aj_pbs_24 - aj_ova_24, aj_6 = aj_pbs_6 - aj_ova_6, c3h_24 = c3h_pbs_24 - c3h_ova_24, c3h_6 = c3h_pbs_6 - c3h_ova_6, levels=design)
# check the contrast matrix
contrast.matrix72517


# Now the contrast matrix is combined with the per-probeset linear model fit.
combine_fits72517 <- contrasts.fit(fit72517, contrast.matrix72517)
combine_ebFit72517 <- eBayes(combine_fits72517)
# return the top 10 results for any given contrast
# coef=1 is pbs_ova
# top10 <- topTable(aj_ebFit, number=10, coef=1) # look at top10 lowest p-values
# probeset.list <- topTable(aj_ebFit, coef=1, number=10000, lfc=4)

# output table
genetable72517 <- topTable(combine_ebFit72517, number=9999999999999999, coef=1)
# list of probes
probe.list72517 <- row.names(genetable72517)

gene.symbols72517 <- getSYMBOL(probe.list72517, "mouse430a2")
results72517 <- cbind(gene.symbols72517, genetable72517)
head(results72517)
write.table(results72517, file="GSE72517_7days.txt", quote=F, sep="\t", col.names=T)