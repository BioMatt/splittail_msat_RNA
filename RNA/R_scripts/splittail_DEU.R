# A script observing differential exon usage in splittail using voom_diff
# Taken from https://github.com/Oshlack/Lace/wiki/Example%3A-Differential-Transcript-Usage-on-a-non-model-organism
setwd("D:/Dropbox/splittail_redone/exon_usage")

library(edgeR)
library(tidyverse)

counts <- read.table("counts.txt", header = TRUE, sep = "\t")


# Going through the process for hour-0 samples
# Defining how many individuals from each of the Central Valley and San Pablo Populations there are, in the data
treatment_0 <- c(rep('C', 4), rep('S', 4))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_0 <- cbind(counts[,c(1:6)], counts[,c(7:10)], counts[,c(23:26)])
dx_0 <- dx_0[rowSums(dx_0[,c(7:14)]==0)<4,]

# Make exon ID
eid_0 = paste0(dx_0$Chr, ":", dx_0$Start)
# Keep gene IDs available
gene_id_0 <- dx_0$Chr

# Create the DGE list data set, then calculate normalization factors
dx_0 <- DGEList(dx_0[,7:14])
dx_0 <- calcNormFactors(dx_0, group = treatment_0)

# Define design matrix
design_0 <- model.matrix(~treatment_0)

# Voom transform
vx_0 <- voom(dx_0, design_0)

# Fit with limma
fx_0 <- lmFit(vx_0, design_0)

ex_0 <- diffSplice(fx_0, geneid = gene_id_0, exonid = eid_0)

# Pull out every cluster's result, regardless of significance
res_0 <- topSplice(ex_0, number = length(gene_id_0))

# Plot splice
plotSplice(ex_0)

# Write a results table
write.table(res_0[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_hour_0.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")


# Going through the process for hour-72 samples
# Defining how many individuals from each of the Central Valley and San Pablo Populations there are, in the data
treatment_72 <- c(rep('C', 6), rep('S', 6))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_72 <- cbind(counts[,c(1:6)], counts[,c(11:16)], counts[,c(27:32)])
dx_72 <- dx_72[rowSums(dx_72[,c(7:18)]==0)<6,]

# Make exon ID
eid_72 = paste0(dx_72$Chr, ":", dx_72$Start)
# Keep gene IDs available
gene_id_72 <- dx_72$Chr

# Create the DGE list data set, then calculate normalization factors
dx_72 <- DGEList(dx_72[,7:18])
dx_72 <- calcNormFactors(dx_72, group = treatment_72)

# Define design matrix
design_72 <- model.matrix(~treatment_72)

# Voom transform
vx_72 <- voom(dx_72, design_72)

# Fit with limma
fx_72 <- lmFit(vx_72, design_72)

ex_72 <- diffSplice(fx_72, geneid = gene_id_72, exonid = eid_72)

# Pull out every cluster's result, regardless of significance
res_72 <- topSplice(ex_72, number = length(gene_id_72))

# Plot splice
plotSplice(ex_72)

# Write a results table
write.table(res_72[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_hour_72.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# Going through the process for hour-168 samples
# Defining how many individuals from each of the Central Valley and San Pablo Populations there are, in the data
treatment_168 <- c(rep('C', 6), rep('S', 6))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_168 <- cbind(counts[,c(1:6)], counts[,c(17:22)], counts[,c(33:38)])
dx_168 <- dx_168[rowSums(dx_168[,c(7:18)]==0)<6,]

# Make exon ID
eid_168 = paste0(dx_168$Chr, ":", dx_168$Start)
# Keep gene IDs available
gene_id_168 <- dx_168$Chr

# Create the DGE list data set, then calculate normalization factors
dx_168 <- DGEList(dx_168[,7:18])
dx_168 <- calcNormFactors(dx_168, group = treatment_168)

# Define design matrix
design_168 <- model.matrix(~treatment_168)

# Voom transform
vx_168 <- voom(dx_168, design_168)

# Fit with limma
fx_168 <- lmFit(vx_168, design_168)

ex_168 <- diffSplice(fx_168, geneid = gene_id_168, exonid = eid_168)

# Pull out every cluster's result, regardless of significance
res_168 <- topSplice(ex_168, number = length(gene_id_168))

# Plot splice
plotSplice(ex_168)

# Write a results table
write.table(res_168[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_hour_168.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# Going through the process for CV 72 versus 0 hour samples, looking at intra-population DEU
# Defining how many individuals from each of each individual is in the data
treatment_CV_3v0 <- c(rep('0', 4), rep('72', 6))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_CV_3v0 <- cbind(counts[,c(1:6)], counts[,c(7:10)], counts[,c(11:16)])
dx_CV_3v0 <- dx_CV_3v0[rowSums(dx_CV_3v0[,c(7:16)]==0)<4,]

# Make exon ID
eid_CV_3v0 = paste0(dx_CV_3v0$Chr, ":", dx_CV_3v0$Start)
# Keep gene IDs available
gene_id_CV_3v0 <- dx_CV_3v0$Chr

# Create the DGE list data set, then calculate normalization factors
dx_CV_3v0 <- DGEList(dx_CV_3v0[,7:16])
dx_CV_3v0 <- calcNormFactors(dx_CV_3v0, group = treatment_CV_3v0)

# Define design matrix
design_CV_3v0 <- model.matrix(~treatment_CV_3v0)

# Voom transform
vx_CV_3v0 <- voom(dx_CV_3v0, design_CV_3v0)

# Fit with limma
fx_CV_3v0 <- lmFit(vx_CV_3v0, design_CV_3v0)

ex_CV_3v0 <- diffSplice(fx_CV_3v0, geneid = gene_id_CV_3v0, exonid = eid_CV_3v0)

# Pull out every cluster's result, regardless of significance
res_CV_3v0 <- topSplice(ex_CV_3v0, number = length(gene_id_CV_3v0))

# Plot splice
plotSplice(ex_CV_3v0)

# Write a results table
write.table(res_CV_3v0[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_CV_hour_72v0.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")


# Going through the process for CV 168 versus 0 hour samples, looking at intra-population DEU
# Defining how many individuals from each of each individual is in the data
treatment_CV_7v0 <- c(rep('0', 4), rep('168', 6))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_CV_7v0 <- cbind(counts[,c(1:6)], counts[,c(7:10)], counts[,c(33:38)])
dx_CV_7v0 <- dx_CV_7v0[rowSums(dx_CV_7v0[,c(7:16)]==0)<4,]

# Make exon ID
eid_CV_7v0 = paste0(dx_CV_7v0$Chr, ":", dx_CV_7v0$Start)
# Keep gene IDs available
gene_id_CV_7v0 <- dx_CV_7v0$Chr

# Create the DGE list data set, then calculate normalization factors
dx_CV_7v0 <- DGEList(dx_CV_7v0[,7:16])
dx_CV_7v0 <- calcNormFactors(dx_CV_7v0, group = treatment_CV_7v0)

# Define design matrix
design_CV_7v0 <- model.matrix(~treatment_CV_7v0)

# Voom transform
vx_CV_7v0 <- voom(dx_CV_7v0, design_CV_7v0)

# Fit with limma
fx_CV_7v0 <- lmFit(vx_CV_7v0, design_CV_7v0)
ex_CV_7v0 <- diffSplice(fx_CV_7v0, geneid = gene_id_CV_7v0, exonid = eid_CV_7v0)

# Pull out every cluster's result, regardless of significance
res_CV_7v0 <- topSplice(ex_CV_7v0, number = length(gene_id_CV_7v0))

# Plot splice
plotSplice(ex_CV_7v0)

# Write a results table
write.table(res_CV_7v0[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_CV_hour_168v0.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")


# Going through the process for CV 168 versus 72 hour samples, looking at intra-population DEU
# Defining how many individuals from each of each individual is in the data
treatment_CV_7v3 <- c(rep('72', 6), rep('168', 6))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_CV_7v3 <- cbind(counts[,c(1:6)], counts[,c(17:22)], counts[,c(11:16)])
dx_CV_7v3 <- dx_CV_7v3[rowSums(dx_CV_7v3[,c(7:18)]==0)<6,]

# Make exon ID
eid_CV_7v3 = paste0(dx_CV_7v3$Chr, ":", dx_CV_7v3$Start)
# Keep gene IDs available
gene_id_CV_7v3 <- dx_CV_7v3$Chr

# Create the DGE list data set, then calculate normalization factors
dx_CV_7v3 <- DGEList(dx_CV_7v3[,7:18])
dx_CV_7v3 <- calcNormFactors(dx_CV_7v3, group = treatment_CV_7v3)

# Define design matrix
design_CV_7v3 <- model.matrix(~treatment_CV_7v3)

# Voom transform
vx_CV_7v3 <- voom(dx_CV_7v3, design_CV_7v3)

# Fit with limma
fx_CV_7v3 <- lmFit(vx_CV_7v3, design_CV_7v3)
ex_CV_7v3 <- diffSplice(fx_CV_7v3, geneid = gene_id_CV_7v3, exonid = eid_CV_7v3)

# Pull out every cluster's result, regardless of significance
res_CV_7v3 <- topSplice(ex_CV_7v3, number = length(gene_id_CV_7v3))

# Plot splice
plotSplice(ex_CV_7v3)

# Write a results table
write.table(res_CV_7v3[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_CV_hour_168v72.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")


# Going through the process for SP 72 versus 0 hour samples, looking at intra-population DEU
# Defining how many individuals from each of each individual is in the data
treatment_SP_3v0 <- c(rep('0', 4), rep('72', 6))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_SP_3v0 <- cbind(counts[,c(1:6)], counts[,c(23:26)], counts[,c(27:32)])
dx_SP_3v0 <- dx_SP_3v0[rowSums(dx_SP_3v0[,c(7:16)]==0)<4,]

# Make exon ID
eid_SP_3v0 = paste0(dx_SP_3v0$Chr, ":", dx_SP_3v0$Start)
# Keep gene IDs available
gene_id_SP_3v0 <- dx_SP_3v0$Chr

# Create the DGE list data set, then calculate normalization factors
dx_SP_3v0 <- DGEList(dx_SP_3v0[,7:16])
dx_SP_3v0 <- calcNormFactors(dx_SP_3v0, group = treatment_SP_3v0)

# Define design matrix
design_SP_3v0 <- model.matrix(~treatment_SP_3v0)

# Voom transform
vx_SP_3v0 <- voom(dx_SP_3v0, design_SP_3v0)

# Fit with limma
fx_SP_3v0 <- lmFit(vx_SP_3v0, design_SP_3v0)

ex_SP_3v0 <- diffSplice(fx_SP_3v0, geneid = gene_id_SP_3v0, exonid = eid_SP_3v0)

# Pull out every cluster's result, regardless of significance
res_SP_3v0 <- topSplice(ex_SP_3v0, number = length(gene_id_SP_3v0))

# Plot splice
plotSplice(ex_SP_3v0)

# Write a results table
write.table(res_SP_3v0[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_SP_hour_72v0.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")


# Going through the process for SP 168 versus 0 hour samples, looking at intra-population DEU
# Defining how many individuals from each of each individual is in the data
treatment_SP_7v0 <- c(rep('0', 4), rep('168', 6))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_SP_7v0 <- cbind(counts[,c(1:6)], counts[,c(23:26)], counts[,c(33:38)])
dx_SP_7v0 <- dx_SP_7v0[rowSums(dx_SP_7v0[,c(7:16)]==0)<4,]

# Make exon ID
eid_SP_7v0 = paste0(dx_SP_7v0$Chr, ":", dx_SP_7v0$Start)
# Keep gene IDs available
gene_id_SP_7v0 <- dx_SP_7v0$Chr

# Create the DGE list data set, then calculate normalization factors
dx_SP_7v0 <- DGEList(dx_SP_7v0[,7:16])
dx_SP_7v0 <- calcNormFactors(dx_SP_7v0, group = treatment_SP_7v0)

# Define design matrix
design_SP_7v0 <- model.matrix(~treatment_SP_7v0)

# Voom transform
vx_SP_7v0 <- voom(dx_SP_7v0, design_SP_7v0)

# Fit with limma
fx_SP_7v0 <- lmFit(vx_SP_7v0, design_SP_7v0)
ex_SP_7v0 <- diffSplice(fx_SP_7v0, geneid = gene_id_SP_7v0, exonid = eid_SP_7v0)

# Pull out every cluster's result, regardless of significance
res_SP_7v0 <- topSplice(ex_SP_7v0, number = length(gene_id_SP_7v0))

# Plot splice
plotSplice(ex_SP_7v0)

# Write a results table
write.table(res_SP_7v0[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_SP_hour_168v0.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")


# Going through the process for SP 168 versus 72 hour samples, looking at intra-population DEU
# Defining how many individuals from each of each individual is in the data
treatment_SP_7v3 <- c(rep('72', 6), rep('168', 6))
# Pulling out just the cluster information, and the exon counts for just the individuals we care about here
dx_SP_7v3 <- cbind(counts[,c(1:6)], counts[,c(27:32)], counts[,c(33:38)])
dx_SP_7v3 <- dx_SP_7v3[rowSums(dx_SP_7v3[,c(7:18)]==0)<6,]

# Make exon ID
eid_SP_7v3 = paste0(dx_SP_7v3$Chr, ":", dx_SP_7v3$Start)
# Keep gene IDs available
gene_id_SP_7v3 <- dx_SP_7v3$Chr

# Create the DGE list data set, then calculate normalization factors
dx_SP_7v3 <- DGEList(dx_SP_7v3[,7:18])
dx_SP_7v3 <- calcNormFactors(dx_SP_7v3, group = treatment_SP_7v3)

# Define design matrix
design_SP_7v3 <- model.matrix(~treatment_SP_7v3)

# Voom transform
vx_SP_7v3 <- voom(dx_SP_7v3, design_SP_7v3)

# Fit with limma
fx_SP_7v3 <- lmFit(vx_SP_7v3, design_SP_7v3)
ex_SP_7v3 <- diffSplice(fx_SP_7v3, geneid = gene_id_SP_7v3, exonid = eid_SP_7v3)

# Pull out every cluster's result, regardless of significance
res_SP_7v3 <- topSplice(ex_SP_7v3, number = length(gene_id_SP_7v3))

# Plot splice
plotSplice(ex_SP_7v3)

# Write a results table
write.table(res_SP_7v3[,c("GeneID","NExons","P.Value","FDR")],"VoomSplice_SP_hour_168v72.txt",col.names = TRUE,row.names=FALSE,quote=FALSE,sep="\t")

# Define groups
treatment = c(rep('C', 16), rep('S', 16))
sample = c('C.0.4B.1',	'C.0.4B.2',	'C.0.4B.4',	'C.0.4B.5',	'C.HD3.2',	'C.HD3.3',	'C.HD3.4',	'C.HD3.5', 'C.HD3.6', 'C.HD3.8',	'C.HD7.2',	'C.HD7.3',	'C.HD7.4',	'C.HD7.7',	'C.HD7.8',	'C.HD7.9',	'S.0.4B.1',	'S.0.4B.2',	'S.0.4B.3',	'S.0.4B.5',	'S.HD3.1',	'S.HD3.2',	'S.HD3.3',	'S.HD3.4',	'S.HD3.5',	'S.HD3.7',	'S.HD7.2',	'S.HD7.3',	'S.HD7.4',	'S.HD7.5',	'S.HD7.8',	'S.HD7.9')

# Make DGE list and normalize
dx <- DGEList(counts[,c(7:38)])
dx <- calcNormFactors(dx, group = treatment)

# Make exon ID
eid = paste0(counts$Chr, ":", counts$Start)

# Define design matrix
design <- model.matrix(~treatment)

# Voom transform
vx <- voom(dx, design)

# Fit with limma
fx <- lmFit(vx, design)

ex <- diffSplice(fx, geneid = counts$Chr, exonid = eid)
res <- topSplice(ex, number = "20")

# Plot splice
plotSplice(ex)
# This is an overall response over the 3 time points between the two populations. Need to do comparisons at the 0, 72, and 168 hour time points specifically.


# Following the dexseq pipeline for DEU as per the supertranscripts paper supplementary code: https://github.com/Oshlack/superTranscript_paper_code/blob/master/differential_isoform_usage/DEXY.R

counts_dex <- read.table("counts.txt", sep = "\t", header = TRUE, stringsAsFactors = F)

clusters <- counts_dex[,2]
counts_only <- counts_dex[,7:38]

# Convert data to integer
cc <- data.matrix(counts_only)
cc = round(cc)

# Breaking up the data into the 0, 72, and 168 hour time points
hour_0 <- cbind(cc[,1:4], cc[,17:20])
hour_0<- hour_0[apply(hour_0, 1, function(x) !all(x==0)),]
hour_0 <- hour_0[rowSums(hour_0==0)<4, ]
hour_72 <- cbind(cc[,5:10], cc[,21:26])
hour_168 <- cbind(cc[,11:16], cc[,27:32])

# Unique exon names, using gene names and the starting base pair for the exon
exon_ids <- paste0(counts_dex[,1], "_", counts_dex[,3], "_", counts_dex[,4])

# Install DEXseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(DEXSeq)
library(Matrix)
library(GenomicRanges)
library(rtracklayer)
metadata <- read.delim("splittail_metadata.txt")

# Breaking up the metadata into the 3 different exposure times
colnames(metadata)[1] <- "sample"
colnames(metadata)[2] <- "condition"
metadata_0 <- metadata[metadata$Exposure_Time == 0,]
metadata_72 <- metadata[metadata$Exposure_Time == 3,]
metadata_168 <- metadata[metadata$Exposure_Time == 7,]


# Creating a hour 0 dexseq data set
dxd_0 <- DEXSeqDataSet(hour_0, design = ~sample + exon + condition:exon, featureID = as.factor(exon_ids), groupID = as.factor(clusters), sampleData = metadata_0)

##Estimate size factors and dispersions DEXseq does this based on a negative bionmial distribution
dxd_0 <- estimateSizeFactors(dxd_0)
dxd_0 <- estimateDispersions(dxd_0)

#Test for DEU
deu_0 <- testForDEU(dxd_0)

#Extract the results
dex_res_0 <- DEXSeqResults(deu_0)

#Get p-value per gene based on the p-values for all exons in the gene/cluster
pgq_0 <- perGeneQValue(dex_res_0, p = "pvalue")

## Save results to a text file and R object
save(deu, dex_res, pgq, file = "DEXY.Rdata")
tmp <- cbind(gene = names(pgq), "adjP" = pgq)
write.table(tmp, file = "DEXY.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Now switching over to the DEXSeq vignette: https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf
# How many exons are significant?


# How many genes do these exons represent?
table ( tapply( dex_res$padj < 0.1, dex_res$groupID, any ) )
#FALSE  TRUE 
#14632  5864 

# Visualizing dispersion
plotDispEsts(dxd)

# Estimating exon fold changes
exon_fold_change <- estimateExonFoldChanges(deu, fitExpToVar = "site")


# Trying exon fold change results
exon_res <- DEXSeqResults(exon_fold_change)


# An exploratory MA plot
DEXSeq::plotMA(exon_res)


# Attempting to make my own GFF using information from counts_dex, in an effort to plot the data
write.table(cbind(exon_ids, rep("SuperTranscripts", 490417), rep("exon", 490417), counts_dex$Start, counts_dex$End, rep(".", 490417), rep(".", 490417), rep("0", 490417), rep(".", 490417)), file = "dexseq.gff", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Importing my GFF for sending plot information to DEXSeq

gffRangedData<-import.gff("dexseq.gff")
myGranges<-as(gffRangedData, "GRanges")

# Adding in the genomic range data to the DEX seq models
exon_res$genomicData <- myGranges
dex_res$genomicData <- myGranges

# With length info included I can now plot the data
DEXSeqHTML(exon_res, FDR = 0.1, fitExpToVar = "site")

# Plotting the gene with the exon with the lowest adjusted p value between sites
plotDEXSeq(exon_res, exon_res$groupID[which.min(exon_res$padj)], legend=TRUE, fitExpToVar = "site", cex.axis=1.2, cex=1.3, lwd=2)

# Plotting the gene with the exon with the highest log fold change between the Red and Dauphin Rivers
plotDEXSeq(exon_res, exon_res$groupID[which.max(exon_res$log2fold_Red_River_Dauphin_River)], legend=TRUE, fitExpToVar = "site", cex.axis=1.2, cex=1.3, lwd=2)

#Get p-value per gene based on the p-values for all exons in the gene/cluster, and write this info to a table
pgq2 <- perGeneQValue(exon_res, p = "pvalue")
write.table(cbind(gene = names(pgq2), "Adj_P" = pgq2), file = "DEXY_per_gene_qvalues.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(cbind(gene = as.character(exon_res$groupID), exon = as.character(exon_res$featureID), RedvDauphin = exon_res$log2fold_Red_River_Dauphin_River, MathvDauphin = exon_res$log2fold_Matheson_Dauphin_River), file = "DEXY_exon_logfoldchanges.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")