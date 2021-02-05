# Install DEBrowser
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("debrowser", version = "3.8")

# This script is pulled from the main walleye differential expression script, and just shows how I read in the corset data.
library(tximport)
library(tidyverse)
library(BiocParallel)
# Set the number of cores for various functions to use, using parallel=TRUE in the function later
register(SnowParam(4))
library(apeglm)
library(DESeq2)
library(pheatmap)
library(edgeR)
#####################################################################
# This section involves importing Corset output data, adding in some metadata about site and year, then creating the DESeq2 table.


# Read in a table with my sample information and other metadata about them
samples <- read.table(file.path(getwd(), "splittail_metadata.txt"), header = TRUE, sep = "\t")


# Reading in a csv relating transcript ID's to gene ID's
#tx2gene <- read_tsv(file.path(getwd(), "splittail-clusters.txt"), col_names = FALSE)

# Reading the corset gene counts
corset_counts <- read_tsv(file.path(getwd(), "splittail-counts.txt"), col_names = TRUE)
corset_counts <- column_to_rownames(corset_counts, var = "gene")

# Adding rownames in to the metadata sample data frame
#rownames(samples) <- colnames(corset_counts)
# Turning the Exposure Time (in Days) into a factor for DESeq2
#samples$Exposure_Time <- as.factor(samples$Exposure_Time)
# Add in a combined variable of population + exposure time, calling it group, for easier contrasts
#samples$Group <- factor(paste0(population, "_", exposure), levels = c("C_0", "C_3", "C_7", "S_0", "S_3", "S_7"))
  
# Creating the DeSeq2 dataset
#corset.dds <- DESeqDataSetFromMatrix(corset_counts, DataFrame(samples), ~ Group)



# Find differential expression using DESeq. This step takes ~20 minutes. Much faster when parallel = TRUE and cores set to 4 at the beginning of the script.
#corset.deseq <- DESeq(corset.dds, betaPrior = FALSE, parallel = TRUE)

# A short line to filter out genes with 1 or fewer fragments per million counts in fewer than 4 samples. Choosing 4 here because that's the minimum number in a treatment group
#corset.deseq <- corset.deseq[rowSums(fpm(corset.deseq)>=1)>=4]

# Running a PCA on the salmon options + corset DESeq2 model. 
#corset.vsd <- varianceStabilizingTransformation(corset.deseq, blind = TRUE, fitType = "parametric")

# Running a PCA on the variance stabilized data, looking at Population and Exposure Time in days.
#pcaData <- plotPCA(corset.vsd, intgroup=c("Population", "Exposure_Time"), returnData=TRUE)
#percentVar <- round(100 * attr(pcaData, "percentVar"))
#corset.pca <- ggplot(pcaData, aes(PC1, PC2, color=Exposure_Time, shape=Population)) +
  #geom_point(size=3) +
  #xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  #ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  #ggtitle("PCA of DESeq2 Model from Corset")
#corset.pca

# Extract results from the deseq analysis
#C_v_S <- results(corset.deseq, name = "Population_S_vs_C", alpha = 0.05)

# Trying EdgeR because DESeq2 does not return useful comparisons
# Taken from the edgeR user's manual: https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
data_edger <- DGEList(counts = corset_counts)


# Create my design matrix for running the EdgeR model
population <- factor(samples$Population, levels = c("C", "S"))
exposure <- factor(samples$Exposure_Time, levels = c("0", "3", "7"))
group <- as.matrix(factor(paste0(population, "_", exposure), levels = c("C_0", "C_3", "C_7", "S_0", "S_3", "S_7")))
design <- model.matrix(~ 0 + group)

# Filter the count data for genes with any expression
keep <- filterByExpr(data_edger, design = design)
table(keep)
data_edger <- data_edger[keep,]


# Estimate dispersion in the data
data_edger <- estimateDisp(data_edger, design)

# Making some contrasts between the two populations at each time point
# Also doing intra-population comparisons, to find any gene showing plasticity
my.contrasts <- makeContrasts(CVvSP.0 = groupC_0 - groupS_0,
                              CVvSP.3 = groupC_3 - groupS_3,
                              CVvSP.7 = groupC_7 - groupS_7,
                              CV.3v0 = groupC_3 - groupC_0,
                              CV.7v0 = groupC_7 - groupC_0,
                              CV.7v3 = groupC_7 - groupC_3,
                              SP.3v0 = groupS_3 - groupS_0,
                              SP.7v0 = groupS_7 - groupS_0,
                              SP.7v3 = groupS_7 - groupS_3,
                              levels = design)

# Use a GLM with quasi-likelihoods to test for differential expression
quasi_glm_fit <- glmQLFit(data_edger, design)
plotQLDisp(quasi_glm_fit)

# Pulling results for the CV vs SP 0 hour treatment
quasi_glm_test_0 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"CVvSP.0"])
topTags(quasi_glm_test_0, p.value = 0.05)
write.table(topTags(quasi_glm_test_0, n=Inf, p.value = 0.05), "edgeR_DGE_hour_0.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Pulling results for the CV vs SP 72 hour treatment
quasi_glm_test_72 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"CVvSP.3"])
topTags(quasi_glm_test_72, p.value = 0.05)
write.table(topTags(quasi_glm_test_72, n=Inf, p.value = 0.05), "edgeR_DGE_hour_72.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Pulling results for the CV vs SP 168 hour treatment
quasi_glm_test_168 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"CVvSP.7"])
topTags(quasi_glm_test_168, p.value = 0.05)
write.table(topTags(quasi_glm_test_168, n=Inf, p.value = 0.05), "edgeR_DGE_hour_168.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Pulling results for the CV 72 hour versus 0 hour treatment
quasi_glm_test_CV3v0 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"CV.3v0"])
topTags(quasi_glm_test_CV3v0, p.value = 0.05)
write.table(topTags(quasi_glm_test_CV3v0, n=Inf, p.value = 0.05), "edgeR_DGE_CV_72v0_hour.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Pulling results for the CV 168 hour versus 0 hour treatment
quasi_glm_test_CV7v0 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"CV.7v0"])
topTags(quasi_glm_test_CV7v0, p.value = 0.05)
write.table(topTags(quasi_glm_test_CV7v0, n=Inf, p.value = 0.05), "edgeR_DGE_CV_168v0_hour.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Pulling results for the CV 168 hour versus 72 hour treatment
quasi_glm_test_CV7v3 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"CV.7v3"])
topTags(quasi_glm_test_CV7v3, p.value = 0.05)
write.table(topTags(quasi_glm_test_CV7v3, n=Inf, p.value = 0.05), "edgeR_DGE_CV_168v72_hour.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Pulling results for the SP 72 hour versus 0 hour treatment
quasi_glm_test_SP3v0 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"SP.3v0"])
topTags(quasi_glm_test_SP3v0, p.value = 0.05)
write.table(topTags(quasi_glm_test_SP3v0, n=Inf, p.value = 0.05), "edgeR_DGE_SP_72v0_hour.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Pulling results for the SP 168 hour versus 0 hour treatment
quasi_glm_test_SP7v0 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"SP.7v0"])
topTags(quasi_glm_test_SP7v0, p.value = 0.05)
write.table(topTags(quasi_glm_test_SP7v0, n=Inf, p.value = 0.05), "edgeR_DGE_SP_168v0_hour.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Pulling results for the SP 168 hour versus 72 hour treatment
quasi_glm_test_SP7v3 <- glmQLFTest(quasi_glm_fit, contrast = my.contrasts[,"SP.7v3"])
topTags(quasi_glm_test_SP7v3, p.value = 0.05)
write.table(topTags(quasi_glm_test_SP7v3, n=Inf, p.value = 0.05), "edgeR_DGE_SP_168v72_hour.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")