# Script to use PCAdapt with the Splittail data, to look for signatures of local adaptation

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue")

library(pcadapt)
library(vcfR)
library(adegenet)
library(qvalue)
library(tidyverse)

# Taking a look at the VCF to see how the individuals are ordered
vcf_splittail <- read.vcfR("splittail_biallelic_maf0.05_snps_q30.recode.vcf")
vcf_splittail <- vcfR2genlight(vcf_splittail)
ind_names <- vcf_splittail@ind.names


splittail_pcadapt <- read.pcadapt("splittail_biallelic_maf0.05_snps_q30.recode.vcf", type = "vcf")

# Using PCAdapt with 20 dimensions in an exploratory capacity
exploratory_pcadapt <- pcadapt(input = splittail_pcadapt, K = 20)
plot(exploratory_pcadapt, option = "screeplot")

poplist.names <- c(rep("Central Valley", 16),rep("San Pablo Bay", 16))

plot(exploratory_pcadapt, option = "scores", pop = poplist.names)

plot(exploratory_pcadapt, option = "scores", i = 3, j = 4, pop = poplist.names)

plot(exploratory_pcadapt, option = "qqplot", threshold = 0.1)

plot(exploratory_pcadapt , option = "manhattan")


# Based on the previous scree plot with K=20, it seems like 2 PCs are most appropriate because only PC1 explains pop structure
final_pcadapt <- pcadapt(input = splittail_pcadapt, K = 2)

plot(final_pcadapt, option = "screeplot")

plot(final_pcadapt, option = "scores", pop = poplist.names)

# Exploring the data a bit
plot(final_pcadapt, option = "qqplot", threshold = 0.1)

plot(final_pcadapt, option = "manhattan")

# How much variance is explained by each PC
final_pcadapt$singular.values

# Applying q values to correct for multiple tests
qval <- qvalue(final_pcadapt$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)

# Relating outliers with principle components to get evolutionary pressure
snp_pc <- get.pc(final_pcadapt, outliers)


# Getting quick counts of how many outliers belong in PC1 and PC2
nrow(subset(snp_pc))
nrow(subset(snp_pc, PC == "1"))
nrow(subset(snp_pc, PC == "2"))

# Adding in locus names to the outlier SNPs and principle components
outlier.snp.names <- matrix(nrow = nrow(snp_pc), ncol = 1, dimnames = list(c(), c("SNP_ID")))
# Pulling the SNP ID's from the genlight object based on which SNP they are in PC Adapt
for (i in 1:nrow(snp_pc)){
  outlier.snp.names[i] <- vcf_splittail$loc.names[snp_pc$SNP[i]]
} 

# Put the SNP identifiers and PCs together
snp_pc <- cbind(outlier.snp.names, snp_pc)
rm(outlier.snp.names)

# Add in q values
snp_pc <- cbind(snp_pc, q.value = qval[outliers])
# Split the SNP ID's to get cluster ID's, for relating to gene annotations
snp_pc <- snp_pc %>% 
  tidyr::separate(SNP_ID, into = c("Cluster", "Position"), sep = "_", remove = FALSE)

# Removing the now unhelpful SNP index from the VCF, since we now have cluster and position information
snp_pc <- select(snp_pc, -SNP)

# Writing out the results to a tsv for downstream analysis
write.table(snp_pc, "splittail_PCadapt.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


##########################################################################################################
# Using Rtsne for fun
library(Rtsne)
tsne <- Rtsne(LD_pruned_matrix, dims = 2, initial_dims = 50, perplexity = 15, theta = 0.5, check_duplicates = FALSE, pca = TRUE, max_iter = 5000, verbose = TRUE)

tsne_metadata <- read_delim("wall_samples_tsne.txt", col_names = TRUE, delim = "\t")
View(as.matrix(rownames(LD_pruned_matrix)))
tsne_metadata <- full_join(tibble::as_tibble(as.matrix(rownames(LD_pruned_matrix))), tsne_metadata, by = c("V1" = "sample"))

plot(tsne$Y)
ggplot(as.data.frame(tsne$Y), aes(x=tsne$Y[,1], y=tsne$Y[,2], color = tsne_metadata$Location)) + 
  geom_point(aes(shape = as.factor(tsne_metadata$Year))) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  labs(color = "Site", shape = "Year Collected") +
  scale_colour_discrete(name = "Site", breaks = c("Dauphin_River", "Matheson", "Red_River"), labels = c("Dauphin River", "Matheson", "Red River"))
  
