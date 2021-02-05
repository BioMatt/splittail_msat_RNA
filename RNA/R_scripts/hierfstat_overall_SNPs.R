# Using hierfstat to look at Fst and some basic pop gen stats with the RNA seq LD-pruned HWE splittal SNPs
library(vcfR)
library(hierfstat)
library(adegenet)
library(tidyverse)

vcf <- read.vcfR("splittail_biallelic_maf0.05_snps_q30.recode.vcf")
vcf <- addID(vcf)


# Reformatting the vcfR object into a genlight object
splittail_snps <- vcfR2genlight(vcf)
splittail_snps

# Reformatting this to a matrix for input into hierfstat
splittail_snps <- as.matrix(splittail_snps)  

# Prepare 0, 1, 2 matrix in hierfstat format
# Use locations to test for population differentiation
splittail_snps[splittail_snps==0] <- 11
splittail_snps[splittail_snps==1] <- 12
splittail_snps[splittail_snps==2] <- 22

splittail_snps <- as.data.frame(cbind(as.numeric(c(rep("1", 16), rep("2", 16))), splittail_snps))
splittail_snps[,1:3]

# Writing an fstat file for use with genepop and NeEstimator
# NeEstimator V2 cannot take missing data, so removing those SNPs first
no_missing <- as.data.frame(t(drop_na(as.data.frame(t(splittail_snps)))))
write.fstat(no_missing, "splittail_allSNPs_nomissing_fstat.dat")

# Taking a look at Weir & Cockerham's Fst with Hierfstat
fst <- wc(splittail_snps, diploid = TRUE)
fst$FST

# Nei's Fst
nei_fst <- pairwise.neifst(splittail_snps)
nei_fst

# Generating confidence intervals for Fst over 1000 iterations
boot_fst <- boot.ppfst(splittail_snps, nboot = 1000)
boot_fst$ll
boot_fst$ul

# Generating confidence intervals for Fis over 1000 iterations
boot_fis <- boot.ppfis(splittail_snps, nboot = 1000)
boot_fis$fis.ci

# Taking a look at coalescent Fst with betas
all_betas <- betas(splittail_snps, nboot = 1000)
all_betas$betaiovl
all_betas$ci

# Looking at basic stats 
hier_stats <- basic.stats(splittail_snps)
hier_stats$overall

# Pulling out heterozygosity per locus from the basic stats
heterozygosity <- hier_stats$Ho
# Double checking that the average heterozygosity from both populations matches the value given by hier_stats$overall. It does!
mean(c(heterozygosity[,1], heterozygosity[,2]))

# Taking a look at population specific heterozygosities
cv_hetero <- mean(heterozygosity[,1])
sp_hetero <- mean(heterozygosity[,2])

# Following the same process but for gene diversity
# Pulling out diversity per locus from the basic stats
diversity <- hier_stats$Hs
# Double checking that the average diversity from both populations matches the value given by hier_stats$overall. It does!
mean(c(diversity[,1], diversity[,2]))

# Taking a look at population specific heterozygosities
cv_diversity <- mean(diversity[,1])
sp_diversity <- mean(diversity[,2])

# Calculating Fis based on Ho and Hs stats from each pop
# Taking a look at population specific heterozygosities
cv_inbreeding <- 1 - (cv_hetero/cv_diversity)
sp_inbreeding <- 1 - (sp_hetero/sp_diversity)

###########################################################
# Re-running the stats for the dataset with no missing values
# Taking a look at Weir & Cockerham's Fst with Hierfstat
nomissing_fst <- wc(no_missing, diploid = TRUE)
nomissing_fst$FST

# Generating confidence intervals for Fst over 1000 iterations
no_missing_boot_fst <- boot.ppfst(no_missing, nboot = 1000)
no_missing_boot_fst$ll
no_missing_boot_fst$ul

# Looking at basic stats 
nomissing_stats <- basic.stats(no_missing)
nomissing_stats$overall


# Generating confidence intervals for fis
no_missing_boot_fis <- boot.ppfis(no_missing, nboot = 1000)



# Running a PCA with adegenet
genlight <- vcfR2genlight(vcf)

pca <- glPca(genlight)
scatter(pca)
# Creating a fraction for variance explained by each principal component
snps_var_frac <- pca$eig/sum(pca$eig)
snps_var_frac[1] * 100
snps_var_frac[2] * 100

# Saving the PCA to rdata, for combination with PCAs from other datasets as a figure
save(pca, snps_var_frac, file = "Overall_SNP_pca.RData")
