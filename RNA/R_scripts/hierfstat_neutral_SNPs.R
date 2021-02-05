# Using hierfstat to look at Fst and some basic pop gen stats with the RNA seq LD-pruned HWE splittal SNPs
library(vcfR)
library(hierfstat)
library(adegenet)
library(tidyverse)

vcf <- read.vcfR("splittail_HWE_LD_prune2.vcf")


# Add unique IDs to make Adegenet happy. addID does not work for some reason.
vcf@fix[, "ID"] <- paste0(vcf@fix[,"CHROM"], "_", vcf@fix[, "POS"])

# Reformatting the vcfR object into a genlight object
splittail_snps <- vcfR2genlight(vcf)
splittail_snps

# Pulling individual names before continuing with analyses, then creating a table for PGD spider to make a genepop file from the VCF that includes pop info
ind.names <- splittail_snps$ind.names
ind.names <- str_split(ind.names, "_", 2, simplify = TRUE)[,2]
ind.pops <- cbind(ind.names, c(rep("Central_Valley", 16), rep("San_Pablo", 16)))

# Writing this table to a txt file for PGD spider to use for population definitions
write.table(ind.pops, "PGD_spider_splittail_pops.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Reformatting this to a matrix for input into hierfstat
splittail_snps <- as.matrix(splittail_snps)  

# Prepare 0, 1, 2 matrix in hierfstat format
# Use locations to test for population differentiation
splittail_snps[splittail_snps==0] <- 11
splittail_snps[splittail_snps==1] <- 12
splittail_snps[splittail_snps==2] <- 22

splittail_snps <- as.data.frame(cbind(as.numeric(c(rep("1", 16), rep("2", 16))), splittail_snps))
splittail_snps[,1:3]

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


# Writing an fstat file for use with genepop and NeEstimator
write.fstat(splittail_snps, "splittail_neutral_fstat.dat")

# Writing a Structure file
write.struct(splittail_snps, ilab = TRUE, pop=FALSE, fname = "splittail_neutral_snps.str")

# Reading the Structure file in, then adding individual names to it
str_file <- read_delim("splittail_neutral_snps.str", delim = " ")
# Remove the pop info, since we want Structure to run a-priori
str_file <- str_file[,-1]
str_names <- rep(ind.names, 1, each = 2)
str_names <- as_tibble(str_names)
str_file_ind.names <- cbind(str_names, str_file)
str_file_ind.names[,1:3]
write_delim(str_file_ind.names, "str_file_ind.names.str", col_names = TRUE)


# Looking at basic stats 
hier_stats <- basic.stats(splittail_snps)
hier_stats$overall

# Generating confidence intervals for overall Fis over 1000 iterations
overall_boot_fis <- boot.ppfis(cbind(rep(1, 32), splittail_snps[,-1]), nboot = 1000)
overall_boot_fis$fis.ci

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

# Taking a look at coalescent Fst with betas
neutral_betas <- betas(splittail_snps, nboot = 1000)
neutral_betas$betaiovl
neutral_betas$ci

# Running a PCA with adegenet
genlight <- vcfR2genlight(vcf)

pca <- glPca(genlight)
scatter(pca)
# Creating a fraction for variance explained by each principal component
snps_var_frac <- pca$eig/sum(pca$eig)
snps_var_frac[1] * 100
snps_var_frac[2] * 100

# Saving the PCA to rdata, for combination with PCAs from other datasets as a figure
save(pca, snps_var_frac, file = "LD_SNP_pca.RData")
