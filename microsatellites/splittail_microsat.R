# Script to play with reading in microsatellite data and put it into hierfstat format
setwd("C:/Users/mattt/Dropbox/splittail_redone/microsatellites")
library(hierfstat)
library(tidyverse)
library(adegenet)
library(pegas)

# Read in the data
microsats <- read_delim("splittail_msat_data.txt", delim = "\t", col_names = FALSE, trim_ws = TRUE)
metadata <- read_delim("splittail_msat_metadata.txt", delim = "\t", col_names = TRUE, trim_ws = TRUE)

rownames(microsats) <- metadata$Sample
metadata %>% 
  group_by(region, location) %>% 
  tally()

# Create an empty matrix with half the columns of the original data, for combining diploid microsats
hierf_microsats <- matrix(data = NA, nrow = nrow(microsats), ncol = ncol(microsats)/2)

# Setting rownames for individual names, and column names to the various microsatellite loci
rownames(hierf_microsats) <- rownames(microsats)
colnames(hierf_microsats) <- c(as.vector(t(read.table("microsat_list.txt"))))

# Combine the microsats, first going by rows, then columns
for (i in 1:nrow(hierf_microsats)){
  for (j in 1:ncol(hierf_microsats)){
    hierf_microsats[i,j] <- as.character(paste0(microsats[i, (j*2)-1], microsats[i, (j*2)]))
  }
  rm(i,j)
}


# Replace NA with 0 in the Hierfstat data frame
hierf_microsats <- replace(hierf_microsats, hierf_microsats == "missingmissing", NA)

# Adding in population info, 1 for San Pablo and 2 for Central Valley fish
Pop <- as.character(as.numeric(as.factor(metadata$region)))

hierf_microsats <- cbind(Pop, hierf_microsats)
hierf_microsats <- as.data.frame(hierf_microsats)

ind.names <- rownames(hierf_microsats)

# Convert that data frame factors into numeric characters for Hierfstat. Add individual names back in afterward since they go away after converting everything to numeric values
hierf_microsats <- mutate_all(hierf_microsats, function(x) as.numeric(as.character(x)))
rownames(hierf_microsats) <- ind.names

# Trying to keep individuals with 80% of microsatellites present. Out of 19 loci, that's 16 present. So we drop rows with more than 3  NA values
hierf_microsats <- hierf_microsats[rowSums( is.na(hierf_microsats) ) <= 3, ]


# Basic pop gen stats
basic_stats <- basic.stats(hierf_microsats, diploid = TRUE)
basic_stats
# Looking at Fst
fst <- wc(hierf_microsats, diploid = TRUE)
fst$FST




# Writing a structure file
write.struct(hierf_microsats, ilab = TRUE, pop = NULL, fname = "splittail_msats.str")
# Reading in the file to take out population labels, because write.struct can't not write those labels and does not write individual labels either. To the structure file add "pop " to the first space in the first line so read_delim can read it in properly
str_file <- read_delim("splittail_msats.str", delim = " ")
# Remove the pop column
str_file <- select(str_file, -pop)

# Create a table of duplicated individual names, for adding to the structure input
str_ind_names <- as.data.frame(rownames(hierf_microsats))
str_ind_names <- str_ind_names[rep(seq_len(nrow(str_ind_names)), each=2),]

# Add those names to the structure file
str_file <- cbind(str_ind_names, str_file)

# Delete the unneeded duplicated individual names, write the structure file for input into Structure. Be sure to delete the first column name
rm(str_ind_names)
write_delim(str_file, "structure_no_pop.str", delim = " ")

# Write an fstat file, for reloading with adegenet
write.fstat(hierf_microsats, fname = "splittailmsats.dat")

# Read in the fstat file, converting it to genind format
genind_microsats <- import2genind("splittailmsats.dat")

# Trying a DAPC, using interactive mode to get a feel for the data
dapc1 <- dapc(genind_microsats)
scatter.dapc(dapc1)

# Trying find.clusters. Here, 3 clusters have been chosen which was reasonable given a steep drop in the Bayesian Information Criterion.
clusters <- find.clusters(genind_microsats, n.pca = 75)

dapc_clusters <- dapc(genind_microsats, clusters$grp, max.n.clust = 40)
scatter.dapc(dapc_clusters)

# Adding the cluster data to a new data frame to see how they stack up
cluster_metadata <- cbind(metadata, cluster = clusters$grp)
# Tallying up the number of individuals in each cluster, within "regions" as defined by Ken's paper to see where these fish are coming from
cluster_tally <- cluster_metadata %>% 
  group_by(cluster, region, location) %>% 
  dplyr::tally()

# Plotting the number of individuals in each river for each of the 3 clusters

barplot(subset(cluster_tally, cluster == 1)$n, names.arg = subset(cluster_tally, cluster == 1)$location, xlab = "River", ylab = "Number of individuals", main = "Cluster 1 Splittail", cex.names = 1)

barplot(subset(cluster_tally, cluster == 2)$n, names.arg = subset(cluster_tally, cluster == 2)$location, xlab = "River", ylab = "Number of individuals", main = "Cluster 2 Splittail")

barplot(subset(cluster_tally, cluster == 3)$n, names.arg = subset(cluster_tally, cluster == 3)$location, xlab = "River", ylab = "Number of individuals", main = "Cluster 3 Splittail")


# No real differences popped out between Clusters 2 and 3 except the Sacramento and San Joaquin Rivers switching in importance, so 2 clusters is probably more reasonable for these data. 


# Using find.clusters with 2 clusters
two_clusters <- find.clusters(genind_microsats, n.pca = 75, n.clust = 2)
dapc_2clusters <- dapc(genind_microsats, two_clusters$grp, n.pca = 75, max.n.clust = 40)
scatter.dapc(dapc_2clusters)

# Adding the cluster data to a new data frame to see how they stack up
two_cluster_metadata <- cbind(metadata, cluster = two_clusters$grp)

# Tallying up the number of individuals in each cluster, within "regions" as defined by Ken's paper to see where these fish are coming from
two_cluster_tally <- two_cluster_metadata %>% 
  group_by(cluster, region, location) %>% 
  dplyr::tally()
barplot(subset(two_cluster_tally, cluster == 1)$n, names.arg = subset(two_cluster_tally, cluster == 1)$location, xlab = "River", ylab = "Number of individuals", main = "Cluster 1 Splittail", cex.names = 1)

barplot(subset(two_cluster_tally, cluster == 2)$n, names.arg = subset(two_cluster_tally, cluster == 2)$location, xlab = "River", ylab = "Number of individuals", main = "Cluster 2 Splittail", cex.names = 1)


# Taking a look at the population assignment plots
assignplot(dapc_2clusters)
compoplot(dapc_2clusters, posi = "bottomright", txt.leg = paste("Cluster", 1:2), lab = "", ncol = 1, xlab = "Individuals")


# Taking a look at pairwise Fst using the two clusters with pegas
# Even though pairwise.fst is apparently a wrapper for Hierfstat Nei's pairwise Fst estimates, the value given here is ~1/2 of what the pairwise.neifst function returns from Hierfstat.
clustered_fst <- pairwise.fst(genind_microsats, pop = two_cluster_metadata$cluster)
clustered_fst

# Taking a look at Nei's Fst using two clusters and Hierfstat
two_clust_neifst <- pairwise.neifst(cbind(two_cluster_metadata$cluster, hierf_microsats[,-1]), diploid = TRUE)
two_clust_neifst

# Taking a look at pairwise Fst using the two clusters, but with Hierfstat to compare to pegas
two_clust_wcfst <- wc(cbind(two_cluster_metadata$cluster, hierf_microsats[,-1]), diploid = TRUE)
two_clust_wcfst$FST

# Creating a two cluster data frame, sorting it by population, then turning population into numeric for boot.ppfst to give Fst confidence intervals
two_clust_msats <- cbind(two_cluster_metadata$cluster, hierf_microsats[,-1])
two_clust_msats <- two_clust_msats[order(two_clust_msats$`two_cluster_metadata$cluster`),]
two_clust_msats <- as.data.frame(two_clust_msats)
two_clust_msats$`two_cluster_metadata$cluster` <- as.numeric(two_clust_msats$`two_cluster_metadata$cluster`)

# Generating confidence intervals for Fst over 1000 iterations
boot_fst <- boot.ppfst(two_clust_msats, nboot = 1000)
boot_fst$ll
boot_fst$ul



# Writing a .dat file using the reassigned msats for use with Ne Estimator V2
write.fstat(two_clust_msats, fname = "adegenet_reassigned_splittailmsats.dat")

# Using basic.stats from Hierfstat and the two reassigned clusters to look at...basic stats 
two_clust_basic.stats <- basic.stats(cbind(two_cluster_metadata$cluster, hierf_microsats[,-1]), diploid = TRUE)
two_clust_basic.stats$Ho
two_clust_basic.stats$Hs
two_clust_basic.stats$Fis
two_clust_basic.stats$overall

# Writing the two-cluster metadata out to a text file 
write_tsv(two_cluster_metadata, "splittail_msat_2clust.txt")


# Pulling out heterozygosity per locus from the basic stats
heterozygosity <- two_clust_basic.stats$Ho
# Double checking that the average heterozygosity from both populations matches the value given by hier_stats$overall. It does!
mean(c(heterozygosity[,1], heterozygosity[,2]))

# Taking a look at population specific heterozygosities
cv_hetero <- mean(heterozygosity[,1])
sp_hetero <- mean(heterozygosity[,2])

# Following the same process but for gene diversity
# Pulling out diversity per locus from the basic stats
diversity <- two_clust_basic.stats$Hs
# Double checking that the average diversity from both populations matches the value given by hier_stats$overall. It does!
mean(c(diversity[,1], diversity[,2]))

# Taking a look at population specific heterozygosities
cv_diversity <- mean(diversity[,1])
sp_diversity <- mean(diversity[,2])

# Calculating Fis based on Ho and Hs stats from each pop
# Taking a look at population specific heterozygosities
cv_inbreeding <- 1 - (cv_hetero/cv_diversity)
sp_inbreeding <- 1 - (sp_hetero/sp_diversity)


# Pulling 16 random individuals from each of cluster 1 (CV) and cluster 2 (SP), then outputting a fstat file for input into Ne Estimator V2. Doing this multiple times to test for how sample size might affect Ne estimates between microsats and RNA SNPs

# First, breaking apart the reassigned msat table into the two populations, Central Valley (CV) and San Pablo Bay (SP)
cv_msats <- two_clust_msats[two_clust_msats[, "two_cluster_metadata$cluster"] == 1,]
sp_msats <- two_clust_msats[two_clust_msats[, "two_cluster_metadata$cluster"] == 2,]

# A function to randomly sample a chosen number of individuals from each of the CV and SP microsatellite matrices, and write fstat input files for use with Ne Estimator V2
# Sampling is without replacement, here
write_fstat_random_msats <- function(iterations, sample_size) {
  for (i in 1:iterations) {
    temp_cv <- cv_msats[sample(nrow(cv_msats),size=sample_size,replace=FALSE),]
    temp_sp <- sp_msats[sample(nrow(sp_msats),size=sample_size,replace=FALSE),]
    
    temp_total <- rbind(temp_cv, temp_sp)
    
    file_name <- paste0(i, "_splittail_random_msats.dat")
    
    write.fstat(temp_total, file_name)
  }
}

# Using the function to write 40 input files with 16 individuals from each cluster, for 32 individuals total in each. Exactly matching the sample sizes in the RNA data set
write_fstat_random_msats(40, 16)

#################################################################################
# Writing Colony offspring input files from the Adegenet reassignments
cv_indivs <- rownames(cv_msats)
sp_indivs <- rownames(sp_msats)

# Pulling the original msat data, but only for the individuals in each of the population lists
cv_colony <- microsats[cv_indivs,]
# Add informative rownames, replace "missing" with 0
cv_colony <- cv_colony %>% mutate_all(~ replace(., . == "missing", 0))
rownames(cv_colony) <- cv_indivs


# Repeat the process for the San Pablo splittail
sp_colony <- microsats[sp_indivs,]
sp_colony <- sp_colony %>% mutate_all(~ replace(., . == "missing", 0))
rownames(sp_colony) <- sp_indivs

# Write the two datasets for Colony
write.table(cv_colony, "cv_colony_offspring.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)

write.table(sp_colony, "sp_colony_offspring.txt", sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
#################################################################################
# Full-sibling results from Colony. Using Inclusive Probabilities >0.80, otherwise the individuals are considered non-sibs
# First, the list of individuals in the San Pablo population
sp_fullsibs <- c("PET2011-106", "PET2011-255", "LO-PET2011-01", "NAPA2012-10", "LO-PET2011-10")
# Then the individuals in the Central Valley population
cv_fullsibs <- c("PET2011-96", "PET2011-55", "PET2011-12")

# Creating a hierfstat data frame without the sibling individuals
no_fullsibs <- subset(two_clust_msats, !(rownames(two_clust_msats) %in% c(sp_fullsibs, cv_fullsibs)))

# Looking at fst with the sibling-filtered data
nosib_wcfst <- wc(no_fullsibs, diploid = TRUE)
nosib_wcfst$FST

# Taking a look at Nei's Fst using two clusters and Hierfstat
nonsib_neifst <- pairwise.neifst(no_fullsibs, diploid = TRUE)
nonsib_neifst

# Generating confidence intervals for Fst over 1000 iterations
nosib_bootfst <- boot.ppfst(no_fullsibs, nboot = 1000)
nosib_bootfst$ll
nosib_bootfst$ul


# Using basic.stats from Hierfstat and the two reassigned clusters to look at...basic stats 
nosib_basic.stats <- basic.stats(no_fullsibs, diploid = TRUE)
nosib_basic.stats$Ho
nosib_basic.stats$Hs
nosib_basic.stats$Fis
nosib_basic.stats$overall

# Looking at the Fis confidence interval
nosib_bootfis <- boot.ppfis(no_fullsibs, nboot = 1000)
nosib_bootfis$fis.ci

# Pulling out heterozygosity per locus from the basic stats
heterozygosity <- nosib_basic.stats$Ho
# Double checking that the average heterozygosity from both populations matches the value given by hier_stats$overall. It does!
mean(c(heterozygosity[,1], heterozygosity[,2]))

# Taking a look at population specific heterozygosities
cv_hetero <- mean(heterozygosity[,1])
sp_hetero <- mean(heterozygosity[,2])

# Following the same process but for gene diversity
# Pulling out diversity per locus from the basic stats
diversity <- nosib_basic.stats$Hs
# Double checking that the average diversity from both populations matches the value given by hier_stats$overall. It does!
mean(c(diversity[,1], diversity[,2]))

# Taking a look at population specific heterozygosities
cv_diversity <- mean(diversity[,1])
sp_diversity <- mean(diversity[,2])


# Calculating Fis based on Ho and Hs stats from each pop
# Taking a look at population specific heterozygosities
cv_inbreeding <- 1 - (cv_hetero/cv_diversity)
sp_inbreeding <- 1 - (sp_hetero/sp_diversity)

# Writing a .dat file using the nonsib msats for use with Ne Estimator V2
# This version writes the .dat file with no missing data. Not used for Ne estimates
#write.fstat(no_fullsibs[complete.cases(no_fullsibs), ], fname = "nosibs_splittailmsats.dat")
# This version writes the Ne estimates with missing data
write.fstat(no_fullsibs, fname = "nosibs_splittailmsats.dat")

# Taking a look at coalescent Fst with betas
msat_betas <- betas(no_fullsibs, nboot = 1000)
msat_betas$betaiovl
msat_betas$ci

# Running a PCA
pca <- indpca(no_fullsibs)
plot(pca)

# Saving the PCA to rdata, for combination with PCAs from other datasets as a figure
save(pca, file = "Msat_pca.RData")