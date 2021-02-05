# Using SNPrelate to LD prune splittail SNPs by cluster
#if (!requireNamespace("BiocManager", quietly=TRUE))
  #install.packages("BiocManager")
#BiocManager::install("SNPRelate")
#BiocManager::install("gdsfmt")
#BiocManager::install("SeqArray")

library(gdsfmt)
library(SNPRelate)
library(tidyverse)
library(SeqArray)



vcf.fn <- "splittail_biallelic_snps_noNA_HWE_maf0.05.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "splittail.gds", method="biallelic.only")

snpgdsSummary("splittail.gds")
snpgdsClose("splittail.gds")

genofile <- SNPRelate::snpgdsOpen("splittail.gds")

# LD Pruning
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, autosome.only = FALSE, ld.threshold = 0.2)

snpset.id <- unlist(snpset)
# Write out the list of SNPs found to be out of LD
write_csv(as.data.frame(snpset.id), "LD_snp_list.csv")


snpgdsCreateGenoSet("splittail.gds", "LD_prune.gds", snp.id = snpset.id, verbose = TRUE)
snpgdsClose(genofile)

LD_prune_genofile <- SNPRelate::snpgdsOpen("LD_prune.gds")
snpgdsSummary("LD_prune.gds")

pca <- snpgdsPCA(LD_prune_genofile, num.thread=2, autosome.only = FALSE)
pc.percent <- pca$varprop*100
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

snpgdsOption(LD_prune_genofile, autosome.start = 1, autosome.end = 69951)
snpgdsSummary(LD_prune_genofile)


# Trying seq array, converting the SNP GDS to a seq array object, then to a VCF. Worked!
# Write the ped file with the pruned SNPs. Not super good, works but difficult to work with the map file in PGD spider.
snpgdsGDS2PED(LD_prune_genofile, ped.fn = "HWE_LD_splittail")

snpgdsClose(LD_prune_genofile)

# Rewriting the map file because it has one extra column and misses chromosome information for PGD spider. Calling all snp's as belonging to chromosome "1" since these are really all genes
map_file <- read_tsv("HWE_LD_splittail.map", col_names = FALSE)
map_file <- map_file %>% 
  select(-X2)
map_file <- map_file %>% 
  add_column(chromosome = c(rep("1", nrow(map_file))), .before = 1)

write_tsv(map_file, "LD_prune_snps_reformat.map", col_names = FALSE)

ped_file <- read_tsv("HWE_LD_splittail.ped", col_names = FALSE)

# With the new .map file and the original .ped file, use PGD spider to conver these into a LD-pruned vcf. Used version 2.1.5 here