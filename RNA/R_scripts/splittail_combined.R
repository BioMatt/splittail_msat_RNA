# Script for combining SNPs and expression data in the Sacramento splittail re-analysis

library(tidyverse)

# Reading in gene annotations, keep only relevant rows (data is pre-filtered for E values), then rename for clarity
annotations <- read_delim("gene_ids.txt", col_names = TRUE, delim = "\t")
annotations <- select(annotations, c(transcript_id, sprot_Top_BLASTX_hit, 'Sprot gene symbols', Uniprot, 'Sprot Blast Hit Species'))
annotations <- rename(annotations, 
                      transcript = transcript_id,
                      sprot_function = sprot_Top_BLASTX_hit,
                      sprot_gene = 'Sprot gene symbols',
                      uniprot_gene = Uniprot,
                      sprot_species = 'Sprot Blast Hit Species')

# Reading in the Corset cluster to Trinity transcript translation file
clusters <- read_delim("splittail-clusters.txt", col_names = FALSE, delim = "\t")
clusters <- rename(clusters,
                   transcript = X1,
                   corset_cluster = X2)

# Putting the cluster and annotation info together
clusters_w_annotations <- left_join(clusters, annotations)
# Keep only the clusters with some known gene function
# Not doing this for counts to be most accurate
#clusters_w_annotations <- filter(clusters_w_annotations,
                                 #is.na(sprot_function) == FALSE,
                                 #sprot_function != 0)
# Keep only distinct clusters
clusters_w_annotations <- distinct(clusters_w_annotations, corset_cluster, .keep_all = TRUE)

# Reading in the PCAdapt SNPs
pcadapt_snps <- read_delim("splittail_PCadapt.txt", col_names = TRUE, delim = "\t")
pcadapt_snps <- rename(pcadapt_snps,
                       pcadapt_SNP_ID = SNP_ID,
                       corset_cluster = Cluster,
                       pcadapt_position = Position,
                       PC_adapt_loading = PC,
                       PC_adapt_qvalue = q.value)

# Taking only the SNPs along PC1, which corresponds to population/site for these data
pcadapt_snps_PC1 <- filter(pcadapt_snps, PC_adapt_loading == 1)

# Reading the Bayescan output
bayescan_output <- read_delim("splittail_bayescan_out_fst.txt", col_names = TRUE, delim = " ")
# Read the Bayescan SNP IDs
bayescan_snps <- read_delim("bayescan_SNPs.txt", col_names = FALSE, delim = "\t")
# Combine the data, filter for q < 0.05
bayescan_output <- cbind(bayescan_snps, bayescan_output)
bayescan_output <- rename(bayescan_output,
                          bayescan_SNP_ID = X1,
                          bayescan_qval = `  qval`,
                          bayescan_alpha = " alpha")
bayescan_output <- separate(bayescan_output, bayescan_SNP_ID, 
                            into = c("corset_cluster", "SNP_position"),
                            sep = "_",
                            remove = FALSE)

bayescan_output <- filter(bayescan_output, bayescan_qval < 0.05)

# Select only q values and the strength of selection as measured by alpha
bayescan_output <- select(bayescan_output, c("bayescan_SNP_ID", "corset_cluster", "bayescan_qval", "bayescan_alpha"))

# Reading in the differential gene expression data, filtering for logCPM over 0 and a False Discovery Rate of 0.05 or lower
CvS_DGE_hour0 <- read_delim("edgeR_DGE_hour_0.txt", col_names = TRUE, delim = "\t")
# CvS_DGE_hour0 <- rename(CvS_DGE_hour0,
#                         CvS_LFC_0 = logFC,
#                         CvS_FDR_0 = FDR)
# CvS_DGE_hour0 <- filter(CvS_DGE_hour0,
#                         CvS_FDR_0 <= 0.05)
# Nothing is significant!


CvS_DGE_hour72 <- read_delim("edgeR_DGE_hour_72.txt", col_names = TRUE, delim = "\t")
# CvS_DGE_hour72 <- rename(CvS_DGE_hour72,
#                         CvS_LFC_72 = logFC,
#                         CvS_FDR_72 = FDR)
# CvS_DGE_hour72 <- filter(CvS_DGE_hour72,
#                         CvS_FDR_72 <= 0.05)
# Nothing is significant!


CvS_DGE_hour168 <- read_delim("edgeR_DGE_hour_168.txt", col_names = TRUE, delim = "\t")
CvS_DGE_hour168 <- rename(CvS_DGE_hour168,
                         CvS_LFC_168 = logFC,
                         CvS_FDR_168 = FDR)
CvS_DGE_hour168 <- filter(CvS_DGE_hour168,
                          CvS_FDR_168 <= 0.05)
CvS_DGE_hour168 <- select(CvS_DGE_hour168, c(corset_cluster, CvS_LFC_168, CvS_FDR_168))


# Reading the differential exon usage data, filtering for FDR of 0.05 or lower
CvS_DEU_hour0 <- read_delim("VoomSplice_hour_0.txt", col_names = TRUE, delim = "\t")
CvS_DEU_hour0 <- rename(CvS_DEU_hour0,
                        corset_cluster = GeneID,
                        N.Exons_0 = NExons,
                        DEU_FDR_0 = FDR)
CvS_DEU_hour0 <- filter(CvS_DEU_hour0,
                        DEU_FDR_0 <= 0.05)
CvS_DEU_hour0 <- select(CvS_DEU_hour0, c(corset_cluster, N.Exons_0, DEU_FDR_0))

CvS_DEU_hour72 <- read_delim("VoomSplice_hour_72.txt", col_names = TRUE, delim = "\t")
CvS_DEU_hour72 <- rename(CvS_DEU_hour72,
                        corset_cluster = GeneID,
                        N.Exons_72 = NExons,
                        DEU_FDR_72 = FDR)
CvS_DEU_hour72 <- filter(CvS_DEU_hour72,
                         DEU_FDR_72 <= 0.05)
CvS_DEU_hour72 <- select(CvS_DEU_hour72, c(corset_cluster, N.Exons_72, DEU_FDR_72))

CvS_DEU_hour168 <- read_delim("VoomSplice_hour_168.txt", col_names = TRUE, delim = "\t")
CvS_DEU_hour168 <- rename(CvS_DEU_hour168, 
                         corset_cluster = GeneID,
                         N.Exons_168 = NExons,
                         DEU_FDR_168 = FDR)
CvS_DEU_hour168 <- filter(CvS_DEU_hour168,
                          DEU_FDR_168 <= 0.05)
CvS_DEU_hour168 <- select(CvS_DEU_hour168, c(corset_cluster, N.Exons_168, DEU_FDR_168))



# Looking at intra-population DGE
# Starting with Central Valley hour 72 vs hour 0
CV_DGE_hour72v0 <- read_delim("edgeR_DGE_CV_72v0_hour.txt", col_names = TRUE, delim = "\t")
CV_DGE_hour72v0 <- rename(CV_DGE_hour72v0,
                          CV_LFC_72v0 = logFC,
                          CV_FDR_72v0 = FDR)
CV_DGE_hour72v0 <- filter(CV_DGE_hour72v0,
                          CV_FDR_72v0 <= 0.05)
CV_DGE_hour72v0 <- select(CV_DGE_hour72v0, c(corset_cluster, CV_LFC_72v0, CV_FDR_72v0))

# Now Central Valley hour 168 vs hour 0
CV_DGE_hour168v0 <- read_delim("edgeR_DGE_CV_168v0_hour.txt", col_names = TRUE, delim = "\t")
CV_DGE_hour168v0 <- rename(CV_DGE_hour168v0,
                          CV_LFC_168v0 = logFC,
                          CV_FDR_168v0 = FDR)
CV_DGE_hour168v0 <- filter(CV_DGE_hour168v0,
                          CV_FDR_168v0 <= 0.05)
CV_DGE_hour168v0 <- select(CV_DGE_hour168v0, c(corset_cluster, CV_LFC_168v0, CV_FDR_168v0))

# Now Central Valley hour 168 vs hour 72
CV_DGE_hour168v72 <- read_delim("edgeR_DGE_CV_168v72_hour.txt", col_names = TRUE, delim = "\t")
CV_DGE_hour168v72 <- rename(CV_DGE_hour168v72,
                           CV_LFC_168v72 = logFC,
                           CV_FDR_168v72 = FDR)
CV_DGE_hour168v72 <- filter(CV_DGE_hour168v72,
                           CV_FDR_168v72 <= 0.05)
CV_DGE_hour168v72 <- select(CV_DGE_hour168v72, c(corset_cluster, CV_LFC_168v72, CV_FDR_168v72))

# Doing San Pablo hour 72 vs hour 0
SP_DGE_hour72v0 <- read_delim("edgeR_DGE_SP_72v0_hour.txt", col_names = TRUE, delim = "\t")
SP_DGE_hour72v0 <- rename(SP_DGE_hour72v0,
                          SP_LFC_72v0 = logFC,
                          SP_FDR_72v0 = FDR)
SP_DGE_hour72v0 <- filter(SP_DGE_hour72v0,
                          SP_FDR_72v0 <= 0.05)
SP_DGE_hour72v0 <- select(SP_DGE_hour72v0, c(corset_cluster, SP_LFC_72v0, SP_FDR_72v0))

# Now San Pablo hour 168 vs hour 0
SP_DGE_hour168v0 <- read_delim("edgeR_DGE_SP_168v0_hour.txt", col_names = TRUE, delim = "\t")
SP_DGE_hour168v0 <- rename(SP_DGE_hour168v0,
                           SP_LFC_168v0 = logFC,
                           SP_FDR_168v0 = FDR)
SP_DGE_hour168v0 <- filter(SP_DGE_hour168v0,
                           SP_FDR_168v0 <= 0.05)
SP_DGE_hour168v0 <- select(SP_DGE_hour168v0, c(corset_cluster, SP_LFC_168v0, SP_FDR_168v0))

# Now San Pablo hour 168 vs hour 72
SP_DGE_hour168v72 <- read_delim("edgeR_DGE_SP_168v72_hour.txt", col_names = TRUE, delim = "\t")
SP_DGE_hour168v72 <- rename(SP_DGE_hour168v72,
                            SP_LFC_168v72 = logFC,
                            SP_FDR_168v72 = FDR)
SP_DGE_hour168v72 <- filter(SP_DGE_hour168v72,
                            SP_FDR_168v72 <= 0.05)
SP_DGE_hour168v72 <- select(SP_DGE_hour168v72, c(corset_cluster, SP_LFC_168v72, SP_FDR_168v72))



# Reading intra-population DEU data, filtering for relevant variables and significance as before
# Starting with Central Valley population, and 72 versus 0 hour comparison at first
CV_DEU_72v0 <- read_delim("VoomSplice_CV_hour_72v0.txt", col_names = TRUE, delim = "\t")
CV_DEU_72v0 <- rename(CV_DEU_72v0,
                        corset_cluster = GeneID,
                        N.Exons_CV_72v0 = NExons,
                        DEU_FDR_CV_72v0 = FDR)
CV_DEU_72v0 <- filter(CV_DEU_72v0,
                      DEU_FDR_CV_72v0 <= 0.05)
CV_DEU_72v0 <- select(CV_DEU_72v0, c(corset_cluster, N.Exons_CV_72v0, DEU_FDR_CV_72v0))

# Reading the Central Valley 168 v 0 hour time point
CV_DEU_168v0 <- read_delim("VoomSplice_CV_hour_168v0.txt", col_names = TRUE, delim = "\t")
CV_DEU_168v0 <- rename(CV_DEU_168v0,
                      corset_cluster = GeneID,
                      N.Exons_CV_168v0 = NExons,
                      DEU_FDR_CV_168v0 = FDR)
CV_DEU_168v0 <- filter(CV_DEU_168v0,
                      DEU_FDR_CV_168v0 <= 0.05)
CV_DEU_168v0 <- select(CV_DEU_168v0, c(corset_cluster, N.Exons_CV_168v0, DEU_FDR_CV_168v0))
# Nothing is significant!

# Reading the Central Valley 72 v 168 hour time point
CV_DEU_72v168 <- read_delim("VoomSplice_CV_hour_72v168.txt", col_names = TRUE, delim = "\t")
CV_DEU_72v168 <- rename(CV_DEU_72v168,
                       corset_cluster = GeneID,
                       N.Exons_CV_72v168 = NExons,
                       DEU_FDR_CV_72v168 = FDR)
CV_DEU_72v168 <- filter(CV_DEU_72v168,
                       DEU_FDR_CV_72v168 <= 0.05)
CV_DEU_72v168 <- select(CV_DEU_72v168, c(corset_cluster, N.Exons_CV_72v168, DEU_FDR_CV_72v168))

# Now doing San Pablo intra-population DEU, starting with 72 v 0 hour time point comparison
SP_DEU_72v0 <- read_delim("VoomSplice_SP_hour_72v0.txt", col_names = TRUE, delim = "\t")
SP_DEU_72v0 <- rename(SP_DEU_72v0,
                      corset_cluster = GeneID,
                      N.Exons_SP_72v0 = NExons,
                      DEU_FDR_SP_72v0 = FDR)
SP_DEU_72v0 <- filter(SP_DEU_72v0,
                      DEU_FDR_SP_72v0 <= 0.05)
SP_DEU_72v0 <- select(SP_DEU_72v0, c(corset_cluster, N.Exons_SP_72v0, DEU_FDR_SP_72v0))

# Reading the San Pablo 168 v 0 hour time point
SP_DEU_168v0 <- read_delim("VoomSplice_SP_hour_168v0.txt", col_names = TRUE, delim = "\t")
SP_DEU_168v0 <- rename(SP_DEU_168v0,
                       corset_cluster = GeneID,
                       N.Exons_SP_168v0 = NExons,
                       DEU_FDR_SP_168v0 = FDR)
SP_DEU_168v0 <- filter(SP_DEU_168v0,
                       DEU_FDR_SP_168v0 <= 0.05)
SP_DEU_168v0 <- select(SP_DEU_168v0, c(corset_cluster, N.Exons_SP_168v0, DEU_FDR_SP_168v0))

# Reading the San Pablo 168 v 72 hour time point
SP_DEU_168v72 <- read_delim("VoomSplice_SP_hour_168v72.txt", col_names = TRUE, delim = "\t")
SP_DEU_168v72 <- rename(SP_DEU_168v72,
                       corset_cluster = GeneID,
                       N.Exons_SP_168v72 = NExons,
                       DEU_FDR_SP_168v72 = FDR)
SP_DEU_168v72 <- filter(SP_DEU_168v72,
                       DEU_FDR_SP_168v72 <= 0.05)
SP_DEU_168v72 <- select(SP_DEU_168v72, c(corset_cluster, N.Exons_SP_168v72, DEU_FDR_SP_168v72))


# Combining all DGE and DEU data sets together with just significant hits to get all genes showing plasticity. Nothing showed plasticity in terms of expression variation, so not including those data
# No DGE between C and S in hour 0 and hour 72, so only doing hour 168
plasticity <- left_join(clusters_w_annotations, CvS_DGE_hour168)

# Adding the intra-population DGE to the list. Starting with the intra-Central Valley comparisons
plasticity <- left_join(plasticity, CV_DGE_hour72v0)
plasticity <- left_join(plasticity, CV_DGE_hour168v0)
plasticity <- left_join(plasticity, CV_DGE_hour168v72)

# More intra-population DGE, with San Pablo this time
plasticity <- left_join(plasticity, SP_DGE_hour72v0)
plasticity <- left_join(plasticity, SP_DGE_hour168v0)
plasticity <- left_join(plasticity, SP_DGE_hour168v72)

# Adding in the intra-population DEU into the plasticity spreadsheet
plasticity <- left_join(plasticity, CV_DEU_72v0)
plasticity <- left_join(plasticity, CV_DEU_168v0)
plasticity <- left_join(plasticity, CV_DEU_72v168)

plasticity <- left_join(plasticity, SP_DEU_72v0)
plasticity <- left_join(plasticity, SP_DEU_168v0)
plasticity <- left_join(plasticity, SP_DEU_168v72)

# Adding inter-population DEU. Hours 0, 72, and 168
plasticity <- left_join(plasticity, CvS_DEU_hour0)
plasticity <- left_join(plasticity, CvS_DEU_hour72)
plasticity <- left_join(plasticity, CvS_DEU_hour168)


# Adding a column called plasticity, with 1 if a transcript shows DGE or DEU under any inter or intra-population comparison
plasticity <- plasticity %>% 
       mutate(plasticity = case_when(
         CvS_FDR_168 != 'NA' | 
         CV_FDR_72v0 != 'NA' |
         CV_FDR_168v0 != 'NA' |
         CV_FDR_168v72 != 'NA' |
         SP_FDR_72v0 != 'NA' |
         SP_FDR_168v0 != 'NA' |
         SP_FDR_168v72 != 'NA' |
         DEU_FDR_CV_72v0 != 'NA' |
         DEU_FDR_CV_168v0 != 'NA' |
         DEU_FDR_CV_72v168 != 'NA' |
         DEU_FDR_SP_72v0 != 'NA' |
         DEU_FDR_SP_168v0 != 'NA' |
         DEU_FDR_SP_168v72 != 'NA' |
           DEU_FDR_0 != 'NA' |
           DEU_FDR_72 != 'NA' |
           DEU_FDR_168 != 'NA'
         ~ 1))
     
# Looking at transcripts showing selection, along PC1 with PCAdapt and significant transcripts from Bayescan
selection <- full_join(clusters_w_annotations, pcadapt_snps_PC1)
selection <- full_join(selection, bayescan_output, by = c("pcadapt_SNP_ID" = "bayescan_SNP_ID"))

# Adding a column called selection, with a 1 if a transcript shows selection in PCAdapt and Bayescan
selection <- selection %>% 
  mutate(selection = case_when(
    PC_adapt_qvalue != 'NA' & bayescan_qval != 'NA'
    ~ 1
  ))

# Put the data together for a chi^2 test
abridged_selection <- select(selection, c(corset_cluster = corset_cluster.x, selection))
abridged_plasticity <- select(plasticity, c(corset_cluster, plasticity))

# Checking how many unique transcripts show signatures of selection in both Bayescan and PCadapt
filtered_selection <- dplyr::filter(abridged_selection, !is.na(selection))
filtered_selection <- dplyr::distinct(filtered_selection)

# Create a new combined dataframe. Using the table with unique transcripts showing selection so that transcripts with multiple SNPs are not weighted any differently, or duplicated for chi^2 tests.
combined_data <- left_join(clusters_w_annotations, filtered_selection)
combined_data <- left_join(combined_data, abridged_plasticity)

# Deleting rows with neither plasticity or selection
# The line is commented out because the number of rows with neither selection or plasticity are important for the chi^2
# combined_data <- combined_data[!with(combined_data,is.na(selection)& is.na(plasticity)),]


# How many transcripts show both selection & plasticity
nrow(subset(combined_data, !is.na(selection) & !is.na(plasticity)))
# 8

# How many transcripts show selection & no plasticity
nrow(subset(combined_data, !is.na(selection) & is.na(plasticity)))
# 67

# How many transcripts show no selection & plasticity
nrow(subset(combined_data, is.na(selection) & !is.na(plasticity)))
# 4880

# How many transcripts show no selection & no plasticity
nrow(subset(combined_data, is.na(selection) & is.na(plasticity)))
# 244021

###################################################################################################
# Create a data frame with the maximum absolute log-fold change on one axis, and -log10 FDR value on the other axis for signatures of selection. 
# First, between Bayescan and PCadapt, taking the higher of the two q-values for the SNPs that overlap. Taking the higher value is a conservative approach here, as the lower q-value is 'more significant'. 
selection <- selection %>% mutate(max_sel_qval = pmax(PC_adapt_qvalue, bayescan_qval))
selection <- selection %>% 
  mutate(log10_FDR = -log10(max_sel_qval))

# Pulling the maximum q value and minimum -log10 q value as conservative estimates, by transcript. 
summarized_selection <- selection %>% 
  group_by(corset_cluster.x) %>% 
  summarize(Sel.Max.q.value = max(max_sel_qval, na.rm = TRUE),
            Sel.Min.log10_FDR = min(log10_FDR, na.rm = TRUE)) %>% 
  mutate(Sel.Max.q.value = na_if(Sel.Max.q.value, "-Inf")) %>% 
  mutate(Sel.Min.log10_FDR = na_if(Sel.Min.log10_FDR, "Inf"))
summarized_selection <- rename(summarized_selection, corset_cluster = corset_cluster.x)

colMax <- function(data) sapply(data, max, na.rm = TRUE)
colMin <- function(data) sapply(data, min, na.rm = TRUE)


LFC_only <- as_tibble(t(cbind(abs(plasticity$CvS_LFC_168), abs(plasticity$CV_LFC_72v0), abs(plasticity$CV_LFC_168v0), abs(plasticity$CV_LFC_168v72), abs(plasticity$SP_LFC_72v0), abs(plasticity$SP_LFC_168v0), abs(plasticity$CV_LFC_168v72))))

summarized_plasticity <- colMax(LFC_only)
plasticity$max_LFC <- summarized_plasticity
plasticity$max_LFC <- na_if(plasticity$max_LFC, "-Inf")

#######################################################################################################################################
# To combine different kinds of plasticity (DGE, DEU, GEV), we'll try the -log10 qvalue of each kind of plasticity, then get the maximum -log10 qvalue for each transcript to compare to the -log10 qvalues of signatures of selection
# First, taking a look at only qvalues and differential gene expression
DGE_FDR_only <- as_tibble(t(cbind(-log10(plasticity$CvS_FDR_168), -log10(plasticity$CV_FDR_72v0), -log10(plasticity$CV_FDR_168v0), -log10(plasticity$CV_FDR_168v72), -log10(plasticity$SP_FDR_72v0), -log10(plasticity$SP_FDR_168v0), -log10(plasticity$CV_FDR_168v72))))
summarized_DGE_FDR <- colMax(DGE_FDR_only)
plasticity$DGE_max_logFDR <- summarized_DGE_FDR
plasticity$DGE_max_logFDR <- na_if(plasticity$DGE_max_logFDR, "-Inf")

# Checking how well LFC for DGE and -log10 Q values correlate 
dge_LFC_qval_lm <- lm(DGE_max_logFDR ~ max_LFC, data = plasticity)
summary(dge_LFC_qval_lm)
plot(plasticity$DGE_max_logFDR, plasticity$max_LFC)


FDR_only <- as.data.frame(cbind(-log10(plasticity$CvS_FDR_168), -log10(plasticity$CV_FDR_72v0), -log10(plasticity$CV_FDR_168v0), -log10(plasticity$CV_FDR_168v72), -log10(plasticity$SP_FDR_72v0), -log10(plasticity$SP_FDR_168v0), -log10(plasticity$CV_FDR_168v72), -log10(plasticity$DEU_FDR_CV_72v0), -log10(plasticity$DEU_FDR_CV_168v0), -log10(plasticity$DEU_FDR_CV_72v168), -log10(plasticity$DEU_FDR_SP_72v0), -log10(plasticity$DEU_FDR_SP_168v0), -log10(plasticity$DEU_FDR_SP_168v72), -log10(plasticity$DEU_FDR_0), -log10(plasticity$DEU_FDR_72), -log10(plasticity$DEU_FDR_168)))
FDR_only <- FDR_only %>% mutate(min_FDR = pmin(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, na.rm = TRUE))
plasticity$Plasticity.min_logFDR <- FDR_only$min_FDR


# Join the combined plasticity information with selection information
max_plasticity_selection <- full_join(plasticity, summarized_selection)

# Doing a simple plot like this misses the many transcripts under selection that do not show plasticity, and vice-versa
plot(max_plasticity_selection$Sel.Min.log10_FDR, max_plasticity_selection$Plasticity.min_logFDR)

max_plasticity_selection <- max_plasticity_selection %>%
  mutate(Sel.Min.log10_FDR = ifelse(!is.na(Plasticity.min_logFDR) & is.na(Sel.Min.log10_FDR), 0, Sel.Min.log10_FDR)) %>%
  mutate(Plasticity.min_logFDR = ifelse(!is.na(Sel.Min.log10_FDR) & is.na(Plasticity.min_logFDR), 0, Plasticity.min_logFDR))

# Doing the same plot again shows how the 0's have been added to each axis
plot(max_plasticity_selection$Sel.Min.log10_FDR, max_plasticity_selection$Plasticity.min_logFDR)
plot(max_plasticity_selection$Sel.Min.log10_FDR, max_plasticity_selection$max_LFC)

whichMax <- function(data) sapply(data, which.max)

max_plasticity_selection <- max_plasticity_selection %>% 
  mutate(plasticity_category = case_when(
    ((is.na(CvS_FDR_168) & is.na(CV_FDR_72v0) & is.na(CV_FDR_168v0) & is.na(CV_FDR_168v72) & is.na(SP_FDR_72v0) & is.na(SP_FDR_168v0) & is.na(SP_FDR_168v72)) & (is.na(DEU_FDR_CV_72v0) & is.na(DEU_FDR_CV_168v0) & is.na(DEU_FDR_CV_72v168) & is.na(DEU_FDR_SP_72v0) & is.na(DEU_FDR_SP_168v0) & is.na(DEU_FDR_SP_168v72) & is.na(DEU_FDR_0) & is.na(DEU_FDR_72) & is.na(DEU_FDR_168))) == TRUE ~ "No_Plast",
    ((!is.na(CvS_FDR_168) | !is.na(CV_FDR_72v0) | !is.na(CV_FDR_168v0) | !is.na(CV_FDR_168v72) | !is.na(SP_FDR_72v0) | !is.na(SP_FDR_168v0) | !is.na(SP_FDR_168v72)) & (!is.na(DEU_FDR_CV_72v0) | !is.na(DEU_FDR_CV_168v0) | !is.na(DEU_FDR_CV_72v168) | !is.na(DEU_FDR_SP_72v0) | !is.na(DEU_FDR_SP_168v0) | !is.na(DEU_FDR_SP_168v72) | !is.na(DEU_FDR_0) | !is.na(DEU_FDR_72) | !is.na(DEU_FDR_168))) == TRUE ~ "DGE+DEU",
    ((!is.na(CvS_FDR_168) | !is.na(CV_FDR_72v0) | !is.na(CV_FDR_168v0) | !is.na(CV_FDR_168v72) | !is.na(SP_FDR_72v0) | !is.na(SP_FDR_168v0) | !is.na(SP_FDR_168v72)) & (is.na(DEU_FDR_CV_72v0) & is.na(DEU_FDR_CV_168v0) & is.na(DEU_FDR_CV_72v168) & is.na(DEU_FDR_SP_72v0) & is.na(DEU_FDR_SP_168v0) & is.na(DEU_FDR_SP_168v72) & is.na(DEU_FDR_0) & is.na(DEU_FDR_72) & is.na(DEU_FDR_168))) == TRUE ~ "DGE",
    ((is.na(CvS_FDR_168) & is.na(CV_FDR_72v0) & is.na(CV_FDR_168v0) & is.na(CV_FDR_168v72) & is.na(SP_FDR_72v0) & is.na(SP_FDR_168v0) & is.na(SP_FDR_168v72)) & (!is.na(DEU_FDR_CV_72v0) | !is.na(DEU_FDR_CV_168v0) | !is.na(DEU_FDR_CV_72v168) | !is.na(DEU_FDR_SP_72v0) | !is.na(DEU_FDR_SP_168v0) | !is.na(DEU_FDR_SP_168v72) | !is.na(DEU_FDR_0) | !is.na(DEU_FDR_72) | !is.na(DEU_FDR_168))) == TRUE ~ "DEU"
  ))

# Plotting selection against plasticity
selection_plasticity_plot <- ggplot(data = max_plasticity_selection, mapping = aes(x = Sel.Min.log10_FDR, y = Plasticity.min_logFDR, colour = plasticity_category, group = plasticity_category)) +
  geom_point(alpha = 0.5) +
  scale_fill_manual(values = c("DGE+DEU" = "#fdae61", "DGE" = "#b2182b", "DEU" = "#4393c3", "No_Plast" = "black"), breaks = c("DGE+DEU", "DGE", "DEU", "No_Plast"), labels = c("Both DGE & DEU", "Differential Gene Expression", "Differential Exon Usage", "No Plasticity")) +
  scale_color_manual(values = c("DGE+DEU" = "#fdae61", "DGE" = "#b2182b", "DEU" = "#4393c3", "No_Plast" = "black"), breaks = c("DGE+DEU", "DGE", "DEU", "No_Plast"), labels = c("Both DGE & DEU", "Differential Gene Expression", "Differential Exon Usage", "No Plasticity")) +
  labs(x = bquote('Minimum -log'[10]~italic(q)*'-value of outlier SNP'), y = bquote('Minimum -log'[10]~italic(q)*'-value of differential gene expression or differential exon usage'), color = "Gene Category") +
  theme_bw() +
  geom_point(mapping = aes(Sel.Min.log10_FDR, y = Plasticity.min_logFDR), data = dplyr::filter(max_plasticity_selection, plasticity_category == "DGE"), alpha = 0.05) +
  geom_point(mapping = aes(Sel.Min.log10_FDR, y = Plasticity.min_logFDR), data = dplyr::filter(max_plasticity_selection, plasticity_category == "DGE+DEU"), alpha = 0.05) +
  theme(text = element_text(size=15, family="serif"), legend.justification = c("right", "top"), legend.position = c(0.99, 1.0), legend.background = element_blank())
selection_plasticity_plot
ggsave("selection_v_plasticity.pdf", selection_plasticity_plot, dpi = 2000)

ggplot(data = max_plasticity_selection, mapping = aes(x = Sel.Min.log10_FDR, y = Plasticity.min_logFDR, colour = plasticity_category, group = plasticity_category)) +
  geom_point(alpha = 0.5) +
  scale_fill_manual(values = c("DGE+DEU" = "#fdae61", "DGE" = "#b2182b", "DEU" = "#4393c3", "No_Plast" = "black"), breaks = c("DGE+DEU", "DGE", "DEU", "No_Plast"), labels = c("Both DGE & DEU", "Differential Gene Expression", "Differential Exon Usage", "No Plasticity")) +
  scale_color_manual(values = c("DGE+DEU" = "#fdae61", "DGE" = "#b2182b", "DEU" = "#4393c3", "No_Plast" = "black"), breaks = c("DGE+DEU", "DGE", "DEU", "No_Plast"), labels = c("Both DGE & DEU", "Differential Gene Expression", "Differential Exon Usage", "No Plasticity")) +
  labs(x = bquote('Minimum -log'[10]~italic(q)*'-value of outlier SNP'), y = bquote('Minimum -log'[10]~italic(q)*'-value of differential gene expression or differential exon usage'), color = "Gene Category") +
  theme_bw() +
  geom_point(mapping = aes(Sel.Min.log10_FDR, y = Plasticity.min_logFDR), data = dplyr::filter(max_plasticity_selection, plasticity_category == "DGE"), alpha = 0.05) +
  geom_point(mapping = aes(Sel.Min.log10_FDR, y = Plasticity.min_logFDR), data = dplyr::filter(max_plasticity_selection, plasticity_category == "DGE+DEU"), alpha = 0.05) +
  theme(text = element_text(size=15, family="serif"), legend.justification = c("right", "top"), legend.position = c(0.99, 1.0), legend.background = element_blank())

# Pulling out the lists of genes under different conditions
# First, some genes that showed both selection and plasticity
sel_plasticity <- dplyr::filter(max_plasticity_selection, Sel.Min.log10_FDR > 0 & Plasticity.min_logFDR > 0)
deu <- dplyr::filter(max_plasticity_selection, plasticity_category == "DEU")

# For gene set enrichment analysis, exporting the combined data to a new R script
save(max_plasticity_selection, file = "plasticity_selection.RData")
