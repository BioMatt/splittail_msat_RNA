# This script takes the combined table looking at plasticity and selection in splittail, and runs gene set enrichment analysis with EnrichR to take a look at summary terms

library(tidyverse)
library(enrichR)
listEnrichrDbs()
split_go <- function(x) {
  # Taking apart the GO descriptions and GO ID terms
  x <- x %>%
    select(-starts_with("Old")) %>%
    separate(Term, c("GO_term", "GO_ID"), sep = "GO")
  
  # Changing all the floating :#### GO terms to be GO:####
  x <- mutate_if(x, is.character, str_replace_all, pattern = ":", replacement = "GO:")
  
  # Removing the last parentheses in the GO ID 
  x$GO_ID <- x$GO_ID %>% 
    str_replace("\\)", "")
  
  # Removing the last parentheses in the GO term to have clean looking cells 
  x$GO_term <- x$GO_term %>% 
    str_replace("\\($", "")
  return(x)
}
load("plasticity_selection.RData")

#Going in order of the manuscript, looking at between-population DEU at hour 0
deu_hour0 <- dplyr::filter(max_plasticity_selection, !is.na(DEU_FDR_0))
# Hour 72 between populations
deu_hour72 <- dplyr::filter(max_plasticity_selection, !is.na(DEU_FDR_72))
# Hour 168 between populations
deu_hour168 <- dplyr::filter(max_plasticity_selection, !is.na(DEU_FDR_168))
# 72 genes have annotations, using EnrichR to summarize the data. No GO term is significant.
deu_hour168_enriched <- enrichr(c(deu_hour168$uniprot_gene), "GO_Biological_Process_2018")
deu_hour168_bioprocess <- split_go(deu_hour168_enriched$GO_Biological_Process_2018)
# Writing out the DEU hour 168 table for inclusion in the supplementary materials
write_delim(deu_hour168, "DEU_hour168.txt", delim = "\t")

# With between-population DEU done, now intra-population DEU. First, Central Valley fish.
deu_CV_72v0 <- dplyr::filter(max_plasticity_selection, !is.na(DEU_FDR_CV_72v0))
deu_CV_168v72 <- dplyr::filter(max_plasticity_selection, !is.na(DEU_FDR_CV_72v168))

# Intra-population DEU in the San Pablo fish
deu_SP_72v0 <- dplyr::filter(max_plasticity_selection, !is.na(DEU_FDR_SP_72v0))
deu_SP_168v0 <- dplyr::filter(max_plasticity_selection, !is.na(DEU_FDR_SP_168v0))
# Filtering for genes with annotations, so removing ones named 0 or NA
deu_SP_168v0 <- dplyr::filter(deu_SP_168v0, !is.na(uniprot_gene) & uniprot_gene != "0")
length(unique(deu_SP_168v0$uniprot_gene))

# Looking at GO terms enriched in this list of genes
deu_SP_168v0_enriched <- enrichr(c(deu_SP_168v0$uniprot_gene), "GO_Biological_Process_2018")
# Use a custom script to take the table from enrichR and split GO terms from descriptions, along with reformatting it into a table
deu_SP_168v0_bioprocess <- split_go(deu_SP_168v0_enriched$GO_Biological_Process_2018)
write_delim(deu_SP_168v0, "DEU_SP_hour168v0.txt", delim = "\t")


# Looking at the 168 v 72 hour comparison of DEU for San Pablo fish, with 2,697 significant Corset-clustered reads (genes)
deu_SP_168v72 <- dplyr::filter(max_plasticity_selection, !is.na(DEU_FDR_SP_168v72))
# Filtering for genes with annotations, so removing ones named 0 or NA
deu_SP_168v72 <- dplyr::filter(deu_SP_168v72, !is.na(uniprot_gene) & uniprot_gene != "0")
#write_delim(deu_SP_168v72, "DEU_SP_hour168v72.txt", delim = "\t")

# Looking at GO terms enriched in this list of genes
deu_SP_168v72_enriched <- enrichr(c(deu_SP_168v72$uniprot_gene), "GO_Biological_Process_2018")
# Use a custom script to take the table from enrichR and split GO terms from descriptions, along with reformatting it into a table
deu_SP_168v72_bioprocess <- split_go(deu_SP_168v72_enriched$GO_Biological_Process_2018)
