if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR", version = "3.8")


#install.packages("gamlss")

library(tidyverse)
library(edgeR)
library(gamlss)

# Looking at variance in walleye expression data, following: a combination of the https://github.com/Oshlack/Corset/wiki/Example, edgeR user's guide, and https://github.com/Vityay/ExpVarQuant/blob/master/ExpVarQuant.R

# Reading the corset gene counts
corset_counts <- read.delim("splittail-counts.txt", row.names = 1)


# This section involves importing Corset output data, adding in some metadata about site and year, then creating the EdgeR table.

# Read in a table with my sample information and other metadata about them
samples <- read.table("splittail_metadata.txt", header = TRUE, sep = "\t")

# Create the EdgeR dge object, at the 0, 3, and 7 day time points of exposure 
Counts_edgeR_0 <- DGEList(counts = corset_counts[,c(1:4,17:20)], group = samples$Population[c(1:4,17:20)])
Counts_edgeR_3 <- DGEList(counts = corset_counts[,c(5:10,21:26)], group = samples$Population[c(5:10,21:26)])
Counts_edgeR_7 <- DGEList(counts = corset_counts[,c(11:16,27:32)], group = samples$Population[c(11:16,27:32)])


# Calculate library size and normalization factor (see edgeR manual). Starting with the time 0 data
Counts_edgeR_0 <- calcNormFactors(Counts_edgeR_0, method="upperquartile")
# Calculate counts per million (CPM).
Counts_edgeR_0$CPM <- cpm.DGEList(Counts_edgeR_0)
# Calculate offset variable for each library as natural log of the product of library size and normalization factor.
ofs_0 <- log(Counts_edgeR_0$samples$lib.size*Counts_edgeR_0$samples$norm.factors)
Counts_edgeR_0$samples$offset <- ofs_0

# Remove genes with 0 counts in any of the sample as analysis of zero inflated samples could introduce a significant bias.
idx_0 <- !apply(Counts_edgeR_0$CPM,1, function(x) any(x == 0) )
Counts_edgeR_0$counts <- Counts_edgeR_0$counts[idx_0,]
Counts_edgeR_0$CPM <- Counts_edgeR_0$CPM[idx_0,]

# Remove lowly expressed genes as for these genes technical variations might have a substantial impact on gene noise estimation.
idx_0 <- rowMeans(Counts_edgeR_0$CPM) > 1
Counts_edgeR_0$counts <- Counts_edgeR_0$counts[idx_0,]
Counts_edgeR_0$CPM <- Counts_edgeR_0$CPM[idx_0,]

# Remove idx and ofs variables from the workspace.
rm(idx_0, ofs_0)

# Estimate exposure time effects on mean and overdispersion parameters with GAMLSS for each gene.
gene_i_0 <- seq_along(Counts_edgeR_0$counts[,1])

# To try algorithm for just some genes, change gene_i variable. For example, set gene_i <- c(1:100) to estimate GAMLSS models for the first hundred genes.
gamlss_NB_0 <- lapply(gene_i_0, function(i) {
  
  # For each gene (i) a table containing: x - a factor of interest (population); y - RNA-seq. counts and offset (ofs) is created.
  dat <- data.frame(
    x = Counts_edgeR_0$samples$group,
    y = Counts_edgeR_0$counts[i,],
    ofs = Counts_edgeR_0$samples$offset
  )
  # x is releveled to use Central Valley splittail as a reference.
  dat$x <- relevel(dat$x, ref = c("C"))
  
  # Fit negative binomial (family = NBI()) GAMLSS model, which accounts for group effects on mean and overdispersion (non-Poisson noise).
  # fo = y~0+x+offset(ofs) specifies model for mean and sigma.fo=~0+x for overdispersion, offset - offset(ofs) normalize counts to library size. sigma.start = 0.1 provides starting value for overdispersion estimation (default is 1). n.cyc – number of fitting algorithm cycles, see help(gamlss).
  # In some cases, fitting of NB model may fail and tryCatch(..., warning= function(w) NULL, error= function(e) NULL) will return NULL as a result.
  m0 <- tryCatch(
    gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 0+x, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit reduced model by omitting group factor from the estimation of overdispersion: sigma.fo = ~ 1. In essence, this model corresponds to the GLM model implemented in edgeR.
  m1 <- tryCatch(
    gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 1, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit reduced model by omitting group factor from the estimation of mean: fo = y ~ offset(ofs).
  m2 <- tryCatch(
    gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 0+x, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit null model.
  m3 <- tryCatch(
    gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 1, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Create data frame res to store the results.
  res <- data.frame(
    cpm.CV = NA,
    cpm.SP=NA,
    LR.cpm = NA,
    p_gamlss.cpm = NA,
    p_glm.cpm = NA,
    CV.CV = NA,
    CV.SP = NA,
    LR.cv = NA,
    p.cv = NA
  )
  
  # Because fitting of the NB model may fail for some genes, check whether all models were fitted successfully. 
  if(!any(sapply(list(m0,m1,m2,m3), is.null))) 
  {
    # Write GAMLSS estimations of gene’s mean (CPM) counts from the m0 model.
    res$cpm.CV = exp(m0$mu.coefficients+log(1e06))[[1]]
    res$cpm.SP = exp(m0$mu.coefficients+log(1e06))[[2]]
    
    # Calculate log2 ratio for changes in CPMs between San Pablo and Central Valley splittail.
    res$LR.cpm = log2(exp(m0$mu.coefficients+log(1e06))[[2]]/exp(m0$mu.coefficients+log(1e06))[[1]])
    
    # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on gene’s mean (CPM) counts.
    # p_gamlss.cpm – p value of LR test statistic: D_μ=-2loga[L(μ_0,α_j  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m2 models.
    res$p_gamlss.cpm = pchisq(2*(logLik(m0)-logLik(m2)), df=m0$df.fit-m2$df.fit, lower=F)[[1]]
    
    # GLM log-likelihood ratio (LR) test for a significance of age effect on gene’s mean (CPM) counts.
    # p_glm.cpm – p value of LR test statistic: D_(μ_GLM )=-2loga[L(μ_0,α_0  ┤|  X_ij)/L(μ_j,α_0  ┤|  X_ij), comparing m1 and m3 models.
    res$p_glm.cpm = pchisq(2*(logLik(m1)-logLik(m3)), df=m1$df.fit-m3$df.fit, lower=F)[[1]]
    
    # Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: cv(μ)=√α.
    res$CV.CV  = sqrt(exp(m0$sigma.coefficients)[[1]])
    res$CV.SP = sqrt(exp(m0$sigma.coefficients)[[2]])
    
    # Calculate log2 ratio for changes in cv(μ) between San Pablo and Central Valley splittail.
    res$LR.cv = log2(sqrt(exp(m0$sigma.coefficients)[[2]])/sqrt(exp(m0$sigma.coefficients)[[1]]))
    
    # GAMLSS log-likelihood ratio (LR) test for a significance of an group effect on non-Poisson noise.
    # p.cv – p value of LR test statistic: D_α=-2loga[L(μ_j,α_0  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
    res$p.cv = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
  }
  res
})

# Transform list gamlss_NB containing GAMLSS estimations for each gene to data frame
gamlss_NB_0 <- do.call(rbind, gamlss_NB_0)
rownames(gamlss_NB_0) <- rownames(Counts_edgeR_0$counts)[gene_i_0]

# Because GAMLSS fitting might fail for some genes or yield inflated estimates of overdispersion, the results have to be cleaned.
# First, remove genes, for which GAMLSS model has failed.
gamlss_NB_clean_0 <- na.omit(gamlss_NB_0)

# Second, remove genes, for which estimates of cv(μ)=√α were either inflated > 3 or close to Poisson < 10-3.
idx_0 <- gamlss_NB_clean_0$CV.CV > 3 | gamlss_NB_clean_0$CV.CV < 1e-03 | gamlss_NB_clean_0$CV.SP > 3 | gamlss_NB_clean_0$CV.SP < 1e-03
gamlss_NB_clean_0 <- gamlss_NB_clean_0[!idx_0,]

# Finally, calculate false discovery rates to account for multiple hypothesis testing.
gamlss_NB_clean_0$padj_gamlss.cpm <- p.adjust(gamlss_NB_clean_0$p_gamlss.cpm, "fdr")
gamlss_NB_clean_0$padj_glm.cpm <- p.adjust(gamlss_NB_clean_0$p_glm.cpm, "fdr")
gamlss_NB_clean_0$padj.cv <- p.adjust(gamlss_NB_clean_0$p.cv, "fdr")

# Save the results, which later can be loaded to R and used for further analysis
rm(idx_0, gene_i_0, gamlss_NB_0)

################################################################################################################################################################################################

# Calculate library size and normalization factor (see edgeR manual). Here, using the 3 day exposure data
Counts_edgeR_3 <- calcNormFactors(Counts_edgeR_3, method="upperquartile")
# Calculate counts per million (CPM).
Counts_edgeR_3$CPM <- cpm.DGEList(Counts_edgeR_3)
# Calculate offset variable for each library as natural log of the product of library size and normalization factor.
ofs_3 <- log(Counts_edgeR_3$samples$lib.size*Counts_edgeR_3$samples$norm.factors)
Counts_edgeR_3$samples$offset <- ofs_3

# Remove genes with 0 counts in any of the sample as analysis of zero inflated samples could introduce a significant bias.
idx_3 <- !apply(Counts_edgeR_3$CPM,1, function(x) any(x == 0) )
Counts_edgeR_3$counts <- Counts_edgeR_3$counts[idx_3,]
Counts_edgeR_3$CPM <- Counts_edgeR_3$CPM[idx_3,]

# Remove lowly expressed genes as for these genes technical variations might have a substantial impact on gene noise estimation.
idx_3 <- rowMeans(Counts_edgeR_3$CPM) > 1
Counts_edgeR_3$counts <- Counts_edgeR_3$counts[idx_3,]
Counts_edgeR_3$CPM <- Counts_edgeR_3$CPM[idx_3,]

# Remove idx and ofs variables from the workspace.
rm(idx_3, ofs_3)

# Estimate age effects on mean and overdispersion parameters with GAMLSS for each gene.
gene_i_3 <- seq_along(Counts_edgeR_3$counts[,1])

# To try algorithm for just some genes, change gene_i variable. For example, set gene_i <- c(1:100) to estimate GAMLSS models for the first hundred genes.
gamlss_NB_3 <- lapply(gene_i_3, function(i) {
  
  # For each gene (i) a table containing: x - a factor of interest (population); y - RNA-seq. counts and offset (ofs) is created.
  dat <- data.frame(
    x = Counts_edgeR_3$samples$group,
    y = Counts_edgeR_3$counts[i,],
    ofs = Counts_edgeR_3$samples$offset
  )
  # x is releveled to use Central Valley splittail as a reference.
  dat$x <- relevel(dat$x, ref = c("C"))
  
  # Fit negative binomial (family = NBI()) GAMLSS model, which accounts for age effects on mean and overdispersion (non-Poisson noise).
  # fo = y~0+x+offset(ofs) specifies model for mean and sigma.fo=~0+x for overdispersion, offset - offset(ofs) normalize counts to library size. sigma.start = 0.1 provides starting value for overdispersion estimation (default is 1). n.cyc – number of fitting algorithm cycles, see help(gamlss).
  # In some cases, fitting of NB model may fail and tryCatch(..., warning= function(w) NULL, error= function(e) NULL) will return NULL as a result.
  m0 <- tryCatch(
    gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 0+x, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit reduced model by omitting age factor from the estimation of overdispersion: sigma.fo = ~ 1. In essence, this model corresponds to the GLM model implemented in edgeR.
  m1 <- tryCatch(
    gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 1, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit reduced model by omitting age factor from the estimation of mean: fo = y ~ offset(ofs).
  m2 <- tryCatch(
    gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 0+x, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit null model.
  m3 <- tryCatch(
    gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 1, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Create data frame res to store the results.
  res <- data.frame(
    cpm.CV = NA,
    cpm.SP=NA,
    LR.cpm = NA,
    p_gamlss.cpm = NA,
    p_glm.cpm = NA,
    CV.CV = NA,
    CV.SP = NA,
    LR.cv = NA,
    p.cv = NA
  )
  
  # Because fitting of the NB model may fail for some genes, check whether all models were fitted successfully. 
  if(!any(sapply(list(m0,m1,m2,m3), is.null))) 
  {
    # Write GAMLSS estimations of gene’s mean (CPM) counts from the m0 model.
    res$cpm.CV = exp(m0$mu.coefficients+log(1e06))[[1]]
    res$cpm.SP = exp(m0$mu.coefficients+log(1e06))[[2]]
    
    # Calculate log2 ratio for changes in CPMs between San Pablo and Central Valley splittail.
    res$LR.cpm = log2(exp(m0$mu.coefficients+log(1e06))[[2]]/exp(m0$mu.coefficients+log(1e06))[[1]])
    
    # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on gene’s mean (CPM) counts.
    # p_gamlss.cpm – p value of LR test statistic: D_μ=-2loga[L(μ_3,α_j  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m2 models.
    res$p_gamlss.cpm = pchisq(2*(logLik(m0)-logLik(m2)), df=m0$df.fit-m2$df.fit, lower=F)[[1]]
    
    # GLM log-likelihood ratio (LR) test for a significance of age effect on gene’s mean (CPM) counts.
    # p_glm.cpm – p value of LR test statistic: D_(μ_GLM )=-2loga[L(μ_3,α_3  ┤|  X_ij)/L(μ_j,α_3  ┤|  X_ij), comparing m1 and m3 models.
    res$p_glm.cpm = pchisq(2*(logLik(m1)-logLik(m3)), df=m1$df.fit-m3$df.fit, lower=F)[[1]]
    
    # Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: cv(μ)=√α.
    res$CV.CV  = sqrt(exp(m0$sigma.coefficients)[[1]])
    res$CV.SP = sqrt(exp(m0$sigma.coefficients)[[2]])
    
    # Calculate log2 ratio for changes in cv(μ) between San Pablo and Central Valley splittail.
    res$LR.cv = log2(sqrt(exp(m0$sigma.coefficients)[[2]])/sqrt(exp(m0$sigma.coefficients)[[1]]))
    
    # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on non-Poisson noise.
    # p.cv – p value of LR test statistic: D_α=-2loga[L(μ_j,α_3  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
    res$p.cv = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
  }
  res
})

# Transform list gamlss_NB containing GAMLSS estimations for each gene to data frame
gamlss_NB_3 <- do.call(rbind, gamlss_NB_3)
rownames(gamlss_NB_3) <- rownames(Counts_edgeR_3$counts)[gene_i_3]

# Because GAMLSS fitting might fail for some genes or yield inflated estimates of overdispersion, the results have to be cleaned.
# First, remove genes, for which GAMLSS model has failed.
gamlss_NB_clean_3 <- na.omit(gamlss_NB_3)

# Second, remove genes, for which estimates of cv(μ)=√α were either inflated > 3 or close to Poisson < 10-3.
idx_3 <- gamlss_NB_clean_3$CV.CV > 3 | gamlss_NB_clean_3$CV.CV < 1e-03 | gamlss_NB_clean_3$CV.SP > 3 | gamlss_NB_clean_3$CV.SP < 1e-03
gamlss_NB_clean_3 <- gamlss_NB_clean_3[!idx_3,]

# Finally, calculate false discovery rates to account for multiple hypothesis testing.
gamlss_NB_clean_3$padj_gamlss.cpm <- p.adjust(gamlss_NB_clean_3$p_gamlss.cpm, "fdr")
gamlss_NB_clean_3$padj_glm.cpm <- p.adjust(gamlss_NB_clean_3$p_glm.cpm, "fdr")
gamlss_NB_clean_3$padj.cv <- p.adjust(gamlss_NB_clean_3$p.cv, "fdr")

# Save the results, which later can be loaded to R and used for further analysis
rm(idx_3, gene_i_3, gamlss_NB_3)

################################################################################################################################################################################################

# Calculate library size and normalization factor (see edgeR manual). Here, using the 7 day exposure data
Counts_edgeR_7 <- calcNormFactors(Counts_edgeR_7, method="upperquartile")
# Calculate counts per million (CPM).
Counts_edgeR_7$CPM <- cpm.DGEList(Counts_edgeR_7)
# Calculate offset variable for each library as natural log of the product of library size and normalization factor.
ofs_7 <- log(Counts_edgeR_7$samples$lib.size*Counts_edgeR_7$samples$norm.factors)
Counts_edgeR_7$samples$offset <- ofs_7

# Remove genes with 0 counts in any of the sample as analysis of zero inflated samples could introduce a significant bias.
idx_7 <- !apply(Counts_edgeR_7$CPM,1, function(x) any(x == 0) )
Counts_edgeR_7$counts <- Counts_edgeR_7$counts[idx_7,]
Counts_edgeR_7$CPM <- Counts_edgeR_7$CPM[idx_7,]

# Remove lowly expressed genes as for these genes technical variations might have a substantial impact on gene noise estimation.
idx_7 <- rowMeans(Counts_edgeR_7$CPM) > 1
Counts_edgeR_7$counts <- Counts_edgeR_7$counts[idx_7,]
Counts_edgeR_7$CPM <- Counts_edgeR_7$CPM[idx_7,]

# Remove idx and ofs variables from the workspace.
rm(idx_7, ofs_7)

# Estimate age effects on mean and overdispersion parameters with GAMLSS for each gene.
gene_i_7 <- seq_along(Counts_edgeR_7$counts[,1])

# To try algorithm for just some genes, change gene_i variable. For example, set gene_i <- c(1:100) to estimate GAMLSS models for the first hundred genes.
gamlss_NB_7 <- lapply(gene_i_7, function(i) {
  
  # For each gene (i) a table containing: x - a factor of interest (population); y - RNA-seq. counts and offset (ofs) is created.
  dat <- data.frame(
    x = Counts_edgeR_7$samples$group,
    y = Counts_edgeR_7$counts[i,],
    ofs = Counts_edgeR_7$samples$offset
  )
  # x is releveled to use Central Valley splittail as a reference.
  dat$x <- relevel(dat$x, ref = c("C"))
  
  # Fit negative binomial (family = NBI()) GAMLSS model, which accounts for age effects on mean and overdispersion (non-Poisson noise).
  # fo = y~0+x+offset(ofs) specifies model for mean and sigma.fo=~0+x for overdispersion, offset - offset(ofs) normalize counts to library size. sigma.start = 0.1 provides starting value for overdispersion estimation (default is 1). n.cyc – number of fitting algorithm cycles, see help(gamlss).
  # In some cases, fitting of NB model may fail and tryCatch(..., warning= function(w) NULL, error= function(e) NULL) will return NULL as a result.
  m0 <- tryCatch(
    gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 0+x, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit reduced model by omitting age factor from the estimation of overdispersion: sigma.fo = ~ 1. In essence, this model corresponds to the GLM model implemented in edgeR.
  m1 <- tryCatch(
    gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 1, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit reduced model by omitting age factor from the estimation of mean: fo = y ~ offset(ofs).
  m2 <- tryCatch(
    gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 0+x, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Fit null model.
  m3 <- tryCatch(
    gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 1, data=dat,
           family = NBI(), sigma.start = 0.1, n.cyc = 100),
    warning= function(w) NULL, error= function(e) NULL
  )
  
  # Create data frame res to store the results.
  res <- data.frame(
    cpm.CV = NA,
    cpm.SP=NA,
    LR.cpm = NA,
    p_gamlss.cpm = NA,
    p_glm.cpm = NA,
    CV.CV = NA,
    CV.SP = NA,
    LR.cv = NA,
    p.cv = NA
  )
  
  # Because fitting of the NB model may fail for some genes, check whether all models were fitted successfully. 
  if(!any(sapply(list(m0,m1,m2,m3), is.null))) 
  {
    # Write GAMLSS estimations of gene’s mean (CPM) counts from the m0 model.
    res$cpm.CV = exp(m0$mu.coefficients+log(1e06))[[1]]
    res$cpm.SP = exp(m0$mu.coefficients+log(1e06))[[2]]
    
    # Calculate log2 ratio for changes in CPMs between San Pablo and Central Valley splittail.
    res$LR.cpm = log2(exp(m0$mu.coefficients+log(1e06))[[2]]/exp(m0$mu.coefficients+log(1e06))[[1]])
    
    # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on gene’s mean (CPM) counts.
    # p_gamlss.cpm – p value of LR test statistic: D_μ=-2loga[L(μ_7,α_j  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m2 models.
    res$p_gamlss.cpm = pchisq(2*(logLik(m0)-logLik(m2)), df=m0$df.fit-m2$df.fit, lower=F)[[1]]
    
    # GLM log-likelihood ratio (LR) test for a significance of age effect on gene’s mean (CPM) counts.
    # p_glm.cpm – p value of LR test statistic: D_(μ_GLM )=-2loga[L(μ_7,α_7  ┤|  X_ij)/L(μ_j,α_7  ┤|  X_ij), comparing m1 and m3 models.
    res$p_glm.cpm = pchisq(2*(logLik(m1)-logLik(m3)), df=m1$df.fit-m3$df.fit, lower=F)[[1]]
    
    # Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: cv(μ)=√α.
    res$CV.CV  = sqrt(exp(m0$sigma.coefficients)[[1]])
    res$CV.SP = sqrt(exp(m0$sigma.coefficients)[[2]])
    
    # Calculate log2 ratio for changes in cv(μ) between San Pablo and Central Valley splittail.
    res$LR.cv = log2(sqrt(exp(m0$sigma.coefficients)[[2]])/sqrt(exp(m0$sigma.coefficients)[[1]]))
    
    # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on non-Poisson noise.
    # p.cv – p value of LR test statistic: D_α=-2loga[L(μ_j,α_7  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
    res$p.cv = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
  }
  res
})

# Transform list gamlss_NB containing GAMLSS estimations for each gene to data frame
gamlss_NB_7 <- do.call(rbind, gamlss_NB_7)
rownames(gamlss_NB_7) <- rownames(Counts_edgeR_7$counts)[gene_i_7]

# Because GAMLSS fitting might fail for some genes or yield inflated estimates of overdispersion, the results have to be cleaned.
# First, remove genes, for which GAMLSS model has failed.
gamlss_NB_clean_7 <- na.omit(gamlss_NB_7)

# Second, remove genes, for which estimates of cv(μ)=√α were either inflated > 3 or close to Poisson < 10-3.
idx_7 <- gamlss_NB_clean_7$CV.CV > 3 | gamlss_NB_clean_7$CV.CV < 1e-03 | gamlss_NB_clean_7$CV.SP > 3 | gamlss_NB_clean_7$CV.SP < 1e-03
gamlss_NB_clean_7 <- gamlss_NB_clean_7[!idx_7,]

# Finally, calculate false discovery rates to account for multiple hypothesis testing.
gamlss_NB_clean_7$padj_gamlss.cpm <- p.adjust(gamlss_NB_clean_7$p_gamlss.cpm, "fdr")
gamlss_NB_clean_7$padj_glm.cpm <- p.adjust(gamlss_NB_clean_7$p_glm.cpm, "fdr")
gamlss_NB_clean_7$padj.cv <- p.adjust(gamlss_NB_clean_7$p.cv, "fdr")

# Save the results, which later can be loaded to R and used for further analysis
rm(idx_7, gene_i_7, gamlss_NB_7)



# Writing a function to speed up the process of finding expression variation from a dataset

expression_variation <- function(data, filename){
  # Calculate library size and normalization factor (see edgeR manual).
  data <- calcNormFactors(data, method="upperquartile")
# Calculate counts per million (CPM).
  data$CPM <- cpm.DGEList(data)
  # Calculate offset variable for each library as natural log of the product of library size and normalization factor.
  ofs <- log(data$samples$lib.size*data$samples$norm.factors)
  data$samples$offset <- ofs
  
  # Remove genes with 0 counts in any of the sample as analysis of zero inflated samples could introduce a significant bias.
  idx <- !apply(data$CPM,1, function(x) any(x == 0) )
  data$counts <- data$counts[idx,]
  data$CPM <- data$CPM[idx,]
  
  # Remove lowly expressed genes as for these genes technical variations might have a substantial impact on gene noise estimation.
  idx <- rowMeans(data$CPM) > 1
  data$counts <- data$counts[idx,]
  data$CPM <- data$CPM[idx,]
  
  # Remove idx and ofs variables from the workspace.
  rm(idx, ofs)
  
  # Estimate age effects on mean and overdispersion parameters with GAMLSS for each gene.
  gene_i <- seq_along(data$counts[,1])
  
  # To try algorithm for just some genes, change gene_i variable. For example, set gene_i <- c(1:100) to estimate GAMLSS models for the first hundred genes.
  gamlss_NB <- lapply(gene_i, function(i) {
    
    # For each gene (i) a table containing: x - a factor of interest (population); y - RNA-seq. counts and offset (ofs) is created.
    dat <- data.frame(
      x = data$samples$group,
      y = data$counts[i,],
      ofs = data$samples$offset
    )
    # x is releveled to use Central Valley splittail as a reference.
    dat$x <- relevel(dat$x, ref = as.character(data$samples$group[1]))
    
    # Fit negative binomial (family = NBI()) GAMLSS model, which accounts for age effects on mean and overdispersion (non-Poisson noise).
    # fo = y~0+x+offset(ofs) specifies model for mean and sigma.fo=~0+x for overdispersion, offset - offset(ofs) normalize counts to library size. sigma.start = 0.1 provides starting value for overdispersion estimation (default is 1). n.cyc – number of fitting algorithm cycles, see help(gamlss).
    # In some cases, fitting of NB model may fail and tryCatch(..., warning= function(w) NULL, error= function(e) NULL) will return NULL as a result.
    m0 <- tryCatch(
      gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 0+x, data=dat,
             family = NBI(), sigma.start = 0.1, n.cyc = 100),
      warning= function(w) NULL, error= function(e) NULL
    )
    
    # Fit reduced model by omitting age factor from the estimation of overdispersion: sigma.fo = ~ 1. In essence, this model corresponds to the GLM model implemented in edgeR.
    m1 <- tryCatch(
      gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 1, data=dat,
             family = NBI(), sigma.start = 0.1, n.cyc = 100),
      warning= function(w) NULL, error= function(e) NULL
    )
    
    # Fit reduced model by omitting age factor from the estimation of mean: fo = y ~ offset(ofs).
    m2 <- tryCatch(
      gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 0+x, data=dat,
             family = NBI(), sigma.start = 0.1, n.cyc = 100),
      warning= function(w) NULL, error= function(e) NULL
    )
    
    # Fit null model.
    m3 <- tryCatch(
      gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 1, data=dat,
             family = NBI(), sigma.start = 0.1, n.cyc = 100),
      warning= function(w) NULL, error= function(e) NULL
    )
    
    # Create data frame res to store the results.
    res <- data.frame(
      cpm.CV = NA,
      cpm.SP=NA,
      LR.cpm = NA,
      p_gamlss.cpm = NA,
      p_glm.cpm = NA,
      CV.CV = NA,
      CV.SP = NA,
      LR.cv = NA,
      p.cv = NA
    )
    
    # Because fitting of the NB model may fail for some genes, check whether all models were fitted successfully. 
    if(!any(sapply(list(m0,m1,m2,m3), is.null))) 
    {
      # Write GAMLSS estimations of gene’s mean (CPM) counts from the m0 model.
      res$cpm.CV = exp(m0$mu.coefficients+log(1e06))[[1]]
      res$cpm.SP = exp(m0$mu.coefficients+log(1e06))[[2]]
      
      # Calculate log2 ratio for changes in CPMs between San Pablo and Central Valley splittail.
      res$LR.cpm = log2(exp(m0$mu.coefficients+log(1e06))[[2]]/exp(m0$mu.coefficients+log(1e06))[[1]])
      
      # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on gene’s mean (CPM) counts.
      # p_gamlss.cpm – p value of LR test statistic: D_μ=-2loga[L(μ_7,α_j  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m2 models.
      res$p_gamlss.cpm = pchisq(2*(logLik(m0)-logLik(m2)), df=m0$df.fit-m2$df.fit, lower=F)[[1]]
      
      # GLM log-likelihood ratio (LR) test for a significance of age effect on gene’s mean (CPM) counts.
      # p_glm.cpm – p value of LR test statistic: D_(μ_GLM )=-2loga[L(μ_7,α_7  ┤|  X_ij)/L(μ_j,α_7  ┤|  X_ij), comparing m1 and m3 models.
      res$p_glm.cpm = pchisq(2*(logLik(m1)-logLik(m3)), df=m1$df.fit-m3$df.fit, lower=F)[[1]]
      
      # Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: cv(μ)=√α.
      res$CV.CV  = sqrt(exp(m0$sigma.coefficients)[[1]])
      res$CV.SP = sqrt(exp(m0$sigma.coefficients)[[2]])
      
      # Calculate log2 ratio for changes in cv(μ) between San Pablo and Central Valley splittail.
      res$LR.cv = log2(sqrt(exp(m0$sigma.coefficients)[[2]])/sqrt(exp(m0$sigma.coefficients)[[1]]))
      
      # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on non-Poisson noise.
      # p.cv – p value of LR test statistic: D_α=-2loga[L(μ_j,α_7  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
      res$p.cv = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
    }
    res
  })
  
  # Transform list gamlss_NB containing GAMLSS estimations for each gene to data frame
  gamlss_NB <- do.call(rbind, gamlss_NB)
  rownames(gamlss_NB) <- rownames(data$counts)[gene_i]
  
  # Because GAMLSS fitting might fail for some genes or yield inflated estimates of overdispersion, the results have to be cleaned.
  # First, remove genes, for which GAMLSS model has failed.
  gamlss_NB_clean <- na.omit(gamlss_NB)
  
  # Second, remove genes, for which estimates of cv(μ)=√α were either inflated > 3 or close to Poisson < 10-3.
  idx <- gamlss_NB_clean$CV.CV > 3 | gamlss_NB_clean$CV.CV < 1e-03 | gamlss_NB_clean$CV.SP > 3 | gamlss_NB_clean$CV.SP < 1e-03
  gamlss_NB_clean <- gamlss_NB_clean[!idx,]
  
  # Finally, calculate false discovery rates to account for multiple hypothesis testing.
  gamlss_NB_clean$padj_gamlss.cpm <- p.adjust(gamlss_NB_clean$p_gamlss.cpm, "fdr")
  gamlss_NB_clean$padj_glm.cpm <- p.adjust(gamlss_NB_clean$p_glm.cpm, "fdr")
  gamlss_NB_clean$padj.cv <- p.adjust(gamlss_NB_clean$p.cv, "fdr")
  write.table(gamlss_NB_clean, filename, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
  gamlss_NB_clean <<- gamlss_NB_clean
}

expression_variation(Counts_edgeR_0, "Counts_edgeR_0_test.txt")

# Creating the 6 intra-population edgeR data sets, running the models, then saving the data to the workspace
Counts_edgeR_CV_3v0 <- DGEList(counts = corset_counts[,c(5:10,1:4)], group = samples$Exposure_Time[c(5:10,1:4)])
expression_variation(Counts_edgeR_CV_3v0, "Counts_edgeR_CV_3v0.txt")
gamlss_NB_clean_CV3v0 <- gamlss_NB_clean

Counts_edgeR_CV_7v0 <- DGEList(counts = corset_counts[,c(11:16,1:4)], group = samples$Exposure_Time[c(11:16,1:4)])
expression_variation(Counts_edgeR_CV_7v0, "Counts_edgeR_CV_7v0.txt")
gamlss_NB_clean_CV7v0 <- gamlss_NB_clean

Counts_edgeR_CV_7v3 <- DGEList(counts = corset_counts[,c(11:16,5:10)], group = samples$Exposure_Time[c(11:16,5:10)])
expression_variation(Counts_edgeR_CV_7v3, "Counts_edgeR_CV_7v3.txt")
gamlss_NB_clean_CV7v3 <- gamlss_NB_clean

Counts_edgeR_SP_3v0 <- DGEList(counts = corset_counts[,c(21:26, 17:20)], group = samples$Exposure_Time[c(21:26, 17:20)])
expression_variation(Counts_edgeR_SP_3v0, "Counts_edgeR_SP_3v0.txt")
gamlss_NB_clean_SP3v0 <- gamlss_NB_clean

Counts_edgeR_SP_7v0 <- DGEList(counts = corset_counts[,c(27:32,17:20)], group = samples$Exposure_Time[c(27:32,17:20)])
expression_variation(Counts_edgeR_SP_7v0, "Counts_edgeR_SP_7v0.txt")
gamlss_NB_clean_SP7v0 <- gamlss_NB_clean

Counts_edgeR_SP_7v3 <- DGEList(counts = corset_counts[,c(27:32, 21:26)], group = samples$Exposure_Time[c(27:32,21:26)])
expression_variation(Counts_edgeR_SP_7v3, "Counts_edgeR_SP_7v3.txt")
gamlss_NB_clean_SP7v3 <- gamlss_NB_clean