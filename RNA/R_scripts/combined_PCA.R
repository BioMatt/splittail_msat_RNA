library(tidyverse)
library(patchwork)

# First, reading in the neutral SNP data to plot the PCA. No PCA is run here, just collecting plots from the different datasets and putting them together.
load("LD_SNP_pca.RData")

LD_pca <- pca
rm(pca)

snps_var_frac[1] * 100
snps_var_frac[2] * 100

# In the neutral SNPs, 6.13% and 4.14% variance is explained for PC1 and PC2, respectively
rm(snps_var_frac)

LD_coords <- as_tibble(LD_pca$scores)

LD_coords$pop <- c(rep("CV", 16), rep("SP", 16))

neutral_plot <- ggplot(data = LD_coords, mapping = aes(x = PC1, y = PC2, group = pop, colour = pop)) +
  geom_point(size = 2) +
  labs(x = "PC1 6.13% Variance Explained", y = "PC2 4.14% Variance Explained", colour = "Population") +
  scale_colour_manual(values = c(CV = "#4575b4", SP = "#d73027"), labels = c("Central Valley", "San Pablo")) +
  theme_bw() +
  ggtitle("Neutral SNPs", subtitle = expression('69,951 SNPs')) +
  theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(colour = "black"),
                    axis.ticks = element_blank(),
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    text=element_text(size=16,  family="serif"),
                    legend.position="none"
  )
neutral_plot
# Following the same process for the overall SNP data (not pruned for LD or for HWE)
load("Overall_SNP_pca.RData")
overall_pca <- pca
rm(pca)

snps_var_frac[1] * 100
snps_var_frac[2] * 100

# In the overall SNPs, 5.28% and 4.01% variance is explained for PC1 and PC2, respectively
rm(snps_var_frac)
overall_coords <- as_tibble(overall_pca$scores)
overall_coords$pop <- c(rep("CV", 16), rep("SP", 16))

overall_plot <- ggplot(data = overall_coords, mapping = aes(x = PC1, y = PC2, group = pop, colour = pop)) +
  geom_point(size = 2) +
  labs(x = "PC1 5.28% Variance Explained", y = "PC2 4.01% Variance Explained", colour = "Population") +
  scale_colour_manual(values = c(CV = "#4575b4", SP = "#d73027"), labels = c("Central Valley", "San Pablo")) +
  theme_bw() +
  ggtitle("Overall SNPs", subtitle = expression('420,626 SNPs')) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text=element_text(size=16,  family="serif"),
        legend.position="none"
  )
overall_plot
# Again loading PCA data, but it's in a different format since it was done with ind.pca in hierfstat as opposed to glPca above
load("Msat_pca.RData")
msat_pca <- pca
rm(pca)

plot(msat_pca$ipca$li[,1], msat_pca$ipca$li[,2])

msat_coords <- as_tibble(msat_pca$ipca$li)
msat_pop <- msat_pca$ipca$rownames
msat_pop[msat_pop=="1"] <- "CV"
msat_pop[msat_pop=="2"] <- "SP"

msat_coords$pop <- msat_pop
rm(msat_pop)
# This line shows that 5.03% and 2.79% variance is explained by PC1 and PC2, respectively
msat_pca$ipca$eig

msat_plot <- ggplot(data = msat_coords, mapping = aes(x = Axis1, y = Axis2, group = pop, colour = pop)) +
  geom_point(size = 2) +
  labs(x = "PC1 5.03% Variance Explained", y = "PC2 2.79% Variance Explained", colour = "Population") +
  scale_colour_manual(values = c(CV = "#4575b4", SP = "#d73027"), labels = c("Central Valley", "San Pablo")) +
  theme_bw() +
  ggtitle("Microsatellites", subtitle = expression('19 Markers')) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text=element_text(size=16,  family="serif"),
        legend.justification = c("left", "bottom"),
        legend.position = c(0.01, 0.01),
        legend.background = element_blank()
  ) 
msat_plot
# Using the patchwork package to put together the three PCAs into one figure
combined_plot <- (msat_plot  | (neutral_plot / overall_plot)) 
ggsave("combined_PCA.pdf", combined_plot, dpi = 2000)
#unlink("combined_PCA.pdf")
