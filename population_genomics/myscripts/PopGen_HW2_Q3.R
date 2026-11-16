# open vacc. nano?
# module load plink
# set $DATA_DIR and VCF
# prune for less strict (0.2 -> 0.5) LD threshold
# run PCA with new LD threshold

# head -n 3:
# ASU01 ASU01 -0.0253265 -0.100784 -0.0462359 0.0543123
# ASU02 ASU02 -0.0200904 -0.102493 -0.0747375 0.0205694
# ASU03 ASU03 -0.0539653 -0.10024 -0.065859 -0.0278861

# cat Lg_LD_0.5.eigenval
# 6.52459
# 3.80923
# 1.92898
# 1.62068
# 1.43148
# 1.38745
# 1.34882
# 1.32945
# 1.31227
# 1.28381
# 1.25412
# 1.22957
# 1.21032
# 1.19421
# 1.17783
# 1.17408
# 1.16319
# 1.15677
# 1.14634
# 1.1327

library(tidyverse)
setwd("~/projects/eco_genomics_2026/population_genomics/mydata")
OUTDIR="/gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myresults"

# Reading in PCA results from plink:
# the .eigenvec file has no column names, so col_names = FALSE (they become X1, X2, ...)
pca <- read_table("Lg_LD_0.5.eigenvec", col_names = FALSE)

# Use scan() to read the .eigenval file into a vector instead of the list of numbers
eigenval <- scan("Lg_LD_0.5.eigenval")

# Remove column 2 (remember: [rows, columns]) since family and individual IDs are identical, and Column 1 has the values already.
pca <- pca[, -2]

# rename the columns: the first is the sample ID, the rest are PC1, PC2, ...
# use paste0() to join text and numbers
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

head(pca)

# How much of the variation does each PC explain? (# of PCs worth looking at)
# convert the eigenvalues to a percentage variance explained (pve)
pve <- data.frame(PC = seq_along(eigenval),
                  pve = eigenval / sum(eigenval) * 100)

# create folder for figures
pve_plot <- ggplot(pve, aes(PC, pve)) +
  geom_col() +
  labs(x = "Principal component", y = "Percentage variance explained") +
  theme_bw()
pve_plot

ggsave(pve_plot, file=paste0(OUTDIR, "/variance_0.5.png"), width = 6, height = 4)
# Save the scree plot


# Join the PCA to sample metadata (match names of IDs in sample_metadata.csv and names in VCF) to color by site
meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0 (and did)
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (_%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()
pca_plot

ggsave(pca_plot, file = paste0(OUTDIR, "/Lg_0.5_pca.png"), width = 8, height = 6)
# Save the PCA plot colored by site

# Is there a pattern with latitude?
# Reshape data to give each PC its own panel
# Plot PC1 and PC2 against latitude side by side
lat_plot <- pca %>%
  select(ID, Site, Latitude, PC1, PC2) %>% ##Select columnns to use
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "score") %>%
  ggplot(aes(Latitude, score, colour = Site)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~PC, ncol = 1, scales = "free_y") + ##plots PC1 and PC2 seperately
  scale_colour_viridis_d(end = 0.9) + 
  labs(x = "Latitude (°N)", y = "PC score", colour = "Site (north to south)") +
  theme_bw()+
  coord_flip()
lat_plot

ggsave(lat_plot, file = paste0(OUTDIR, "/Lg_0.5_pca_latitude.png"), width = 8, height = 6)



# Map of the sites

library(patchwork) ##combined plots

world <- map_data("world")
sites <- pca %>% distinct(Site, Latitude, Longitude)

map_plot <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey90", colour = "grey60") +
  geom_point(data = sites, aes(Longitude, Latitude, colour = Site), size = 3) +
  scale_colour_viridis_d(end = 0.9) +
  # zoom to our sites, with a 2-degree margin all round
  coord_quickmap(xlim = range(sites$Longitude) + c(-2, 2),
                 ylim = range(sites$Latitude) + c(-2, 2)) +
  labs(x = "Longitude", y = "Latitude", colour = "Site (north to south)") +
  theme_bw()
map_plot

# with patchwork, + puts plots side by side
# we hide the map's legend since PCA's legend already shows the sites
map_pca <- map_plot + theme(legend.position = "none") + pca_plot
map_pca

ggsave(map_pca, file = paste0(OUTDIR, "/Lg_0.5_pca_map.png"), width = 8, height = 6)



# ADMIXTURE:
# do not run on full VCF!
# For running admixture, K is defined for the min/max populations estimated to run (it looks like a max of 5 populations)
# when doing this and creating bash script, it says cd ~/popgen/class 2 in the code.
# however, you want to add your full working directory instead (to _)

# run bash script? wasn't working. won't generate Q values
cp Lg_LD_0.5.bim Lg_LD_0.5.bim.orig
awk 'BEGIN{OFS="\t"} {$1="0"; print}' Lg_LD_0.5.bim.orig > Lg_LD_0.5.bim
#...? preparing the input

cv <- read_tsv("cv_errors.tsv", col_names = c("K", "CV_error"))

cv_plot <- ggplot(cv, aes(K, CV_error)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = cv$K) +
  labs(x = "K (number of ancestral populations)", y = "Cross-validation error") +
  theme_bw()
cv_plot

ggsave(cv_plot, file = paste0(OUTDIR, "/cv_error_hw2.png"), width = 6, height = 4)

# Plotting ancestry proportions
fam  <- read_table("Lg_LD_0.5.fam", col_names = FALSE)
meta <- read_csv("sample_metadata.csv")

# a function that reads the .Q file for one value of K and reshapes it
# into one row per individual per cluster
read_q <- function(K) {
  q <- read_table(paste0("Lg_LD_0.5.", K, ".Q"), col_names = FALSE)
  names(q) <- paste0("Cluster", 1:K)
  q %>%
    mutate(ID = fam$X2, K = K) %>%
    pivot_longer(starts_with("Cluster"), names_to = "Cluster", values_to = "Ancestry")
}

# map_dfr() runs read_q() for each K and binds the results together row-wise
# (change 2:5 to look at other values of K)
q_all <- map_dfr(1:5, read_q) %>%
  left_join(meta, by = "ID") %>%
  # order sites north to south, as in the PCA plots
  mutate(Site = fct_reorder(Site, -Latitude),
         K = paste0("K = ", K))
########## ERROR HERE IN map() because Lg_LD_0.5.1.Q does not exist. do i need to run

# check every individual matched the metadata: this should print 0
sum(is.na(q_all$Site))

admix_plot <- ggplot(q_all, aes(ID, Ancestry, fill = Cluster)) +
  geom_col(width = 1) +
  facet_grid(K ~ Site, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_minimal() +
  theme(axis.text.x     = element_blank(),         # too many sample names to read
        panel.grid      = element_blank(),
        panel.spacing.x = unit(0.1, "lines"),      # thin gap between sites
        panel.spacing.y = unit(0.8, "lines"),      # gap between values of K
        strip.text.x    = element_text(angle = 90, hjust = 0),
        legend.position = "none")

ggsave(admix_plot, file = paste0(OUTDIR, "/admixture_K1-5_hw2.png"), width = 12, height = 7
# Notes on new admixture:
# ...

### THESE NOTES ARE FROM THE IN-CLASS ASSIGNMENT! OLD ADMIXTURE:
# Looking at this admixture (which is by latitude: north/south), we could surmise that...
# this agrees with our PCA. K=1 and K=3, which were outliers/genetically diff, are actually separate from each other(?)
# K=5 has very bad error rate (too spiky/random forced pink/messy), so data may be inaccurate
# can gage K=2 and K=3 for population structure
# in a perfect world, K=3 (3 pops aka 3 colors) would be completely green in one area, pea green in the other section, and all red in the other
# consistently across K=2 to K=5, the colors in Bahlia Asuncion for example are different
# Why are these groups distinctively different? Is it environmentally driven? We will get to this after looking into genetic diversity.
# The Southernmost groups are most different