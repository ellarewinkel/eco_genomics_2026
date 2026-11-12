# All of this is the scripts copy-pasted



library(tidyverse)
setwd("~/projects/eco_genomics_2026/population_genomics/mydata")


# Now read in the PCA results from plink:

# the .eigenvec file has no column names, so col_names = FALSE (they become X1, X2, ...)
pca <- read_table("Lg_LD_0.2.eigenvec", col_names = FALSE)

# the .eigenval file is just a list of numbers, so scan() reads it into a vector
eigenval <- scan("Lg_LD_0.2.eigenval")

# columns 1 and 2 are the family and individual IDs, which are identical
# thanks to --double-id, so remove column 2 (remember: [rows, columns])
pca <- pca[, -2]

# rename the columns: the first is the sample ID, the rest are PC1, PC2, ...
# paste0() joins text and numbers, e.g. paste0("PC", 1:3) gives "PC1" "PC2" "PC3"
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

head(pca)



# How much variation does each PC explain?
# First let's see how much of the variation each PC explains (helps decide # of PCs worth looking at)
# convert the eigenvalues to a percentage variance explained (pve)
pve <- data.frame(PC = seq_along(eigenval),
                  pve = eigenval / sum(eigenval) * 100)

# make a folder for our figures

pve_plot <- ggplot(pve, aes(PC, pve)) +
  geom_col() +
  labs(x = "Principal component", y = "Percentage variance explained") +
  theme_bw()
pve_plot

OUTDIR="/users/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myresults/population_structure"

ggsave(pve_plot, file=paste0(OUTDIR, "/variance_0.2.png"), width = 6, height = 4)
# This is a scree plot




# To colour the PCA by site we need to join it to our sample metadata.
# Joining matches the rows of the two tables using the ID column, so the IDs in sample_metadata.csv must match the sample names in the VCF exactly.

meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (4.2%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()

OUTDIR="/users/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myresults/population_structure"

ggsave(pca_plot, file = paste0(OUTDIR, "/Lg_0.2_pca.png"), width = 8, height = 6)
# This is a colored PCA plot




# Is there a pattern with latitude?
# To plot PC1 and PC2 against latitude side by side, we first reshape the data with pivot_longer(). Instead of one column per PC, we get one row per individual per PC, which lets facet_wrap() give each PC its own panel:

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

OUTDIR="/users/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myresults/population_structure"

ggsave(lat_plot, file = paste0(OUTDIR, "/Lg_0.2_pca_latitude.png"), width = 8, height = 6)





# Bonus: A map of our sites
# We can also map our sites and put the map next to the PCA (might be useful for your project!)
# map_data() needs the maps package to be installed (we don’t need to load it with library())
# patchwork lets us combine ggplots
# If either isn’t installed, run install.packages(c("maps", "patchwork")) first.
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

# with patchwork, + puts plots side by side
# we hide the map's legend, as the PCA's legend already shows the sites
map_pca <- map_plot + theme(legend.position = "none") + pca_plot
map_pca

OUTDIR="/users/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myresults/population_structure"

ggsave(map_pca, file = paste0(OUTDIR, "/Lg_0.2_pca_map.png"), width = 8, height = 6)








###############################################################
ADMIXTURE now

# Do not run this on the full VCF!
# ADMIXTURE takes a long time to run even with batchjobs, so it has already been run for us

# Note on our data:
# For running admixture, K is defined for the min/max populations estimated to run (it looks like a max of 5 populations)



# We're not running the PCA and admixture as a batch job right now, but might be helpful to do after class
# when doing this and creating bash script, it says cd ~/popgen/class 2 in the code.
# however, you want to add your full working directory instead (to _)

library(tidyverse)
setwd("~/projects/eco_genomics_2026/population_genomics/mydata")

cv <- read_tsv("cv_errors.tsv", col_names = c("K", "CV_error"))

cv_plot <- ggplot(cv, aes(K, CV_error)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = cv$K) +
  labs(x = "K (number of ancestral populations)", y = "Cross-validation error") +
  theme_bw()

OUTDIR="~/projects/eco_genomics_2026/population_genomics/myresults/population_structure"

ggsave(cv_plot, file = paste0(OUTDIR, "/cv_error.png"), width = 6, height = 4)




# Plotting ancestry proportions
fam  <- read_table("Lg_LD_0.2.fam", col_names = FALSE)
meta <- read_csv("sample_metadata.csv")

# a function that reads the .Q file for one value of K and reshapes it
# into one row per individual per cluster
read_q <- function(K) {
  q <- read_table(paste0("Lg_LD_0.2.", K, ".Q"), col_names = FALSE)
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

OUTDIR="~/projects/eco_genomics_2026/population_genomics/myresults/population_structure"

ggsave(admix_plot, file = paste0(OUTDIR, "/admixture_K1-5.png"), width = 12, height = 7)
# Looking at this admixture (which is by latitude: north/south), we could surmise that...
# this agrees with our PCA. K=1 and K=3, which were outliers/genetically diff, are actually separate from each other(?)
# K=5 has very bad error rate (too spiky/random forced pink/messy), so data may be inaccurate
# can gage K=2 and K=3 for population structure
# in a perfect world, K=3 (3 pops aka 3 colors) would be completely green in one area, pea green in the other section, and all red in the other
# consistently across K=2 to K=5, the colors in Bahlia Asuncion for example are different
# Why are these groups distinctively different? Is it environmentally driven? We will get to this after looking into genetic diversity.
# The Southernmost groups are most different