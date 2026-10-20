# Admixture
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