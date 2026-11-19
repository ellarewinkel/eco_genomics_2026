# set to OWN working directory
rm(list=ls())
setwd("~/projects/eco_genomics_2026/final_project/mydata")

# Install R packages
.libPaths(c("/users/d/s/dsadler1/programmes/rpackages_441", .libPaths()))
libpath <- .libPaths()[1]
.libPaths() 

# setup environment
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
library(phyloseq)

library(BiocManager)
BiocManager::install("microbiome")
library(microbiome)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
library(decontam)

library(ANCOMBC) #v2.9.1 #didn't run yet (must install)
library(phangorn) #to convert the midpoint rooted tree
library(MicrobiomeStat)
library(picante)
library(ape)
library(qiime2R) #didn't run yet (must install)
#library(microbiomeutilities)

library(ggplot2)
library(ggpubr)
library(ggsci)
#library(ggprism)
#library(gghalves)

library(vegan)
library(tidyverse) 
library(devtools)
library(knitr)
library(plyr)
library(dplyr)
library(rstatix)

library(RColorBrewer)

#library(scico)

# IMPORT OTU TABLE
otu <- read.delim(
  "otu_table.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

otu <- as.matrix(otu)

# create phyloseq OTU object
OTU <- otu_table(
  otu,
  taxa_are_rows = TRUE
)

# IMPORT TAXONOMY TABLE
tax <- read.delim(
  "taxonomy_1.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# OTUID column contains the same IDs as the rows of the OTU table -> make OTUID the row names
rownames(tax) <- tax$OTUID

# remove OTUID from the taxonomy columns
tax$OTUID <- NULL

# convert to a phyloseq taxonomy table
TAX <- tax_table(
  as.matrix(tax)
)

# IMPORT METADATA
metadata <- read.table(
  "METADATA_HYP_1.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# create phyloseq sample-data object
MET <- sample_data(metadata)

# IMPORT PHYLOGENETIC TREE
phy_tree <- read.tree("tree.nwk")

# 5. CHECK THAT THE NAMES MATCH
cat("DATA DIMENSIONS\n")

cat("OTU table:", nrow(otu), "taxa x", ncol(otu), "samples\n") # 36756 taxa x 473 samples
cat("Taxonomy:", nrow(tax), "taxa\n") # 73171 taxa
cat("Metadata:", nrow(metadata), "samples\n") # 477 samples
cat("Tree:", length(phy_tree$tip.label), "tips\n") # 36756 tips

cat("NAME MATCHING\n")

cat(
  "OTU taxa found in taxonomy:",
  sum(taxa_names(OTU) %in% taxa_names(TAX)),
  "of",
  ntaxa(OTU),
  "\n"
)
# OTU taxa found in taxonomy: 36756 of 36756

cat(
  "OTU taxa found in tree:",
  sum(taxa_names(OTU) %in% phy_tree$tip.label),
  "of",
  ntaxa(OTU),
  "\n"
)
# OTU taxa found in tree: 36756 of 36756 

cat(
  "OTU samples found in metadata:",
  sum(sample_names(OTU) %in% rownames(metadata)),
  "of",
  nsamples(OTU),
  "\n"
)
# OTU samples found in metadata: 473 of 473 

################################################################################
# KEEP ONLY TAXA PRESENT IN ALL COMPONENTS
# find taxa shared by OTU table, taxonomy and tree
#common_taxa <- Reduce(
#  intersect,
#  list(
#    taxa_names(OTU),
#    taxa_names(TAX),
#    phy_tree$tip.label
#  )
#)

#cat("\nCommon taxa:", length(common_taxa), "\n")

# keep common taxa
#OTU <- prune_taxa(
#  common_taxa,
#  OTU
#)

#TAX <- prune_taxa(
#  common_taxa,
#  TAX
#)

#phy_tree <- prune_taxa(
#  common_taxa,
#  phy_tree
#)

# KEEP ONLY SAMPLES PRESENT IN BOTH OTU TABLE AND METADATA
#common_samples <- intersect(
#  sample_names(OTU),
#  rownames(metadata)
#)

#cat("Common samples:", length(common_samples), "\n")

#OTU <- prune_samples(
#  common_samples,
#  OTU
#)

#MET <- prune_samples(
#  common_samples,
#  MET
#)
################################################################################


# CREATE PHYLOSEQ OBJECT
ps <- phyloseq(
  OTU,
  TAX,
  MET,
  phy_tree
)

# CHECK
cat("PHYLOSEQ OBJECT CREATED\n")
print(ps)

cat("\nNumber of taxa:", ntaxa(ps), "\n")
cat("Number of samples:", nsamples(ps), "\n")
cat("\nFirst 10 taxon names:\n")
print(head(taxa_names(ps), 10))

cat("\nFirst 10 sample names:\n")
print(head(sample_names(ps), 10))

cat("\nTaxonomy:\n")
print(tax_table(ps)[1:5, ])

microbiome::summarize_phyloseq(ps)

################################################################################

# CLEAN PHYLOSEQ TAXONOMY TABLE
# convert phyloseq taxonomy table to data frame
tax <- data.frame(
  tax_table(ps),
  stringsAsFactors = FALSE
)

# make sure everything is character
tax[] <- lapply(tax, as.character)

# remove prefixes (e.g. d__, p__, c__, o__, f__, g__, s__)
tax[] <- lapply(
  tax,
  function(x) gsub("^[a-z]__", "", x)
)

# remove "uncultured..." / "unculture..." annotations
tax[] <- lapply(
  tax,
  function(x) gsub("uncultured.*", "", x, ignore.case = TRUE)
)

# convert NA to blank
tax[is.na(tax)] <- ""

# remove leading/trailing whitespace
tax[] <- lapply(
  tax,
  trimws
)

# fill missing lower taxonomy with "unclassified_PARENT"

for (i in seq_len(nrow(tax))) {
  if (tax[i, "Phylum"] == "") {
    
    tax[i, "Phylum"]  <- paste0("unclassified_", tax[i, "Kingdom"])
    tax[i, "Class"]   <- paste0("unclassified_", tax[i, "Kingdom"])
    tax[i, "Order"]   <- paste0("unclassified_", tax[i, "Kingdom"])
    tax[i, "Family"]  <- paste0("unclassified_", tax[i, "Kingdom"])
    tax[i, "Genus"]   <- paste0("unclassified_", tax[i, "Kingdom"])
    tax[i, "Species"] <- paste0("unclassified_", tax[i, "Kingdom"])
    
  } else if (tax[i, "Class"] == "") {
    tax[i, "Class"]   <- paste0("unclassified_", tax[i, "Phylum"])
    tax[i, "Order"]   <- paste0("unclassified_", tax[i, "Phylum"])
    tax[i, "Family"]  <- paste0("unclassified_", tax[i, "Phylum"])
    tax[i, "Genus"]   <- paste0("unclassified_", tax[i, "Phylum"])
    tax[i, "Species"] <- paste0("unclassified_", tax[i, "Phylum"])
    
  } else if (tax[i, "Order"] == "") {
    tax[i, "Order"]   <- paste0("unclassified_", tax[i, "Class"])
    tax[i, "Family"]  <- paste0("unclassified_", tax[i, "Class"])
    tax[i, "Genus"]   <- paste0("unclassified_", tax[i, "Class"])
    tax[i, "Species"] <- paste0("unclassified_", tax[i, "Class"])
    
  } else if (tax[i, "Family"] == "") {
    tax[i, "Family"]  <- paste0("unclassified_", tax[i, "Order"])
    tax[i, "Genus"]   <- paste0("unclassified_", tax[i, "Order"])
    tax[i, "Species"] <- paste0("unclassified_", tax[i, "Order"])
    
  } else if (tax[i, "Genus"] == "") {
    tax[i, "Genus"]   <- paste0("unclassified_", tax[i, "Family"])
    tax[i, "Species"] <- paste0("unclassified_", tax[i, "Family"])
    
  } else if (tax[i, "Species"] == "") {
    tax[i, "Species"] <- paste0("unclassified_", tax[i, "Genus"])
  }
}

# return cleaned taxonomy -> phyloseq
tax_table(ps) <- tax_table(
  as.matrix(tax)
)

# check
tax_table(ps)[1:10, ]
ntaxa(ps) # [1] 36756
nsamples(ps) # [1] 473

####################################################################################
## DECONTAM - prevalence method
####################################################################################

# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# identify contaminants - prevalence
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
head(which(contamdf.prev$contaminant))
contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
contamdf.prev04 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.4)
contamdf.prev02 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.2)

# TRUE = number of sequences listed as a contaminant based on decontaminant prevalence threshold
# threshold 0.5 is null 5, 0.4 is null 4, and 0.2 is null 2 (higher value removes more taxa)
table(contamdf.prev05$contaminant)
# FALSE 36728, TRUE 28

table(contamdf.prev04$contaminant)
# FALSE 36728, TRUE 28

table(contamdf.prev02$contaminant)
# FALSE 36731, TRUE 25

ps.decontam <- prune_taxa(!contamdf.prev05$contaminant, ps)
ps.decontam

ps.decontam
sample_sums(ps.decontam)
sort(sample_sums(ps.decontam))

# remove blanks
ps.clean <- subset_samples(ps.decontam, Sample_or_Control == "Sample")
ps.clean
# 36728 taxa, 465 samples

save(ps.clean,file="ps.clean.Rdata")

load("ps.clean.Rdata")

# skin samples
ps.skin <- subset_samples(ps.decontam, SampleType == "Skin")
ps.skin

save(ps.skin,file="ps.skin.Rdata")
load("ps.skin.Rdata")

# 36728 taxa, 435 samples

################################################################################

# RAREFACTION DEPTH ASSESSMENT
# BASIC INFORMATION 
cat("Number of samples:", nsamples(ps.skin), "\n") # 435 samples
cat("Number of features:", ntaxa(ps.skin), "\n") # 36728 features

reads <- sample_sums(ps.skin)

cat("Minimum reads/sample:", min(reads), "\n") # 33 minimum reads
cat("Median reads/sample:", median(reads), "\n") # 39395 median reads
cat("Maximum reads/sample:", max(reads), "\n") # 76697 max reads 

# DISTRIBUTION OF READ DEPTHS
df_reads <- data.frame(
  sample = names(reads),
  reads = as.numeric(reads)
)

library(tidyverse)

p1<-ggplot(df_reads, aes(x = reads)) +
  geom_histogram(
    bins = 25,
    fill = "#A6C9D9",
    colour = "black"
  ) +
  labs(
    x = "\nNumber of reads",
    y = "Frequency\n"
  ) +
  theme_classic(base_size = 18)
p1

ggsave(p1, file = "../myresults/rare.png")

# We want to cut off samples with low reads because they're less accurate.
# Choose which to cut out

# CANDIDATE RAREFACTION DEPTHS
# NB: a depth above  minimum library size cannot be used for rarefaction (without losing samples).
depths <- c(
  1000,
  2500,
  5000,
  10000,
  15000,
  20000,
  30000,
  60000
)

# Decided cutoff at 10,000 depth since we'll lose 26 samples.

# SAMPLE RETENTION AT EACH DEPTH
sample_retention <- data.frame(
  depth = depths,
  samples_retained = sapply(
    depths,
    function(d) sum(reads >= d)
  )
) %>%
  mutate(
    samples_total = nsamples(ps.skin),
    proportion_retained = samples_retained / samples_total,
    percent_retained = 100 * proportion_retained,
    samples_lost = samples_total - samples_retained
  )

cat("SAMPLE RETENTION\n")
print(sample_retention)

# PLOT SAMPLE RETENTION
ggplot(
  sample_retention,
  aes(x = depth, y = percent_retained)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 80,
    linetype = "dashed",
    colour = "orange"
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "darkgreen"
  ) +
  scale_y_continuous(
    limits = c(0, 100)
  ) +
  labs(
    x = "\nRarefaction depth (reads/sample)",
    y = "Samples retained (%)\n"
  ) +
  theme_classic(base_size = 18)

# Rarefaction curve above: 0 samples at 60,000 depth. 90% of samples remain around 10,000 depth

# CONVERT OTU TABLE TO SAMPLE x FEATURE MATRIX
otu <- as(otu_table(ps.skin), "matrix")

if (taxa_are_rows(ps.skin)) {
  otu <- t(otu)
}

# Check dimensions
cat("\nOTU matrix dimensions:\n")
cat("Samples:", nrow(otu), "\n")
cat("Features:", ncol(otu), "\n")

# RAREFACTION DEPTH ASSESSMENT
# GET OTU TABLE AS SAMPLE x FEATURE MATRIX
otu <- as(otu_table(ps.skin), "matrix")

if (taxa_are_rows(ps.skin)) {
  otu <- t(otu)
}

# check dimensions
cat("DATASET\n")
cat("Samples:", nrow(otu), "\n")
cat("Features:", ncol(otu), "\n")

# CHECK LIBRARY SIZES
sample_depths <- rowSums(otu)

cat("LIBRARY SIZE SUMMARY\n")
cat("Minimum reads:", min(sample_depths), "\n")
cat("1st quartile:", quantile(sample_depths, 0.25), "\n")
cat("Median reads:", median(sample_depths), "\n")
cat("3rd quartile:", quantile(sample_depths, 0.75), "\n")
cat("Maximum reads:", max(sample_depths), "\n")

cat("\nLowest 26 samples:\n") # Lowest 26 based on rarefaction score we picked (10,000)
print(sort(sample_depths)[1:min(20, length(sample_depths))])

# CANDIDATE RAREFACTION DEPTHS
# NB: rarefaction CURVES extend toward the maximum library size.
# RETENTION analysis tells how many samples  have enough reads for each depth
# final maximum is included even if it is not exactly divisible by interval

#max_depth <- max(sample_depths)
max_depth <- max(sample_depths/2)

curve_depths <- seq(
  from = 100,
  to = max_depth,
  by = 500
)

curve_depths <- sort(
  unique(
    c(
      curve_depths,
      max_depth
    )
  )
)

cat("RAREFACTION CURVE DEPTHS\n")
cat(
  "Curve range:",
  min(curve_depths),
  "to",
  max(curve_depths),
  "reads\n"
)

cat(
  "Number of curve depths:",
  length(curve_depths),
  "\n"
)

# CALCULATE RAREFACTION CURVES
# samples with fewer reads than a particular depth will have NA at that depth.
rare_df <- lapply(
  rownames(otu),
  function(samp) {
    
    richness <- sapply(
      curve_depths,
      function(d) {
        
        if (sum(otu[samp, ]) >= d) {
          
          vegan::rarefy(
            otu[samp, , drop = TRUE],
            sample = d
          )
          
        } else {
          
          NA_real_
        }
      }
    )
    
    data.frame(
      sample = samp,
      depth = curve_depths,
      richness = richness
    )
  }
) %>%
  bind_rows()

# PLOT RAREFACTION CURVES
ggplot(
  rare_df,
  aes(
    x = depth,
    y = richness,
    group = sample
  )
) +
  geom_line(
    alpha = 0.25,
    linewidth = 0.4,
    na.rm = TRUE
  ) +
  labs(
    x = "\nRead depth",
    y = "Expected observed features\n"
  ) +
  theme_classic(base_size = 18)

# SAMPLE RETENTION AT CANDIDATE DEPTHS
depths <- sort(
  unique(
    c(
      seq(
        1000,
        max_depth,
        by = 1000
      ),
      max_depth
    )
  )
)

sample_retention <- data.frame(
  depth = depths,
  samples_retained = sapply(
    depths,
    function(d) {
      sum(sample_depths >= d)
    }
  )
) %>%
  mutate(
    samples_total = nrow(otu),
    samples_lost = samples_total - samples_retained,
    percent_retained =
      100 * samples_retained / samples_total
  )

cat("SAMPLE RETENTION\n")
print(sample_retention)

# PLOT SAMPLE RETENTION
ggplot(
  sample_retention,
  aes(
    x = depth,
    y = percent_retained
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    colour = "darkgreen"
  ) +
  geom_hline(
    yintercept = 80,
    linetype = "dashed",
    colour = "orange"
  ) +
  scale_y_continuous(
    limits = c(0, 100)
  ) +
  labs(
    x = "\nRarefaction depth (reads/sample)",
    y = "Samples retained (%)\n"
  ) +
  theme_classic(base_size = 18)

# FEATURE RETENTION AFTER RAREFACTION
# NB: feature retention is stochastic because rarefaction randomly selects reads. Use a fixed seed -> reproduciblity.
set.seed(1)

feature_retention <- lapply(
  depths,
  function(d) {
    
    # Keep only samples with enough reads
    keep <- sample_depths >= d
    
    otu_keep <- otu[
      keep,
      ,
      drop = FALSE
    ]
    
    # If no samples remain, return NA
    if (nrow(otu_keep) == 0) {
      
      return(
        data.frame(
          depth = d,
          features_retained = NA_integer_
        )
      )
    }
    
    # Rarefy every retained sample
    otu_rare <- t(
      apply(
        otu_keep,
        1,
        function(x) {
          
          vegan::rrarefy(
            matrix(
              x,
              nrow = 1
            ),
            sample = d
          )[1, ]
        }
      )
    )
    
    # Make sure dimensions remain correct if only one sample
    if (nrow(otu_keep) == 1) {
      otu_rare <- matrix(
        otu_rare,
        nrow = 1
      )
    }
    
    # Number of features represented by >=1 read
    n_features <- sum(
      colSums(otu_rare) > 0
    )
    
    data.frame(
      depth = d,
      features_retained = n_features
    )
  }
) %>%
  bind_rows()

head(rare_df)
# if we were to cut the depth, how rich is the data?
#    sample depth  richness
# 1 T1.BB.A   100  32.09711
# 2 T1.BB.A   600  87.61957
# 3 T1.BB.A  1100 116.29874
# 4 T1.BB.A  1600 135.60161
# 5 T1.BB.A  2100 149.81670
# 6 T1.BB.A  2600 160.83206

# COMBINE SAMPLE + FEATURE RETENTION
retention_summary <- sample_retention %>%
  left_join(
    feature_retention,
    by = "depth"
  ) %>%
  mutate(
    original_features = ntaxa(ps.skin),
    
    feature_percent_retained =
      100 * features_retained /
      original_features,
    
    feature_loss =
      original_features -
      features_retained
  )
ggplot

cat("RETENTION SUMMARY\n")
print(retention_summary)

# DEPTHS RETAINING >=90% OF SAMPLES
acceptable_depths_90 <- retention_summary %>%
  filter(
    percent_retained >= 90
  )

cat("DEPTHS RETAINING >=90% OF SAMPLES\n")
print(acceptable_depths_90)

# DEPTHS RETAINING >=80% OF SAMPLES
acceptable_depths_80 <- retention_summary %>%
  filter(
    percent_retained >= 80
  )

cat("DEPTHS RETAINING >=80% OF SAMPLES\n")
print(acceptable_depths_80)

# PLOT SAMPLE AND FEATURE RETENTION
ggplot(
  retention_summary,
  aes(x = depth)
) +
  
  geom_line(
    aes(
      y = percent_retained,
      colour = "Samples"
    ),
    linewidth = 1
  ) +
  
  geom_point(
    aes(
      y = percent_retained,
      colour = "Samples"
    ),
    size = 2.5
  ) +
  
  geom_line(
    aes(
      y = feature_percent_retained,
      colour = "Features"
    ),
    linewidth = 1
  ) +
  
  geom_point(
    aes(
      y = feature_percent_retained,
      colour = "Features"
    ),
    size = 2.5
  ) +
  
  geom_hline(
    yintercept = 80,
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  scale_y_continuous(
    limits = c(0, 100)
  ) +
  
  labs(
    x = "\nRarefaction depth (reads/sample)",
    y = "Retention (%)\n",
    colour = ""
  ) +
  
  theme_classic(base_size = 18)



# IDENTIFY THE HIGHEST DEPTH RETAINING >=80% OF SAMPLES
if (nrow(acceptable_depths_80) > 0) {
  
  candidate_depth_80 <- max(
    acceptable_depths_80$depth
  )
  
  cat("CANDIDATE DEPTH: >=80% SAMPLE RETENTION\n")

  print(
    retention_summary %>%
      filter(
        depth == candidate_depth_80
      )
  )
  
} else {
  
  cat(
    "\nNo candidate depth retains >=80% of samples.\n"
  )
}

# IDENTIFY THE HIGHEST DEPTH RETAINING >=90% OF SAMPLES
if (nrow(acceptable_depths_90) > 0) {
  
  candidate_depth_90 <- max(
    acceptable_depths_90$depth
  )
  
  cat("CANDIDATE DEPTH: >=90% SAMPLE RETENTION\n")

  print(
    retention_summary %>%
      filter(
        depth == candidate_depth_90
      )
  )
  
} else {
  
  cat(
    "\nNo candidate depth retains >=90% of samples.\n"
  )
}

# inspect retention_summary table and plots -> rarefaction depth
# examine library read depths to further refine rarefaction depth
sample_sums(ps.skin)
sort(sample_sums(ps.skin))

## 10000 -> 26 samples lost (we had 435) - can do _ or _
ps.skin.rare <- rarefy_even_depth(
   ps.skin,
   sample.size = 10000,
   replace = FALSE,
   rngseed = 1,
   verbose = FALSE
 )
save.image(file = "ps.clean.RData")
load("ps.clean.RData")

#save.image(file = "ps.clean.RData")
#load("ps.clean.RData")

################################################################################
## COMMUNITY COMPOSITION
################################################################################

# abundance in %
# work on pooled samples
# Check that Treatment_main exists in ps.skin
sample_variables(ps.skin)

# See the treatment assigned to every sample
table(
  sample_data(ps.skin)$Treatment_main,
  useNA = "ifany"
)

# Check how many samples you have before pooling
nsamples(ps.skin)

# Pool samples
ps.pool <- merge_samples(
  ps.skin,
  "Treatment_main"
)

# Check how many samples remain
nsamples(ps.pool)

# Check the names of the pooled samples
sample_names(ps.pool)

# Check their metadata
sample_data(ps.pool)

# restore metadata
sample_data(ps.pool)$Treatment_main <- sample_names(ps.pool)
print(sample_data(ps.pool)$Treatment_main)
#  [1] "1_Hyp_25"  "1_Hyp_35"  "1_Norm_25" "1_Norm_35" "3_Hyp_25"  "3_Hyp_35"  "3_Norm_25"
#  [8] "3_Norm_35" "4_Hyp_25"  "4_Hyp_35"  "4_Norm_25" "4_Norm_35"

# check that there are correct number of pooled treatments
print(nsamples(ps.pool))
# 12 samples. This is correct

# check names
print(sample_names(ps.pool))

# convert to percentages
ps.pool.perc <- transform_sample_counts(ps.pool, function(x) x/sum(x))

# PHYLA
install.packages("scico")
library(ggplot2)
library(scico)
glom.p <- tax_glom(ps.pool.perc, taxrank = 'Phylum', NArm = FALSE)
glom.p # should list # taxa as # phyla
data_glom.p <- psmelt(glom.p) # create data frame from phyloseq object
data_glom.p$Phylum <- as.character(data_glom.p$Phylum) #convert to character

# rename phyla with <X% abundance
data_glom.p$Phylum[data_glom.p$Abundance < 0.05] <- "<5% abundance features"

# count # phyla to set colour palette
Count = length(unique(data_glom.p$Phylum))
Count # returns 7 phyla

unique(data_glom.p$Phylum)

# create a list of these Phyla
phyla.levels <- unique(data_glom.p$Phylum)

phyla.levels <- sort(
  phyla.levels[
    phyla.levels != "<5% abundance features"
  ]
)

phyla.levels <- c(
  phyla.levels,
  "<5% abundance features"
)

cat(
  'levels = c("',
  paste(
    phyla.levels,
    collapse = '",\n"'
  ),
  '")\n',
  sep = ""
)

# reorder alphabetically, with pooled samples last
data_glom.p$Phylum <- factor(data_glom.p$Phylum, 
                             levels = c("Actinobacteriota",
                                        "Bacillota",
                                        "Bacteroidota",
                                        "Deinococcota",
                                        "Fusobacteriota",
                                        "Proteobacteria",
                                        "Verrucomicrobiota",
                                        "<5% abundance features")
                             )

# plotting
# check the treatment combinations that actually exist
unique(data_glom.p$Treatment_main)

# ORDER THE POOLED TREATMENT GROUPS
treatment_order <- c(
  "1_Norm_25",
  "1_Norm_35",
  "1_Hyp_25",
  "1_Hyp_35",
  "3_Norm_25",
  "3_Norm_35",
  "3_Hyp_25",
  "3_Hyp_35",
  "4_Norm_25",
  "4_Norm_35",
  "4_Hyp_25",
  "4_Hyp_35"
)

data_glom.p$Treatment_main <- factor(
  data_glom.p$Treatment_main,
  levels = treatment_order
)

# PHYLA BAR PLOT
phyla <- ggplot(
  data = data_glom.p,
  aes(
    x = Treatment_main,
    y = Abundance,
    fill = Phylum
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.9
  ) +
  labs(
    x = "\nTreatment",
    y = "\nRelative abundance\n"
  ) +
  scale_fill_scico_d(
    palette = "vik"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 18
    ),
    axis.text.y = element_text(
      size = 18
    ),
    axis.title = element_text(
      size = 18
    ),
    legend.text = element_text(
      size = 18
    ),
    legend.title = element_text(
      size = 18
    ),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

phyla

# CLASS
glom.c <- tax_glom(ps.pool.perc, taxrank = 'Class', NArm = FALSE)
glom.c # should list # taxa as # class
data_glom.c <- psmelt(glom.c) # create data frame from phyloseq object
data_glom.c$Class <- as.character(data_glom.c$Class) # convert to character

# simple way to rename class with <4% abundance
data_glom.c$Class[data_glom.c$Abundance < 0.04] <- "<4% abundance features"

# count # class to set colour palette
Count.c = length(unique(data_glom.c$Class))
Count.c # returns 10 classes

unique(data_glom.c$Class)
sort(unique(as.character(data_glom.c$Class)))

# reorder alphabetically, with pooled samples last
class.levels <- sort(
  unique(
    as.character(data_glom.c$Class)
  )
)

class.levels <- c(
  class.levels[
    class.levels != "<4% abundance features"
  ],
  "<4% abundance features"
)

data_glom.c$Class <- factor(
  data_glom.c$Class,
  levels = class.levels
)

# plotting
class <- ggplot(data=data_glom.c, aes(x = factor(Sample, level=treatment_order), 
                                      y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", width=0.9) +
  labs(x = "\nsample\n", y = "\nrelative abundance\n") +
  scale_x_discrete(labels=treatment_order) + 
  scale_fill_scico_d(palette = "vik") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(angle = 360, hjust = 1, vjust = 0.5, size = 18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + #kills the facet headers (= strip)
  theme(panel.background = element_blank())

class

# ORDER
glom.o <- tax_glom(ps.pool.perc, taxrank = 'Order', NArm = FALSE)
glom.o # should list # taxa as # order
data_glom.o <- psmelt(glom.o) # create data frame from phyloseq object
data_glom.o$Order <- as.character(data_glom.o$Order) #convert to character

# simple way to rename order with <4% abundance
data_glom.o$Order[data_glom.o$Abundance < 0.04] <- "<4% abundance features"

# count # order to set color palette
Count.o = length(unique(data_glom.o$Order))
Count.o # returns 13 orders

unique(data_glom.o$Order)
sort(unique(as.character(data_glom.o$Order)))

# reorder alphabetically, with pooled samples last
order.levels <- sort(
  unique(
    as.character(data_glom.o$Order)
  )
)

order.levels <- c(
  order.levels[
    order.levels != "<4% abundance features"
  ],
  "<4% abundance features"
)

data_glom.o$Order <- factor(
  data_glom.o$Order,
  levels = order.levels
)

# plotting
order <- ggplot(data=data_glom.o, aes(x = factor(Sample, level=treatment_order), 
                                      y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", width=0.9) +
  labs(x = "\nsample\n", y = "\nrelative abundance\n") +
  scale_x_discrete(labels=treatment_order) + 
  #  scale_fill_scico_d(palette = "batlow") +
  scale_fill_scico_d(palette = "vik") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
        axis.text.y = element_text(angle = 360, hjust = 1, vjust = 0.5, size = 18),
        axis.title = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) + #kills the facet headers (= strip)
  theme(panel.background = element_blank())

order

#single plot
#plot_grid(
#  phyla + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()),
#  class + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()),
#  order,
#  ncol = 1,
#  align = 'v', 
#  rel_heights = c(1, 1, 1.6)  # equalize bar height
#)

#save.image(file = "bar.RData")
#load("bar.RData")


################################################################################
## ALPHA DIVERSITY
################################################################################

# get an overview of diversity in samples for all available indices of alpha diversity
alpha_table <- microbiome::alpha(ps.skin.rare, index = "all")
kable(head(alpha_table))

######################## to include FAITH's PD
# extract OTU table and transpose to samples × taxa
otu_table_pd <- as(otu_table(ps.skin.rare), "matrix")

# transpose for PICANTE
otu_table_for_pd <- t(otu_table_pd)

# extract phylogenetic tree
phy_tree_pd <- phy_tree(ps.skin.rare)
head(phy_tree_pd$tip.label)  # Should match rownames of `otu_table_pd`

# Ensure taxa names match
if (!all(rownames(otu_table_pd) %in% phy_tree_pd$tip.label)) {
  stop("Some taxa in the OTU table are missing from the phylogenetic tree")
}

# calculate Faith's PD
faith_pd <- pd(otu_table_for_pd, phy_tree_pd, include.root = TRUE)

# add to alpha_table
faith_pd_output <- as.data.frame(faith_pd)

# convert rownames (sample names) to a column
faith_pd_output <- faith_pd %>%
  rownames_to_column(var = "SampleID")

# convert alpha_table rownames to a column (if sample names are rownames there too)
alpha_table_output <- alpha_table %>%
  rownames_to_column(var = "SampleID")

# join by SampleID
alpha_table_output <- alpha_table_output %>%
  left_join(faith_pd_output, by = "SampleID")

# optional: convert SampleID back to rownames if you want
# remove existing rownames first by converting to a plain data frame
alpha_table <- alpha_table_output %>% 
  as.data.frame() %>%           # ensure it's a plain data frame (no rownames)
  remove_rownames() %>%        # remove rownames
  column_to_rownames(var = "SampleID")  # convert SampleID column to rownames

# prepare alpha_table with SampleID as a column
alpha_df <- alpha_table %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

# prepare meta with proper SampleID column name - follow this or the metadata as as a phyloseq object and not a data frame and dplyr cannot merge data
meta_df <- as(sample_data(ps.skin.rare), "data.frame")
class(meta_df)
colnames(meta_df)[colnames(meta_df) == "sample"] <- "SampleID"

# prepare meta with proper SampleID column name - follow this or the metadata stas as a phyloseq object and not a data frame and dplyr cannot merge data
# fix because 1st column was rowname and not sampleID
meta_df <- as(sample_data(ps.skin.rare), "data.frame")
class(meta_df)
meta_df <- meta_df %>%
  tibble::rownames_to_column(var = "SampleID")

# join into a new data frame for plotting
plot_data <- left_join(alpha_df, meta_df %>% 
                         select(SampleID, Treatment_main, Tank, Week, Ox, Temp, Temp_1), 
                       by = "SampleID")

#custom  colours
my_colors_12 <- scico(n = 12, palette = "vik")
print(my_colors_12)  # just to check
#"#001260" "#013678" "#095B90" "#4088AE" "#89B7CD" "#D1E0E8" "#EBD2C4"
#"#D9A587" "#C77B50" "#B04F1B" "#811F06" "#590007"

# CUSTOM colours for publication
# "#E7C6B2", "#C98157", "#9D3709",
# example
#color_map.0 <- c(
#  "sponge" = my_colors_6[1],
#  "sediment" = my_colors_6[3],
#  "water" = my_colors_6[6]
#)

color_map <- c(
  "1_Norm_25" = "#001260",
  "1_Norm_35" = "#A6C9D9",
  "1_Hyp_25" = "#590007"
)

obs <- ggplot(plot_data, aes(
#  x = factor(type_1, levels = names(color_map)),
  x = factor(Treatment_main),
    y = observed,
  fill = Treatment_main,
  color = Treatment_main)) +    # also map color for dots outline
#  scale_x_discrete(labels=c('sponge' = 'sponge', 
#                            'sediment' = 'sediment', 
#                            'water'  = 'water')) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 1.0, size = 2.2, shape = 21) +  # shape 21 has fill + border
#  scale_fill_manual(values = color_map) +
#  scale_color_manual(values = color_map) +  # keep dot border same as fill
  labs(x = "\nsample\n", y = "\nobserved features\n") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),   # removes grey grid
    axis.line = element_line(color = "black", size = 0.8),  # black axis lines with thickness 0.8
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

obs

shan <- ggplot(plot_data, aes(
  x = factor(Treatment_main),
  y = diversity_shannon,
  fill = Treatment_main,
  color = Treatment_main)) +    # also map color for dots outline
#  scale_x_discrete(labels=c('sponge' = 'sponge', 
#                            'sediment' = 'sediment', 
#                            'water'  = 'water')) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 1.0, size = 2.2, shape = 21) +  # shape 21 has fill + border
#  scale_fill_manual(values = color_map) +
#  scale_color_manual(values = color_map) +  # keep dot border same as fill
  labs(x = "\nsample\n", y = "\nShannon index\n") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),   # removes grey grid
    axis.line = element_line(color = "black", size = 0.8),  # black axis lines with thickness 0.8
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

shan

faith <- ggplot(plot_data, aes(
  x = factor(Treatment_main),
  y = PD,
  fill = Treatment_main,
  color = Treatment_main)) +    # also map color for dots outline
#  scale_x_discrete(labels=c('sponge' = 'sponge', 
#                            'sediment' = 'sediment', 
#                            'water'  = 'water')) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 1.0, size = 2.2, shape = 21) +  # shape 21 has fill + border
#  scale_fill_manual(values = color_map) +
#  scale_color_manual(values = color_map) +  # keep dot border same as fill
  labs(x = "\nsample\n", y = "\nFaith's PD\n") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),   # removes grey grid
    axis.line = element_line(color = "black", size = 0.8),  # black axis lines with thickness 0.8
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

faith

#save.image(file = "alpha.RData")
#load("alpha.RData")

## SIGNIFICANCE TESTING
## pairwise Kruskal Wallis
res.kruskal <- plot_data %>% kruskal_test(observed ~ Treatment_main)
res.kruskal

# pairwise comparisons
# DUNNS test
dunn_pwc <- plot_data %>% dunn_test(observed ~ Treatment_main, p.adjust.method = "bonferroni") 
print(n = 100, dunn_pwc)

# Shannon
res.kruskal.s <- plot_data %>% kruskal_test(diversity_shannon ~ Treatment_main)
res.kruskal.s
dunn_pwc.s <- plot_data %>% dunn_test(diversity_shannon ~ Treatment_main, p.adjust.method = "bonferroni") 
print(n = 100, dunn_pwc.s)

# Faith
res.kruskal.pd <- plot_data %>% kruskal_test(PD ~ Treatment_main)
res.kruskal.pd
dunn_pwc.pd <- plot_data %>% dunn_test(PD ~ Treatment_main, p.adjust.method = "bonferroni") 
print(n = 100, dunn_pwc.pd)

#save.image(file = "enclosure.bar.alpha.97.RData")
#load("enclosure.bar.alpha.97.RData")

################################################################################
## ALPHA DIVERSITY
## EFFECTS OF OXYGEN, TEMPERATURE AND THEIR INTERACTION
## ACCOUNTING FOR TANK
################################################################################

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

################################################################################
## PREPARE DATA
################################################################################

alpha_data <- plot_data %>%
  mutate(
    Tank = factor(Tank),
    Ox = factor(Ox),
    Temp = factor(Temp),
    Week = factor(Week)
  )

################################################################################
## CHECK THE DATA
################################################################################

table(alpha_data$Week, useNA = "ifany")

table(
  alpha_data$Week,
  alpha_data$Ox,
  alpha_data$Temp,
  useNA = "ifany"
)

table(alpha_data$Tank, useNA = "ifany")


################################################################################
## SET WEEK ORDER
################################################################################

alpha_data$Week <- factor(
  alpha_data$Week,
  levels = c("A", "R", "S")
)

################################################################################
## NOTE
##
## In your data:
##
## A = Week 1 / acclimation
## R = Week 3 / during oxygen treatment
## S = Week 4 / after oxygen treatment
##
## Check this with:
################################################################################

table(alpha_data$Week)
#   A   R   S 
# 143 130 136

################################################################################
## MODEL 1:
## CHANGE FROM WEEK 1 -> WEEK 3
################################################################################

alpha_13 <- alpha_data %>%
  filter(
    Week %in% c("A", "R")
  ) %>%
  droplevels()

alpha_13$Week <- factor(
  alpha_13$Week,
  levels = c("A", "R")
)

################################################################################
## LMM
################################################################################

model_13 <- lmer(
  diversity_shannon ~ Week * Ox * Temp + (1 | Tank),
  data = alpha_13,
  REML = TRUE
)

summary(model_13)

anova(model_13)


################################################################################
## ESTIMATED MEANS
## Useful for interpreting the Week x Ox x Temp interaction
################################################################################

emmeans(
  model_13,
  ~ Week * Ox * Temp
)

################################################################################
## TEST CHANGE BETWEEN WEEK 1 AND WEEK 3
## WITHIN EACH OXYGEN x TEMPERATURE COMBINATION
################################################################################

emmeans(
  model_13,
  pairwise ~ Week | Ox * Temp
)


################################################################################
## MODEL 2:
## CHANGE FROM WEEK 3 -> WEEK 4
################################################################################

alpha_34 <- alpha_data %>%
  filter(
    Week %in% c("R", "S")
  ) %>%
  droplevels()

alpha_34$Week <- factor(
  alpha_34$Week,
  levels = c("R", "S")
)

################################################################################
## LMM
################################################################################

model_34 <- lmer(
  diversity_shannon ~ Week * Ox * Temp + (1 | Tank),
  data = alpha_34,
  REML = TRUE
)

summary(model_34)

anova(model_34)


################################################################################
## ESTIMATED MEANS
################################################################################

emmeans(
  model_34,
  ~ Week * Ox * Temp
)

################################################################################
## TEST CHANGE BETWEEN WEEK 3 AND WEEK 4
## WITHIN EACH OXYGEN x TEMPERATURE COMBINATION
################################################################################

emmeans(
  model_34,
  pairwise ~ Week | Ox * Temp
)


################################################################################
## OTHER ALPHA-DIVERSITY MEASURES
################################################################################
##
## You can repeat the same models for other diversity measures.
##
## Example: inverse Simpson
################################################################################

model_13_invSimpson <- lmer(
  diversity_inverse_simpson ~ Week * Ox * Temp + (1 | Tank),
  data = alpha_13,
  REML = TRUE
)

summary(model_13_invSimpson)

anova(model_13_invSimpson)


model_34_invSimpson <- lmer(
  diversity_inverse_simpson ~ Week * Ox * Temp + (1 | Tank),
  data = alpha_34,
  REML = TRUE
)

summary(model_34_invSimpson)

anova(model_34_invSimpson)


################################################################################
## EXAMPLE: OBSERVED RICHNESS
################################################################################

model_13_observed <- lmer(
  observed ~ Week * Ox * Temp + (1 | Tank),
  data = alpha_13,
  REML = TRUE
)

summary(model_13_observed)

anova(model_13_observed)


model_34_observed <- lmer(
  observed ~ Week * Ox * Temp + (1 | Tank),
  data = alpha_34,
  REML = TRUE
)

summary(model_34_observed)

anova(model_34_observed)


###############################################################
## DISTANCE MATRICES
###############################################################

jaccard_dist = phyloseq::distance(ps.skin.rare, method = "bray", binary = T)
bray_dist    = phyloseq::distance(ps.skin.rare, method = "bray", weighted = T)
unif_dist    = phyloseq::distance(ps.skin.rare, method = "unifrac", weighted = F)
wunif_dist   = phyloseq::distance(ps.skin.rare, method = "unifrac", weighted = T)

meta_df <- as(sample_data(ps.skin.rare), "data.frame")
head(meta_df)


###############################################################
## (1a) ADDITIVE MODEL
###############################################################

## JACCARD
jacc.adonis.simple <- adonis2(
  #  jaccard_dist ~ type_1 + stream_1 + location,
  jaccard_dist ~ Week + Temp + Ox,
  data = meta_df,
  permutations = 9999,
  by = "margin"
)

## BRAY
set.seed(123)
bray.adonis.simple <- adonis2(
  bray_dist ~ Week + Temp + Ox,
  data = meta_df,
  permutations = 9999,
  by = "margin"
)

## UNIFRAC
unif.adonis.simple <- adonis2(
  #unif_dist ~ type_1 + stream_1 + location,
  unif_dist ~ Week + Temp + Ox,
  data = meta_df,
  permutations = 999,
  by = "margin"
)

## WUNIFRAC
wunif.adonis.simple <- adonis2(
  #  wunif_dist ~ type_1 + stream_1 + location,
  wunif_dist ~ Week + Temp + Ox,
  data = meta_df,
  permutations = 999,
  by = "margin"
)

## JACCARD
ordination.jacc = ordinate(ps.skin.rare, method = "PCoA", distance = jaccard_dist)

jaccard <- plot_ordination(ps.skin.rare, ordination.jacc, color = "Treatment_main") + 
  theme(aspect.ratio=1)

jaccard <- plot_ordination(
  ps.skin.rare,
  ordination.jacc,
  color = "Treatment_main"
) +
  facet_wrap(
    ~ Ox
  ) +
  theme(
    aspect.ratio = 1
  )

jaccard


################################################################################
## JACCARD PCoA WITH TREATMENT CENTROIDS
################################################################################

# extract coordinates from the ordination
ord_data <- plot_ordination(
  ps.skin.rare,
  ordination.jacc,
  justDF = TRUE
)

# calculate centroids for each Treatment_main
centroids <- ord_data %>%
  dplyr::group_by(
    Treatment_main,
    Ox,
    Temp
  ) %>%
  dplyr::summarise(
    Axis.1 = mean(Axis.1, na.rm = TRUE),
    Axis.2 = mean(Axis.2, na.rm = TRUE),
    .groups = "drop"
  )

################################################################################
## PLOT
################################################################################

jaccard <- plot_ordination(
  ps.skin.rare,
  ordination.jacc,
  color = "Treatment_main"
) +
  
  geom_point(
    data = centroids,
    aes(
      x = Axis.1,
      y = Axis.2
    ),
    shape = 4,
    size = 5,
    stroke = 1.5,
    colour = "black",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = centroids,
    aes(
      x = Axis.1,
      y = Axis.2,
      label = Treatment_main
    ),
    colour = "black",
    size = 4,
    fontface = "bold",
    vjust = -1,
    inherit.aes = FALSE
  ) +
  
  facet_grid(
    Ox ~ Temp
  ) +
  
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

jaccard



jaccard_custom <- jaccard + 
  geom_point(size = 2.5, alpha = 0.5, shape = 21,
             colour = "black", stroke = 0.4,
             aes(fill = Treatment_main)) +
  stat_ellipse(geom = "path", type = "norm",
               linetype = "dashed", linewidth = 0.6,
               aes(colour = Treatment_main)) +
#  scale_fill_manual(values = color_map) +
#  scale_colour_manual(values = color_map) +
  theme_minimal(base_size = 18) +
  theme(
    aspect.ratio = 1,
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_text(size = 18, vjust = -0.5),
    axis.title.y = element_text(size = 18, vjust = 1),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)
  ) +
  coord_cartesian(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.6)) +
  labs(x = "\nPCoA axis 1 (25.2%)", y = "PCoA axis 2 (9.6%)\n")


################################################################################
## JACCARD PCoA WITH TREATMENT CENTROIDS
## FACET BY TIME
################################################################################

# extract coordinates from the ordination
ord_data <- plot_ordination(
  ps.skin.rare,
  ordination.jacc,
  justDF = TRUE
)

# calculate centroids for each Treatment_main
centroids <- ord_data %>%
  dplyr::group_by(
    Treatment_main,
    Week
  ) %>%
  dplyr::summarise(
    Axis.1 = mean(Axis.1, na.rm = TRUE),
    Axis.2 = mean(Axis.2, na.rm = TRUE),
    .groups = "drop"
  )

################################################################################
## PLOT
################################################################################

jaccard <- plot_ordination(
  ps.skin.rare,
  ordination.jacc,
  color = "Treatment_main"
) +
  
  geom_point(
    data = centroids,
    aes(
      x = Axis.1,
      y = Axis.2
    ),
    shape = 4,
    size = 5,
    stroke = 1.5,
    colour = "black",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    data = centroids,
    aes(
      x = Axis.1,
      y = Axis.2,
      label = Treatment_main
    ),
    colour = "black",
    size = 4,
    fontface = "bold",
    vjust = -1,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(
    ~ Week
  ) +
  
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

jaccard
