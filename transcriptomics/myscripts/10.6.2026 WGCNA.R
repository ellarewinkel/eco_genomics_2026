#############################################################################
# WGCNA Analyses
#############################################################################
setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")

# now in terminal: cd /gpfs1/cl/biol3990/Transcriptomics/WGCNA
# now: ll
# now: cp * ~/projects/eco_genomics_2026/transcriptomics/mydata
# the files Ahud_trait_data.txt and bwnet.rds should now be in mydata

# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

library(DESeq2)
library(ggplot2)

library(tidyverse)

library(CorLevelPlot) 
library(gridExtra)

library(Rmisc) 



# 1. Import the counts matrix and metadata and filter using DESeq2

countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1) # we give this intercept of 1 this time bc it's a blind design.

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(countsTable)) #gsg is good sample genes in DNA; pulls out outliers that would've skewed data
summary(gsg)
gsg$allOK # FALSE: not all outliers were removed (otherwise would've said TRUE)

table(gsg$goodGenes) # FALSE: 3874; TRUE: 64042
table(gsg$goodSamples) # TRUE: 38


# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) # shows them clustering differently



# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
# point is that you can't see any patterns when you're taking all the genes w/o organization

# All of our samples look good; the genes that don't look good shouldn't cluster well and 
# should end up in the "grey", catch all module for uncorrelated genes.

# 3. Normalization ----------------------------------------------------------------------

colData <- row.names(sample_metadata) # WGCNA wants all column names to be the row names

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
   # "TRUE" = yes, all samples are present in both columns/rows
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order
   # "TRUE" = yes, all samples are in same order



# perform variance stabilization
dds_norm <- vst(dds)
# dds_norm <- vst(normalized_counts)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices 

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 6 # should try running with 7 based on R^2
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package
