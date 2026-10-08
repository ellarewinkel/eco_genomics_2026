# objective: show how the differences in AM and OWA change across generations (F0, F4, F11)
# method: in DESeq2, subset by generation and pull out specific contrasts, compare results
# after, compare results: number of DE genes (possibly using Venn diagram for visualization)
# answer this Q: Does GO functional enrichment changes?

setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")


## Import the libraries that we're likely to need in this session

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)


####################################################

### Import our data

####################################################


# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)


####################################################################
# now filter out the OA/OW groups HERE and check if F11 subset works (before DESeq!!)
# script similar to subset function below (select = treatment = _)

####################################################################
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment) # normalizes data

dim(dds)
# [1] 67916     38

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)
# [1] 25260    12

dds_F4 <- subset(dds, select = generation == 'F4')
dim(dds_F4)
# [1] 25260    12

dds_F11 <- subset(dds, select = generation == 'F11')
dim(dds_F11)
# [1] 25260     6

# Perform DESeq2 analysis on the subset
# This looks at the effect of generation, NOT treatment
# because F0 is the first generation so there is no effect of treatment
dds_F0 <- DESeq(dds_F0)

dds_F4 <- DESeq(dds_F4)

dds_F11 <- DESeq(dds_F11)
# this one won't run (ERROR) right now; the other generations have 12 variables to compare
# but this one only has 6 (3 AM; 6 OWA)


####################################################

### Check on the DE results from the DESeq 

####################################################

resultsNames(dds_F0)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM <- res_OWAvsAM[order(res_OWAvsAM$padj),] # sort by significance
head(res_OWAvsAM)  
summary(res_OWAvsAM)
# interpreted from summary:
# 9.3% (or 2343) of 25260 genes running were UPREGULATED (responding to new environment)
# baseMean shows average counts per transcript
# log2FoldChange (LFC) shows pos or neg value, showing up/down-regulation
# lfcSE = LFC standard error
# stat = significance (the sign of it doesn't matter; just indicates LFC)
# p value = significance
# padj = adjusted p value using correction factor (slightly less sig)