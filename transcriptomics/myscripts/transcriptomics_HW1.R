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
library(topGO)


####################################################
# Import our data
####################################################
# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable)
# ^ bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)


####################################################################
# script similar to subset function below (select = treatment = _)
####################################################################
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~  treatment)

dim(dds)
# [1] 67916    38

# Filter out genes that lack enough reads (less than 15 counts in over 75% of samples, so ~28)
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]

nrow(dds) 
# [1] 25260

# Subset by generation before observing OWA vs AM results

#######
# F0 #
#######
dds_AM_OWA_F0 <- subset(dds, select = treatment %in% c("AM", "OWA") & generation == "F0")
dds_AM_OWA_F0$treatment <- droplevels(dds_AM_OWA_F0$treatment)
dds_AM_OWA_F0$treatment <- relevel(dds_AM_OWA_F0$treatment, ref = "AM")

table(dds_AM_OWA_F0$treatment)
# AM OWA
#  3   3

dds_AM_OWA_F0 <- DESeq(dds_AM_OWA_F0)
resultsNames(dds_AM_OWA_F0)
# [1] "Intercept"           "treatment_OWA_vs_AM"

res_OWAvsAMF0 <- results(dds_AM_OWA_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAMF0 <- res_OWAvsAMF0[order(res_OWAvsAMF0$padj),] # sort by significance
head(res_OWAvsAMF0)  
summary(res_OWAvsAMF0)
# notes: interpreted from summary of OWA vs AM (NOT subsetted by generation. must update this):
# 9.3% (or 2343) of 25260 genes running were UPREGULATED (responding to new environment)
# baseMean shows average counts per transcript
# log2FoldChange (LFC) shows pos or neg value, showing up/down-regulation
# lfcSE = LFC standard error
# stat = significance (the sign of it doesn't matter; just indicates LFC)
# p value = significance
# padj = adjusted p value using correction factor (slightly less sig)

write.csv(
  res_OWAvsAMF0.df,
  "myresults/OWAvsAMF0_results.csv",
  row.names = FALSE
)

#######
# F4 #
#######

dds_AM_OWA_F4 <- subset(dds, select = treatment %in% c("AM", "OWA") & generation == "F4")
dds_AM_OWA_F4$treatment <- droplevels(dds_AM_OWA_F4$treatment)
dds_AM_OWA_F4$treatment <- relevel(dds_AM_OWA_F4$treatment, ref = "AM")

table(dds_AM_OWA_F4$treatment)
# AM OWA
#  3   3

dds_AM_OWA_F4 <- DESeq(dds_AM_OWA_F4)
resultsNames(dds_AM_OWA_F4)
# [1] "Intercept"           "treatment_OWA_vs_AM"

res_OWAvsAMF4 <- results(dds_AM_OWA_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAMF4 <- res_OWAvsAMF4[order(res_OWAvsAMF4$padj),]
head(res_OWAvsAMF4)
summary(res_OWAvsAMF4)

#? until before F11
res_OWAvsAMF4.df <- as.data.frame(res_OWAvsAM)
res_OWAvsAM.df$fullID <- rownames(res_OWAvsAM.df)

parts <- strsplit(res_OWAvsAM.df$fullID, "::")

res_OWAvsAMF4.df$shortID <- sapply()
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAMF4.df,
  "myresults/OWAvsAMF4_results.csv",
  row.names = FALSE
)

#######
# F11 #
#######
dds_AM_OWA_F11 <- subset(dds, select = treatment %in% c("AM", "OWA") & generation == "F11")
dds_AM_OWA_F11$treatment <- droplevels(dds_AM_OWA_F11$treatment)
dds_AM_OWA_F11$treatment <- relevel(dds_AM_OWA_F11$treatment, ref = "AM")

table(dds_AM_OWA_F11$treatment)
#  AM OWA 
#   3   3

dds_AM_OWA_F11 <- DESeq(dds_AM_OWA_F11)
resultsNames(dds_AM_OWA_F11)
# [1] "Intercept"           "treatment_OWA_vs_AM"

res_OWAvsAMF11 <- results(dds_AM_OWA_F11, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAMF11 <- res_OWAvsAMF11[order(res_OWAvsAMF11$padj),]
head(res_OWAvsAMF11)
summary(res_OWAvsAMF11)

### GO Analysis on OWA vs AM Overall (not by generation)
library(topGO)

run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){
  
  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )
  
  deseq <- subset(
    deseq,
    !is.na(padj)
  )
  
  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )
  
  names(geneList) <- deseq$shortID
  
  cat("\nGenes tested:",
      length(geneList))
  
  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  resultWeight <- runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )
  
  GOresults <- GenTable(
    GOdata,
    weightFisher = resultWeight,
    orderBy = "weightFisher",
    topNodes = 100
  )
  
  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )
  
  return(GOresults)
}

GO_OWA <- run_topGO_contrast(
  "myresults/F0_OWAvsAM_results.csv",
  "myresults/GO_F0_OWAvsAM_BP.csv"
)
# Genes tested: 24259
# Significant genes: 3918


go <- read.csv(
  "myresults/GO_F0_OWAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  slice(1:10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )

### GO Contrasts
GO_OA <- run_topGO_contrast(
  "myresults/OAvsAM_results.csv",
  "myresults/GO_OAvsAM_BP.csv"
)

GO_OW <- run_topGO_contrast(
  "myresults/OWvsAM_results.csv",
  "myresults/GO_OWvsAM_BP.csv"
)

GO_OWA <- run_topGO_contrast(
  "myresults/OWAvsAM_results.csv",
  "myresults/GO_OWAvsAM_BP.csv"
)
