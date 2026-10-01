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


dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment) # normalizes data

dim(dds)
# [1] 130580     38

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)
# [1] 25260    12

# Perform DESeq2 analysis on the subset
# This looks at the effect of generation, NOT treatment
# because it's the first generation so there is no effect of treatment
dds_F0 <- DESeq(dds_F0)




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


res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)
res_OWvsAM <- res_OWvsAM[order(res_OWvsAM$padj),]
head(res_OWvsAM) 
summary(res_OWvsAM)

# And one more... ?!
res_OAvsAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAsAM <- res_OAvsAM[order(res_OAvsAM$padj),]
head(res_OAvsAM) 
summary(res_OAvsAM)

# results show that OA has less of an inmpact



### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)

# make ggplot
# we keep redefining p because _
d <-plotCounts(dds_F0, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment")), returnData=TRUE)
d

# you can plot other genes if you switch the gene name above, but we didn't do any more here.

p <-ggplot(d, aes(x=treatment, y=count, color=treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p




# Make an MA plot
# x-axis shows mean value of normalized counts of that gene in the data set
# y-axis shows log fold change from ambient, which is the thick grey line at y=0.
# blue color shows significant genes; grey shows nonsignificant genes
plotMA(res_OWvsAM, ylim=c(-5,5))



# make volcano plot
volcano_df <- as.data.frame(res_OWvsAM) # changes into a data frame object

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "steelblue",
    "NS"   = "grey70",
    "Up"   = "firebrick"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )


# MA and volcano plots are different by: 
# _ helps with interpretation and visualization



########################################

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

# By environment
vsd <- vst(dds_F0, blind=FALSE)
# vst command groups genes that are expressed similarly (not necessarily related)

topgenes <- head(rownames(res_OWvsAM),100)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_F0)[,c("treatment", "generation")])
pheatmap(mat, annotation_col=df)
pheatmap(mat, annotation_col=df, cluster_cols = F)

# Can you read that? Try this... How/why is it better?
pheatmap(mat, annotation_col=df, cluster_cols = F, show_rownames = F)
# shows that OWA group has way more genes expressed in 


# Make a Euler plot, which is better than a Venn diagram because it scales the size of the circles to sizes representative of data

#################################################################

#### PLOT OVERLAPPING DEGS IN VENN EULER DIAGRAM

#################################################################

# For OW vs AM
res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05) # pull out the results for the contrast of interest
res_OWvsAM <- res_OWvsAM[order(res_OWvsAM$padj),] # order them by significance
res_OWvsAM <- res_OWvsAM[!is.na(res_OWvsAM$padj),] # get rid of any NAs
degs_OWvsAM <- row.names(res_OWvsAM[res_OWvsAM$padj < 0.05,]) # make a list of significant differentially expressed genes for this contrast

length(degs_OWvsAM) # 5517 significant differentially expressed genes for OW vs AM

# For OA vs AM
res_OAvsAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM <- res_OAvsAM[order(res_OAvsAM$padj),]
res_OAvsAM <- res_OAvsAM[!is.na(res_OAvsAM$padj),]
degs_OAvsAM <- row.names(res_OAvsAM[res_OAvsAM$padj < 0.05,])

length(degs_OAvsAM) # 602 sig diff expressed genes for OA vs AM

# For OWA vs AM
res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM <- res_OWAvsAM[order(res_OWAvsAM$padj),]
res_OWAvsAM <- res_OWAvsAM[!is.na(res_OWAvsAM$padj),]
degs_OWAvsAM <- row.names(res_OWAvsAM[res_OWAvsAM$padj < 0.05,])

length(degs_OWAvsAM) # 3918 sig diff expressed genes for OWA vs AM

library(eulerr)

# Total
length(degs_OAvsAM)  # 602
length(degs_OWvsAM)  # 5517 
length(degs_OWAvsAM)  # 3918

# Intersections
length(intersect(degs_OAvsAM,degs_OWvsAM))  # 444
length(intersect(degs_OAvsAM,degs_OWAvsAM))  # 380
length(intersect(degs_OWAvsAM,degs_OWvsAM))  # 2743

# Shared across all
intWA <- intersect(degs_OAvsAM,degs_OWvsAM)
length(intersect(degs_OWAvsAM,intWA)) # 338

# Number unique to each treatment

602-444-380+338 # 116 OA
5517-444-2743+338 # 2668 OW 
3918-380-2743+338 # 1133 OWA

# Number shared in pairs of treatments

444-338 # 106 OA & OW
380-338 # 42 OA & OWA
2743-338 # 2405 OWA & OW

# Now assemble the results
# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("OA" = 116, "OW" = 2668, "OWA" = 1133, "OA&OW" = 106, "OA&OWA" = 42, "OW&OWA" = 2405, "OA&OW&OWA" = 338))

# And make the plot!
plot(fit1,  lty = 1:3, quantities = TRUE)
# lty changes the lines

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))


#cross check with above lengths of DEGS: the four values, unique, shared with one other, shared with the second other, shared across all treatments, should sum to the length of DEGs for each treatment contrast to AM
2668+2405+338+106 # 5517 total OW
1133+2405+338+42  # 3918 total OWA
116+42+106+338    # 602  total OA






# An Upset plot
# Here’s an upset plot… a bit easier
# Note: A new bit of code below! %in% This asks, “is this member of that group?”
install.packages("UpSetR")
library(UpSetR)

all_genes <- unique(c(
  degs_OAvsAM,
  degs_OWvsAM,
  degs_OWAvsAM
))

upset_df <- data.frame(
  gene = all_genes,
  OA = all_genes %in% degs_OAvsAM,
  OW = all_genes %in% degs_OWvsAM,
  OWA = all_genes %in% degs_OWAvsAM
)

head(upset_df)

deg.list <- list(
  OA  = degs_OAvsAM,
  OW  = degs_OWvsAM,
  OWA = degs_OWAvsAM
)

upset(
  fromList(deg.list),
  order.by = "freq",
  mainbar.y.label = "Number of DEGs",
  sets.x.label = "Total DEGs"
)

####################### A bit prettier data.
upset(
  fromList(deg.list),
  order.by = "freq",
  main.bar.color = "grey30",
  sets.bar.color = c("#00A08A", "#CC3333", "#F2AD00"), # had to manually adjust the order
  mainbar.y.label = "Number of DEGs",
  sets.x.label = "Total DEGs"
)





# Scatter plots to compare expression responses
# Building from your script from last class focused on just the F0 generation, we will pull out the results from the OW vs AM contrast and compare to the OWA vs AM contrast.
# Make sure you load your libraries, import the data, filter the data, run the DESeq model, and define your results dataframes before proceeding below…
# merge and tidyverse’s mutate are the core functions for doing what we are about to do -> create a new data frame of the data we want to plot and create a new column of data based on the values in the data frame!

#################################################################

#### Scatter plot to assess how correlated are responses to OWA vs OW?

#################################################################


# Create merged data frame - need to use rownames because differences in filtering
plot_OWA <- data.frame(
  gene = rownames(res_OWAvsAM),
  LFC_OWA = res_OWAvsAM$log2FoldChange,
  padj_OWA = res_OWAvsAM$padj
)
# this pulls out LFC and adjusted p-value


# now do the same for OW
plot_OW <- data.frame(
  gene = rownames(res_OWvsAM),
  LFC_OW = res_OWvsAM$log2FoldChange,
  padj_OW = res_OWvsAM$padj
)


# merge these two dataframes with merge function. merges them by gene (sorts/syncs)
plot_df <- merge(plot_OWA,
                 plot_OW,
                 by = "gene")

# Remove genes with missing LFC values
plot_df <- plot_df %>%
  filter(!is.na(LFC_OWA),
         !is.na(LFC_OW))

# Classify significance
# tidyverse function %>% tells you a sequence of events: take the dataframe and mutate it
# mutate function creates a new variable
plot_df <- plot_df %>%
  mutate(
    SigGroup = case_when(
      padj_OWA < 0.05 & padj_OW < 0.05 ~ "Both",
      padj_OWA < 0.05 ~ "OWA only",
      padj_OW < 0.05 ~ "OW only",
      TRUE ~ "Neither"
    )
  )

# Correlation for noting on the plot 
r <- cor(plot_df$LFC_OWA,
         plot_df$LFC_OW,
         use = "complete.obs")

# Arrange the genes by significant to make the plotting easier/more interesting to see
# ggplot plots in the order of the df, so random
# rearranging the order of the groups (ex: switching where it says Neither and Both) shows the layers of data differently
plot_df$SigGroup <- factor(
  plot_df$SigGroup,
  levels = c("Neither", "OWA only", "OW only", "Both")
)

plot_df <- plot_df %>%
  arrange(SigGroup)

# Now make the plot!
# aes within plot means "aesthetics"
# alpha means transparency/opacity (0 = fully transparent; 1 = fully opaque)

ggplot(plot_df,
       aes(x = LFC_OW,
           y = LFC_OWA,
           color = SigGroup)) +
  
  geom_point(alpha = 0.6, size = 1.5) +
  
  # the dashed line, slope of 1 running through
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              color = "black") +
  
  # hline is horizontal
  geom_hline(yintercept = 0,
             color = "grey70") +
  
  # vline is vertical
  geom_vline(xintercept = 0,
             color = "grey70") +
  
  # the r in here comes from the correlation above (at #370) that we saved as r
  # we're telling it to put the correlation at xmin and ymax, then round to the 3rd decimal place
  annotate("text",
           x = min(plot_df$LFC_OW, na.rm = TRUE),
           y = max(plot_df$LFC_OWA, na.rm = TRUE),
           hjust = 0,
           label = paste0("r = ", round(r, 3))) +
  
  scale_color_manual(values = c(
    "Both" = "purple",
    "OWA only" = "#CC3333",
    "OW only" = "#00A08A",
    "Neither" = "grey80"
  )) +
  
  coord_fixed() + # forces the same scaling on x and y axes (so the slope is properly represented)
  
  labs(
    x = "Log2 Fold Change: OW vs AM",
    y = "Log2 Fold Change: OWA vs AM",
    color = "",
    title = "GE Responses to OW relative to OWA"
  ) +
  
  # base_size is making the text size 14
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )



# Let’s test for functional enrichment using annotated GO categories for each gene and the TopGO program.
# First, you’ll need two files transcript_universe.csv and trinotate_annotation_GOblastx_forTopGO.txt that you can find in our class directory, /gpfs1/cl/biol3990/Transcriptomics/GOenrichment, and cp over to your mydata directory
# The first step is to create the saved results files with the correct trinity ids

