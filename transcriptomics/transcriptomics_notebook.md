# Transcriptomics Notebook

**Course:** Intro to Ecological Genomics Fall 2026

**Name:** Ella Rewinkel

------------------------------------------------------------------------

## 9.15.2026 - Setting up lab notebook and learning markdown

-   Setting up transcriptomics notebook

-   Learn how to take notes in markdown

-   Push notes to github

**Working Directory:**

`~/projects/eco_genomics_2026/transcriptomics`

**Copy Paste this to set working directory:** setwd("\~/projects/eco_genomics_2026/transcriptomics/mydata")

**Input Files:**

`none`

**Output Files:**

`~/projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

-   `add them here`

**Scripts:**

`none`

**Code:**

``` r
print("Hello World")
```

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |

![](images/images.jpg)

**Notes / Observations:**

-   Oh cool graph! it makes sense! and some interpretation of your data

**Next Steps?**

-   Let's dive into the project!

------------------------------------------------------------------------

## 9.17.2026 - Transcriptomics Tutorial 1: Intro to the copepod study system & dataset

# Commands that I got to use on the VACC

``` r
zcat filename | head -n 4
```

This command opened a gzipped file. But don't run it without "piping" (\|) it to a head command.

``` r
wc
```

This command lets you view the word count. If you want the subcommand -l, say "wc -l" to get the number of lists.

<div>

#### Notes:

</div>

\`\`\`\`

##### Questions we could ask:

-   Multiplicative interaction: Additive response, Antagonistic, or Synergistic

    -   Graph

        -   X-axis: Expected LFC vs AM. Ocean Acidification (OA) vs Ocean Warming & Acidification (OWA)

        -   Y-axis: Observe OWA

        -   Additive: along y=x line

        -   Synergistic: above y=x

        -   Antagonistic: below y=x

-   Change in genetic expression (GE) in response to stressor, and which is more/most impactful across generations?

-   Change in GE thru time w/in a treatment –\> which are stable?

-   Number of genes, magnitude of DGE

-   Overlap across treatments

-   Plasticity vs Adaptive evolution

-   Phenotypes: which are correlated with which genes?

##### General workflow for analyzing gene exp data:

-   Clean raw seq data (completed)

    -   fastp (on github) on raw reads –\> cleaned reads

-   Generate/annotate de novo reference transcriptome assembly (completed)

    -   Use Trinity to evaluate for quality (length + completeness) using BUSCO

-   Map clean reads to ref assembly (completed)

    -   Use Salmon to ref transcriptome + quantify abundance

-   Test for differential expression among groups (START HERE)

    -   Import data into DESeq2 in R for data normalization, visualization, stat tests for differential GE

    -   "Model" definition: Gene X - What's the effect of OA + OW + OA:OW interaction? Look at values of expression within gene + whether it's characterized as a sample from a certain origin (OA, OW...). You're testrunning this for every transcript + must do false discovery rate correction

-   Perform more advanced analyses (do this too!)

    -   Use Weighted Gene Correlation Network Analysis (WGCNA) to identify gene clusters w/ correlated expression

    -   Use TopGO or GO Mann-Whitney U GO_MWU test for fxn'l enrichment (skill: fxn'l enrichment analysis) among genes differentially expressed btwn grps

##### What raw Illumina seq data looks like:

-   Seq data files have suffix .fastq (.fq)

-   With paired end sequencing: 2 files per sample

    -   Left read: \_1.fq.gz

    -   Right read: \_2.fq.gz

-   Each file has 4 lines for each read that include:

    -   Seq identifier (read name)

    -   Nucleotide sequence (A, T, G, C, sometimes N [bad])

    -   Separator line (usually just "+")

    -   Quality scores for each base

        -   Want Q/Phred Quality Score of at LEAST 30 (prob of incorrect base cell: 1/1000; base call accuracy: 99.9%)

        -   Letters in quality encoding mean good quality data; symbols are horrible and numbers aren't good either

        </div>

<div>

## 9.22.2026 - Day 3 of transcriptomics

Today we set up our R working env and copied the data to import into DESeq2. The class got to the end of saving the pretty ggplot, but I got to the end of making the MA plot.

### Data

PCA to visualize global gene expression patterns (ggplot)

![](myresults/PCA_allGens.png)

</div>

</div>

# 9.24.2026 Day 4 of Transcriptomics:

**Overview**

Today, we explored and overviewed basic R functions. We also reviewed pathways in R. I had issues with line 370 onward in ahud_DESeq2_inclass.R and plan to seek solutions next week.

**R Functions**

| Function Names      | Meaning/Purpose                          |
|---------------------|------------------------------------------|
| =                   | defines something                        |
| \<--                | use when creating \_                     |
| c(data, data, data) | c() encapsulates all data into one thing |

</div>

# 9.29.2026 Day 4 of Transcriptomics in Tutorials:

**Overview**

Today, we completed the Pespeni lab tutorial for Day 4 of Transcriptomics and analyzed data from plots we created.

**New Code**

``` r
%in%
```

Asks: "is this a member of that group?"

**Plots Created**

![GGplot for counts of specific top interaction gene](mydata/9.29.26%20GGplot%20for%20counts%20of%20specific%20top%20interaction%20gene.png){width="533"}

![Make an MA plot. The x-axis shows mean value of normalized counts of that gene in the data set. The y-axis shows log fold change from ambient, which is the thick grey line at y=0. \# blue color shows significant genes; grey shows insignificant genes](mydata/9.29.26%20MA%20plot.png){width="681"}

![Volcano plot helps with differentiation and visualization.](mydata/9.29.26%20Volcano%20plot.png)

![Heatmap of top 20 genes sorted by p-value.](mydata/9.29.26%20Heatmap%20of%2020%20genes%20sorted%20by%20pvalue.png)

![Euler plot, which is better than a Venn diagram because it scales the size of the circles to sizes representative of data. This type of plot can get complicated with more treatments.](mydata/9.29.26%20Euler%20plot.png)

![Upset plot](mydata/9.29.26%20Upset%20plot.png)

</div>

# 10.01.2026 Day 5 of Transcriptomics

**Overview**

Today, we used a scatterplot to compare expression response to OW relative to OWA (and each vs AM control). We filtered the data, annotated/classified genes, ordered results, and visualized with ggplot.

**Code**

We used four tidyverse (dplyr) functions:

-   filter() to remove rows

-   mutate() to add a new variable

-   case_when() to classify genes into categories

-   arrange() to sort the rows

**Plots**

![GGplot showing GE responses to OW relative to OWA. Looks at logfold change (LFC). The section for Both (purple), aligning with expectations, is mainly only upregulated or downregulated, following along the 1:1 line well, but not much near 0 (We filtered out the responses that were not significant.). For OW only (green), there is a magnitude of change spread out across x-axis for LFC between OW vs AM. For OWA only, it's vertical following y-axis which shows LFC for OWA vs AM. The "Neither" (gray) is plotted, despite filtering for significance, and shows massive LFC. "Neither" responses show up due to consistency across biological replicates (They have a lot of variation among replicates, so they are not significant due to not being consistently different.).](mydata/9.29.26 Day 5 ggplot.png)
