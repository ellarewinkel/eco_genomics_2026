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
