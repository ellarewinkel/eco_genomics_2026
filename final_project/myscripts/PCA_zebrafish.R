
setwd("~/projects/eco_genomics_2026/final_project/mydata")

install.packages("readr")
library(readr)

# Reading TSV file using read_tsv()
df <- read_tsv('/gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/final_project/mydata/Qiimeoutput/otu_table.txt')
print(df)

