# this is in the vacc/terminal but i'm putting my script here for later reference.
# SCROLL DOWN FOR LOOKING AT DATA IN R!



# SETUP 2.1 Files for today

##change to your directory
cd /gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/mydata
  
##lets call it DIR for ease
DIR=/gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/mydata
    
###copy files across, at this stage should just be pop_level_biooracle.csv
cp /gpfs1/cl/biol3990/PopulationGenomics/Selection/pop_level_biooracle.csv $DIR
  
#check you have the data needed from last week
ll

# sample IDs from the VCF header line (#CHROM ...), one per line
zgrep -m 1 "^#CHROM" Lottia_gigantea_flt3_sub0.01.vcf.gz | cut -f 10- | tr '\t' '\n' > vcf_samples.txt
# it says broken pipe because you're only asking for 10, not to unzip the whole thing, so it only unzipped part of it

wc -l vcf_samples.txt #should be 553
head vcf_samples.txt # is it as expected? it's ASU01 through ASU10
  
  
  
#  3.2 Genome-wide FST between two populations

##first make sure you are in the right directory
## should be </users/d/s/dsadler1/projects/eco_genomics_2026/PopulationGenomics/mydata>
cd $DIR

###load vcftools
module load gcc
module load vcftools
###remember to note versions
###lets do the usual assigning of VCF:
VCF=$DIR/Lottia_gigantea_flt3_sub0.01.vcf.gz
POP1=ASU     # This can be any population, but lets do most southern ...
POP2=KRU     # ... and most northern.
# we would expect F_ST to be higher between populations because of isolation by distance
# the farther apart that populations are, the higher the F_ST factor due to lower breeding
# isolation by distance is the driving factor when other factors are neutral

mkdir -p fst #lets make a FST directory 

###Ok now lets run vcftools 
vcftools --gzvcf ${VCF} \
--weir-fst-pop pops/${POP1}.txt \
--weir-fst-pop pops/${POP2}.txt \
--out fst/${POP1}_vs_${POP2} \
2> fst/${POP1}_vs_${POP2}.log ##prints what would be displayed in the terminal to summarise mean FST

grep "Fst estimate" fst/${POP1}_vs_${POP2}.log
# grep is a search command to find that exact term (where there's an Fst estimate in that file)
# Weir and Cockerham mean Fst estimate: 0.016245
# Weir and Cockerham weighted Fst estimate: 0.020162

head fst/${POP1}_vs_${POP2}.weir.fst



# 3.3 Genome-wide FST between all pairs of populations
##
POPS=( $(cut -f2 popmap.txt | sort -u) ) ##f extracts column, 2 extracts second column, sort -u keeps only unique values (i.e., gives you one row per pop)
mkdir -p fst/pairwise
echo -e "pop1\tpop2\tmean_fst\tweighted_fst" > fst/pairwise_fst_summary.txt #makes empty txt file headed as pop1 etc

for (( i=0; i<${#POPS[@]}; i++ )); do ### starts at first row
  for (( j=i+1; j<${#POPS[@]}; j++ )); do ###starts at second row 
    p1=${POPS[$i]} ##i allocated to p1 (first pairwise)
    p2=${POPS[$j]} ##j allocated to p2 (second pairwise)
    vcftools --gzvcf ${VCF} \
    --weir-fst-pop pops/${p1}.txt \
    --weir-fst-pop pops/${p2}.txt \
    --out fst/pairwise/${p1}_vs_${p2} \
    2> fst/pairwise/${p1}_vs_${p2}.log ##same vcftools code as above but looped through
    mean=$(grep "mean Fst" fst/pairwise/${p1}_vs_${p2}.log | awk '{print $NF}') ##extracts the mean FST through each loop using grep
    weighted=$(grep "weighted Fst" fst/pairwise/${p1}_vs_${p2}.log | awk '{print $NF}') ##also the weighted
    echo -e "${p1}\t${p2}\t${mean}\t${weighted}" >> fst/pairwise_fst_summary.txt ###adds each loop result to the txt file
    done
    done ##note two done's as we have a for loop within a for loop!
    
    
    column -t fst/pairwise_fst_summary.txt
    
##########################################  
# NOW! we can look at the results in R
    
    
    library(tidyverse)
    
    setwd("/gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/mydata")
    
    ###point to your results folder
    OUTDIR="/gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/myresults"
    
    env <- read_csv("sample_metadata.csv")
    pw  <- read_tsv("fst/pairwise_fst_summary.txt")
    
    # order the populations from south to north
    env <- env %>%
      distinct(Population, Latitude, Longitude)
    pop_order <- env$Population[order(env$Latitude)]
    
    # add each pair in both directions so the heatmap is symmetrical
    pw_both <- bind_rows(pw, select(pw, pop1 = pop2, pop2 = pop1, mean_fst, weighted_fst)) %>%
      mutate(pop1 = factor(pop1, levels = pop_order),
             pop2 = factor(pop2, levels = pop_order))
    
    ggplot(pw_both, aes(x = pop1, y = pop2, fill = weighted_fst)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = round(weighted_fst, 3)), size = 2.5) +
      scale_fill_gradient(low = "white", high = "red4", name = "Weighted\nFST") +
      labs(x = NULL, y = NULL,
           title = "Genome-wide pairwise FST (populations ordered south to north)") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0(OUTDIR,"/pairwise_fst_heatmap.png"), width = 7, height = 6)
    
  
  # Does differentiation increase with the geographic distance between populations?
    pw_dist <- pw %>%
      left_join(select(env, Population, lat1 = Latitude), by = c("pop1" = "Population")) %>%
      left_join(select(env, Population, lat2 = Latitude), by = c("pop2" = "Population")) %>%
      mutate(lat_diff = abs(lat1 - lat2))
    
    ggplot(pw_dist, aes(x = lat_diff, y = weighted_fst)) +
      geom_point(size = 2) +
      geom_smooth(method = "lm", se = FALSE, colour = "grey40") +
      labs(x = "Difference in latitude (degrees)", y = "Weighted FST",
           title = "Isolation by distance?") +
      theme_classic()
    
    ggsave(paste0(OUTDIR,"/fst_vs_distance.png"), width = 5, height = 4)
  


##################################
#GO BACK TO VACC/TERMINAL (NOT R!)

# 3.4: Scanning the genome in windows

###lets use the example two populations again.
##remember you assigned VCF and POP1/2 earlier, check above if needed
    
vcftools --gzvcf ${VCF} \
  --weir-fst-pop pops/${POP1}.txt \
  --weir-fst-pop pops/${POP2}.txt \
  --fst-window-size 10000 \
  --fst-window-step 10000 \
  --out fst/${POP1}_vs_${POP2}_FST_10kb \
  2> fst/${POP1}_vs_${POP2}_FST_10kb.log

cd fst

ll
    
head ASU_vs_KRU_FST_10kb.windowed.weir.fst

############################################################### WE STOPPED HERE!

# 3.5: Plotting the scan and finding outlier windows
##in case you closed R
OUTDIR="/gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics/mydata"

#otherwise:
cd .. #to get to mydata

pop1 <- "ASU"   #the same two populations as in the bash code
pop2 <- "KRU"

fst_win <- read_tsv(paste0("fst/", pop1, "_vs_", pop2, "_10k.windowed.weir.fst"))