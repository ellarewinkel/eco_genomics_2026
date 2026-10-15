##Lets first create a new folder for our new module

cp /gpfs1/cl/biol3990/PopulationGenomics/VCF/*.vcf* /gpfs1/home/e/r/erewinke/projects/eco_genomics_2026/population_genomics

### go to your ecol genomics folder something like: "~/projects/eco_genomics_2026/" and make your new directory
mkdir -p population_genomics 

##copy the VCF files to your directory (in vacc shell)

cd projects/eco_genomics_2026/population_genomics/mydata/
  
zcat Lottia_gigantea_biallelic.vcf.gz | head -n 20


module load gcc
module load bcftools

bcftools --version   # record the versions you use in your methods
# bcftools 1.19
# Using htslib 1.19.1

bcftools -h          # -h or --help usually gets you the manual
# note: - is usually a shortened version of -- e.g. h vs help but both do the same thing.

bcftools head Lottia_gigantea_biallelic.vcf.gz

# Again we get a big string of contig lines. Let’s pipe (|) the output into grep to remove them:
bcftools head Lottia_gigantea_biallelic.vcf.gz | grep -v "##contig="
# grep (“global regular expression print”) searches text for a pattern and prints the lines that match.
# The -v option inverts this, printing only the lines that don’t match

# Now you can see how the file was made and what each field means. Each ##INFO and ##FORMAT line defines one of the tags used in the body of the file:
bcftools view -h Lottia_gigantea_biallelic.vcf.gz | grep "^##INFO"
bcftools view -h Lottia_gigantea_biallelic.vcf.gz | grep "^##FORMAT"

# bcftools also writes every command it runs into the header (in lines containing Command=), so you can often reconstruct how a file was made and filtered:
bcftools view -h Lottia_gigantea_biallelic.vcf.gz | grep "Command=":
  
# The VCF also stores our sample names, which we can list with bcftools query -l:
bcftools query -l Lottia_gigantea_biallelic.vcf.gz

# count the samples by piping into wc -l (counts lines)
bcftools query -l Lottia_gigantea_biallelic.vcf.gz | wc -l
# 553 samples

# Let's look at the first few SNPs to see what our data looks like
# The -H option drops the header so we only see variant records
# just the first record
bcftools view -H Lottia_gigantea_biallelic.vcf.gz | head -n 1

# the first 5 records; less -S stops long lines wrapping (press q to exit)
bcftools view -H Lottia_gigantea_biallelic.vcf.gz | head -n 5 | less -S

# -v selects variant types: here SNPs only
bcftools view -H -v snps Lottia_gigantea_biallelic.vcf.gz | head -n 5 | less -S

# what do you notice when we ask for indels?
bcftools view -H -v indels Lottia_gigantea_biallelic.vcf.gz | head -n 5 | less -S

# if the shell says (END) after this, type q and enter to get out of it.