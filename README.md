# Repository for code used in Radcliffe et al. (2025). 
The following is code used for the GATK pipeline (from raw reads to biallelic vcf format for the GWAS for sex), MAF and MT/nuc analyses in Radcliffe et al. (2025) *Sex chromosome turnover and mitonuclear conflict drive reproductive isolation*.
########################################################################################################

*GATK PIPELINE*

Step 1: Processing raw WGS reads. 

1a) Trimmomatic and 1b) CutAdapt remove Illumina adapters.

1c) Combine a and c for paired end reads. > files should end in _R1.fastq.gz and _R2.fastq.gz

*script(s)*: 

a. use_trimmomatic.sh

b. run_cutadapt.sh

c. combine_reads.sh

########################################################################################################

Step 2: Align reads to *Etheostoma spectabile* reference genome.

2a) Make commands file that will run the alignment on multiple samples at once. 

2b) Run the alignment array.

files end in _raw.bam

*script(s)*:

a. ./make_align_jobs.sh > raw_align_commands.txt

b. run_array.sh (*standard command for running an array job, using the raw_align_commands.txt file generated in 1.)*

########################################################################################################

Step 3: Process aligned reads. 

3a) Add read groups > files end in _rgadd.bam

3b) Mark duplicates > files end in _dupsmarked.bam

3c) Split into mapped and unmapped reads > files end in either _mapped.bam or _unmapped.bam

script(s): *all these were run with array jobs, so follow format of make commands.txt file, then run job with the run_array.sh from above*.

a. ./make_align_rgadd.sh > rgadd_commands.txt

b. ./make_jobs_dups.sh > dups_commands.txt

c. split_reads.sh *does not need to be run into an array, will simply work on every file in given directory*

########################################################################################################

Step 4: HaplotypeCaller - calling variants (SNPs and indels).

files will end in .g.vcf.gz

script: haplotype_caller.sh

########################################################################################################

Step 4: Genotype GVCFs. *Sorts variants into chromosome-level vcf (variant call format) file*

4a) Must combine all .g.vcf files into a single "cohort" .g.vcf for genotyping.

4b) Genotyping.

files will end in .vcf

*script(s)*: 

a. combine_gvcfs.sh

b. genotype_gvcf.sh

########################################################################################################

Step 5: VCF Processing 

5a) Subset 

5b) Hard-filter for variants that pass filters.

5c) Biallelic filtering

*script(s)*:

a. subset.sh

b. filter_gatk4.5.sh

c. biallelic_maf.sh

########################################################################################################

Step 6: Sex association and population genetics analysis (FST, SNP Density, Heterozygosity)

########################################################################################################

*Plotting scripts for step 6 are in the "plots" folder.*

########################################################################################################

*MT-Nuc Analyses*
