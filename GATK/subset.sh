#!/bin/bash
#SBATCH --job-name=subset
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=5:00:00
#SBATCH --mem=320gb
#SBATCH --output=subset.out
#SBATCH --error=subset.err
#SBATCH --array=21-31

#   Load modules and dependencies
module load gcc/12.1.0 gatk

#   Define paths for reference genome, genotyped data, and where our filtered data will go
REF="/scratch/group/moranlab/shared/darters/refs/Espec_sym_ref_NCBI/GCF_008692095.1_UIUC_Espe_1.0_genomic.fa"
RAW="/scratch/group/moranlab/shared/darters/RAW_WGS/Oligocephalus_gvcf/genotyped" 
FILTERED="/scratch/group/moranlab/shared/darters/RAW_WGS/Oligocephalus_gvcf/genotyped/filtered" 


#   We define the output region based on the basename of the intervals file that specifies the chromosome or group of unscaffolded reads
#   This was read from the file-of-files
REGFOF="/scratch/group/moranlab/shared/darters/Espectabile_Genome_Regions/GATK_Regions.fof"
REGION=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${REGFOF})
REG_BNAME=$(basename ${REGION})
REG_OUT=${REG_BNAME/.intervals/}

##################################################
###### Subset NO_VARIATION, SNPS, MIXED/INDELS #######
##################################################
#   T = Program to run
#   R = Reference genome
#   V = input
#   o = output
#   [selectTypeToExclude = what to get rid of] or [selectType = what to keep]
#   nt = number of threads

#   Subset monomorphic sites
gatk SelectVariants \
     -R ${REF} \
     -V ${RAW}/${REG_OUT}.vcf \
     --selectTypeToExclude INDEL \
     --selectTypeToExclude MIXED \
     --selectTypeToExclude MNP \
     --selectTypeToExclude SYMBOLIC \
     --selectTypeToInclude NO_VARIATION \
     -o ${FILTERED}/subset_mono/${REG_OUT}_mono_subset.vcf

echo -n "Done: subset mono"
date


#   Subset SNPs
java -Djava.io.tmpdir=/scratch/user/dkimfish -jar ${GATK} \
    -T SelectVariants \
    -R ${REF} \
    -V ${RAW}/${REG_OUT}.vcf \
    -selectType SNP \
    -nt 32 \
    -o ${FILTERED}/subset_snps/${REG_OUT}_snps_subset.vcf

echo -n "Done: subset snps"
date


#   Subset mixed/indels
java -Djava.io.tmpdir=/scratch/user/dkimfish -jar ${GATK} \
    -T SelectVariants \
    -R ${REF} \
    -V ${RAW}/${REG_OUT}.vcf \
    -selectType MIXED \
    -selectType INDEL \
    -nt 32 \
    -o ${FILTERED}/subset_mixed/${REG_OUT}_mixed_indels_subset.vcf

echo -n "Done: subset mixed/indels"
date
