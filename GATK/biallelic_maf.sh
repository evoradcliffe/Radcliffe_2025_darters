#!/bin/bash
#SBATCH --job-name=biallelic_maf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=2:00:00
#SBATCH --mem=16gb
#SBATCH --output=biallelic_maf.out
#SBATCH --error=biallelic_maf.err

module load GCC/11.2.0
module load BCFtools/1.14

bcftools view /scratch/group/moranlab/shared/darters/RAW_WGS/merged_ChrOnly_filtered_RepsRem_Indels10bpRem.vcf -M2 -e'MAF<0.01' -o /scratch/group/moranlab/shared/darters/RAW_WGS/merged_renamedChrOnly_filtered_RepsRem_Indels10bpRem_biallelic_snps.vcf
