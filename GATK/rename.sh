#!/bin/bash
#SBATCH --job-name=rename
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00

module load gcc/12.1.0 bcftools/1.21

bcftools annotate --rename-chrs scaffold_to_chromosome_map.txt -Oz -o renamed_squamosum.vcf.gz esquamosum_biallelic_2025.vcf.gz

TMP="${SLURM_TMPDIR:-${SCRATCH:-$PWD}}/bcftools_sort_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMP"

bcftools sort -T "$TMP" -Oz -o renamed_squamosum.sorted.vcf.gz renamed_squamosum.vcf.gz
tabix -p vcf renamed_squamosum.sorted.vcf.gz
