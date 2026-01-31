#!/bin/bash
#SBATCH --job-name=make_input
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=4:00:00
#SBATCH --mem=250gb
#SBATCH --output=make_input1_allo.out
#SBATCH --error=make_input1_allo.err

#load modules and dependencies
module load gcc/12.1.0 vcftools/0.1.16

#make .ped .map files
## if file is gzipped or bgzipped:
vcftools --gzvcf ./renamed_elepidum_biallelic.vcf.gz --out lepidum --plink
