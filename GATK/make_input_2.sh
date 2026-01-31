#!/bin/bash
#SBATCH --job-name=make_input
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=4:00:00
#SBATCH --mem=250gb
#SBATCH --output=make_input2.out
#SBATCH --error=make_input2.err

#module load gcc/11.3.0 plink/1.07
module load gcc/11.3.0 plink/1.9

#edit pheno file
awk 'BEGIN{OFS="\t"} NR==1{print "FID","IID","Sex"; next} {print $1,$2,$3}' lepidum_old.pheno > lepidum.Sex.pheno
#make .bed .bim .fam files
plink --file fragi --pheno fragi.Sex.pheno --pheno-name Sex --make-bed --out fragi

#check the output files: awk '{print $2, $6}' ot.fam
#male = 2 and female = 1
