#!/bin/bash
#SBATCH --job-name=sex_association
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=124G
#SBATCH --time=24:00:00

module load gcc/12.1.0 gemma/0.98.5

#in_dir="/home/wrad07/rmoran-lab/projects/darter_WGS/gvcfs_2025/sdr"
#output_dir="${in_dir}/output"

#gemma -bfile ${in_dir}/lepidum_sexed -gk 2 -o lepidum_sexed
#gemma -bfile ${in_dir}/lepidum -k ./output/lepidum.sXX.txt -lmm 4 -o lepidum

#gemma -bfile ${in_dir}/fragi_thin -gk 2 -o fragi_thin
#gemma -bfile ${in_dir}/fragi_thin -k ./output/fragi_thin.sXX.txt -lmm 4 -o fragi_thin

gemma -bfile ./fragi_input/fragi -k /home/wrad07/rmoran-lab/projects/darter_WGS/gvcfs_2025/sdr/fragi_input/output/fragi_all.sXX.txt -lmm 4 -o fragi_all
