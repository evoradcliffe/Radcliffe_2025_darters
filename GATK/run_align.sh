#!/bin/bash
#SBATCH --job-name=align_reads
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=36:00:00
#SBATCH --mem=62gb
#SBATCH --output=align_reads.out
#SBATCH --error=align_reads.err
#SBATCH --array=1-105

#Load modules
module load GCCcore/12.3.0
module load parallel/20230722

#Specify the job command list
CMD_LIST="raw_align_commands.txt"

#This is standard command to run an array!
CMD="$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${CMD_LIST})"
eval ${CMD}
