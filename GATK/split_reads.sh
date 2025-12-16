#!/bin/bash
#SBATCH --job-name=split_reads_3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=64gb
#SBATCH --output=split_reads_3.out
#SBATCH --error=split_reads_3.err
#SBATCH --array=0-55

module load gcc/12.1.0 samtools/1.18

#set directory where the bam files are
export BAM_DIR="/home/wrad07/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125/Aligned_Reads/dups_marked/"
#set directories where the reads will go
export MAP_DIR="/home/wrad07/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125/Aligned_Reads/dups_marked/mapped"
export UNMAP_DIR="/home/wrad07/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125/Aligned_Reads/dups_marked/unmapped"

#get a list of the input alignments
ALNS=($(find ${BAM_DIR} -maxdepth 1 -name '*_dupsmarked.bam' | sort -V))

#use array to run the commands below
C_ALN=${ALNS[${SLURM_ARRAY_TASK_ID}]}

#get the sample name from the file name. We will use the sample name for the
#sample flags
bname=$(basename ${C_ALN})
sample=${bname/.bam}

#run sam tools to split the reads into mapped and unmapped
samtools view -hb -@ 16 -f 4 ${C_ALN} > ${UNMAP_DIR}/${sample}_unmapped.bam
samtools view -hb -@ 16 -F 260 ${C_ALN} > ${MAP_DIR}/${sample}_mapped.bam
