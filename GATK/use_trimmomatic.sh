#!/bin/bash
#SBATCH --job-name=use_trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=72:00:00
#SBATCH --mem=62gb
#SBATCH --output=trimmomatic.out
#SBATCH --error=trimmomatic.err

module load openjdk/17.0.2 trimmomatic/0.39
module load gcc/12.1.0 perl/5.36.0 parallel/20240622

export HOME=~/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125

paralleltrim() {
  #Get each sample name from the list SampleID.txt (supplied to the parallel command with the -a flag)
  SAMP=${1}
  
  #trim low quality
  java -Xmx2g -jar /home/wrad07/rmoran-lab/projects/gonad_RNAseq_2025/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 ${HOME}/${SAMP}_R1_001.fastq.gz ${HOME}/${SAMP}_R2_001.fastq.gz ${HOME}/Trimmomatic_Out/${SAMP}_trim_pair_R1.fastq.gz ${HOME}/Trimmomatic_Out/${SAMP}_trim_unpair_R1.fastq.gz ${HOME}/Trimmomatic_Out/${SAMP}_trim_pair_R2.fastq.gz ${HOME}/Trimmomatic_Out/${SAMP}_trim_unpair_R2.fastq.gz SLIDINGWINDOW:6:30 MINLEN:40
}

#   Export function so we can call it with parallel
export -f paralleltrim

#NAvigate into the reads directory
cd ${HOME}

parallel --joblog Trimmomatic_parallel_logfile.txt -a Sample_IDs.txt paralleltrim
