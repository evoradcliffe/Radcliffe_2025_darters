#!/bin/bash
#SBATCH --job-name=count_kmers_ec
#SBATCH --nodes=1
#SBATCH --partition=tier1q
#SBATCH --ntasks-per-node=12
#SBATCH --time=168:00:00
#SBATCH --mem=60gb
#SBATCH --output=count_kmers_ec.out
#SBATCH --error=count_kmers_ec.err
module load gcc/12.1.0 jellyfish/2.3.0

K=31
HASH=200M     # adjust upward if you see "hash table full"
THREADS=12

for i in $(cat caeruleum_males.list caeruleum_females.list); do
  zcat ${i}_R1.fastq.gz ${i}_R2.fastq.gz \
  | jellyfish count -m $K -C -s $HASH -t $THREADS -o ${i}.jf /dev/fd/0

  jellyfish dump -c ${i}.jf > ${i}.kmers.tsv
done
