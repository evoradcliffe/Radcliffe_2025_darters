#!/bin/bash
#SBATCH --job-name=map_female
#SBATCH --nodes=1
#SBATCH --partition=tier2q
#SBATCH --ntasks-per-node=24
#SBATCH --time=168:00:00
#SBATCH --mem=124gb
#SBATCH --output=mapfemale_es.out
#SBATCH --error=mapfemale_es.err

module load gcc/11.3.0 bedtools/2.30.0 samtools/1.18 bwa/0.7.17

#awk '{print ">"NR"\n"$1}' Y_SHARED_80pct.kmers > Y_SHARED_80pct.fa

awk '{print ">"NR"\n"$1}' FEMALE_UNIVERSE.kmers > FEMALE_UNIVERSE.fa

#bwa mem -t 12 /home/wrad07/rmoran-lab/projects/ps6/reference/GCF_008692095.1_UIUC_Espe_1.0_genomic.fna FEMALE_UNIVERSE.fa \
#| samtools view -bS - \
#| samtools sort -o female.bam
#samtools index female.bam

#bedtools genomecov -ibam female.bam -bg > femaele.bedgraph
# or directly count kmer hits per window:
#bedtools makewindows -g /home/wrad07/rmoran-lab/projects/ps6/reference/GCF_008692095.1_UIUC_Espe_1.0_genomic.fna.fai -w 50000 > windows_50kb.bed
#bedtools coverage -a windows_50kb.bed -b female.bam > female_density_50kb.tsv

##bins the genome into fixed-size windows (50 kb here) and counts how many Y-mers fall into each window, instead of base-by-base coverage, you get regional density.

REF=/home/wrad07/rmoran-lab/projects/ps6/reference/GCF_008692095.1_UIUC_Espe_1.0_genomic.fna
FAI=${REF}.fai

#bedtools coverage -sorted -g "$FAI" \
#  -a windows_50kb.sorted.bed \
#  -b female.q20.primary.bam \
#> female_density_50kb.tsv

#wc -l female_density_50kb.tsv
#head female_density_50kb.tsv

bedtools makewindows -g "$FAI" -w 50000 > windows_50kb.bed
bedtools sort -faidx "$FAI" -i windows_50kb.bed > windows_50kb.faidxsorted.bed

bedtools coverage -sorted -g "$FAI" \
  -a windows_50kb.faidxsorted.bed \
  -b female.q20.primary.bam \
> female_density_50kb.tsv
