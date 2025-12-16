#!/bin/bash
#SBATCH --job-name=haplotype_caller
#SBATCH --nodes=1
#SBATCH --partition=tier2q
#SBATCH --ntasks-per-node=24
#SBATCH --time=168:00:00
#SBATCH --mem=124gb
#SBATCH --output=haplotype_caller.out
#SBATCH --error=haplotype_caller.err
#SBATCH --array=0-15 #### change based on number of samples

echo -n "Ran: "
date

# Load required module
module load gcc/12.1.0 gatk4/4.4.0

# Define paths
REF="/home/wrad07/rmoran-lab/projects/ps6/reference/GCF_008692095.1_UIUC_Espe_1.0_genomic.fna"
ALN_DIR="/home/wrad07/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125/Aligned_Reads/dups_marked/mapped/fragi"
OUTDIR="/home/wrad07/rmoran-lab/projects/darter_WGS/gvcfs_2025"

export TMPDIR=/gpfs/data/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125/tmp
mkdir -p ${TMPDIR}

# Create output directory if needed
mkdir -p ${OUTDIR}

# Get sample list and current sample
SAMPLE_LIST=($(find ${ALN_DIR} -name '*_dupsmarked_mapped.bam' | sort -V))
CURRENT_SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

# Derive sample name
FNAME=$(basename ${CURRENT_SAMPLE})
SAMPLENAME=${FNAME/_dupsmarked_mapped.bam/}

# Run HaplotypeCaller (GATK 4 syntax)
gatk --java-options "-Xmx124g -Djava.io.tmpdir=${TMPDIR}" HaplotypeCaller \
  -R ${REF} \
  -I ${CURRENT_SAMPLE} \
  -O ${OUTDIR}/${SAMPLENAME}.g.vcf.gz \
  -ERC GVCF \
  --heterozygosity 0.005 \
  --native-pair-hmm-threads 24

echo -n "Done: "
date
