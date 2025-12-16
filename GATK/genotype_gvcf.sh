#!/bin/bash
#SBATCH --job-name=genotype_gvcf_colletteii
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=72:00:00
#SBATCH --mem=60gb
#SBATCH --output=genotype_coll.out
#SBATCH --error=genotype_coll.err
#SBATCH --array=0-29 ### CHANGE DEPENDING ON NUMBER OF CHROMOSOMES / REGIONS IN REGFOF

#   Prior to running this step we generated a file-of-files that lists the regions to operate over
#   We wanted variant calling by chromosome.
#   Each of chromosomes 1-(24) have their own file and unmapped contigs have their own file

#   Print when the script began running 
echo -n "Ran: "
date

#   Load modules and dependencies
module load gcc/12.1.0 gatk4/4.4.0

#   Set paths for GATK program file, reference genome, GATK regions file of files, combined GCVF files, and output directory
REF="/home/wrad07/rmoran-lab/projects/ps6/reference/GCF_008692095.1_UIUC_Espe_1.0_genomic.fna"
REGFOF="/home/wrad07/rmoran-lab/projects/darter_WGS/refs/rainbow_chrom_regions"
GVCF_DIR="/home/wrad07/rmoran-lab/projects/darter_WGS/gvcfs_2025/collettei"
OUTPUT_DIR="/home/wrad07/rmoran-lab/projects/darter_WGS/gvcfs_2025/genotyped"

export TMPDIR=/gpfs/data/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125/tmp
mkdir -p ${TMPDIR}

#   Build the sample list
SAMPLE_LIST=($(find ${GVCF_DIR} -maxdepth 1 -name '*.g.vcf.gz' | sort -V))
#   Put them into a format that will be accepted by the GATK command line
GATK_IN=()
for s in "${SAMPLE_LIST[@]}"
do
    GATK_IN+=("-V $s")
done

#   We define the output region based on the basename of the intervals file that specifies the chromosome or group of unscaffolded reads
#   This was read from the file-of-files
REGION=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${REGFOF})
REG_BNAME=$(basename ${REGION})
REG_OUT=${REG_BNAME/.intervals/}

#   Print these to check to make sure it is working
echo ${REGION}
echo ${REG_BNAME}
echo ${REG_OUT}
echo ${GATK_IN}

#   R = reference genome
#   L = intervals
#   variant = input the sample list
#   includeNonVariantSite = flag to include non-variant sites
#   heterozygosity = set the heterozygosity
#   O = output file

gatk --java-options "-Xmx60g -Djava.io.tmpdir=${TMPDIR}"\
    GenotypeGVCFs \
    -R ${REF} \
    -L ${REGION} \
    ${GATK_IN[@]} \
    --include-non-variant-sites \
    --heterozygosity 0.005 \
    --allow-old-rms-mapping-quality-annotation-data \
    -O ${OUTPUT_DIR}/${REG_OUT}.vcf


echo -n "Done: "
date
