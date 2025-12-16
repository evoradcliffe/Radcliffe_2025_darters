#!/bin/bash
#SBATCH --job-name=use_cutadapt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH --mem=62gb
#SBATCH --output=cutadapt.out
#SBATCH --error=cutadapt.err

#   Set the folder where the analysis is being done
WORKDIR=/gpfs/data/rmoran-lab/projects/darter_WGS/Raw_Lanes_Concatenated_2024/
INDIR=${WORKDIR}
OUTDIR=${WORKDIR}/CutAdapt_Out

export WORKDIR INDIR OUTDIR

# Make sure output directory exists
mkdir -p "${OUTDIR}"

#   Specify the adapters used
#   Library prep completed using PerkinElmer NEXTFLEX Rapid XP DNA-Seq Kit HT and Unique Dual Index Barcodes
#   Sequenced using: Illumina NovaSeq 6000 S4 - 2x150 v1.5
#   Barcoded Adapters manual: https://perkinelmer-appliedgenomics.com/wp-content/uploads/2021/11/514150-NEXTflex-Unique-Dual-Index-Barcodes-21-06-SetA_1121.pdf
#   Library Prep Kit manual: https://perkinelmer-appliedgenomics.com/wp-content/uploads/marketing/NEXTFLEX/rapid_xp/NEXTFLEX-Rapid-XP-v2-DNA-Seq-Kit-AG072205_24_MA-V23.08_PE.pdf
#   These are sequences that match the Truseq Dual Index Library (https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html)
export R1_adapter="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
export R2_adapter="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"

#Load the software needed (parallel, cutadapt, and dependencies)
module load gcc/12.1.0 parallel/20240622 cutadapt/4.5

#   Define a bash function (a small script within a script) for doing the processing
#   We want to do this because we need to do the same thing for each file
#   Later, we basically tell the computer to do this script for every file
#   Our function is called parallelcut
parallelcut() {
    #Get sample ID from SampleID (first argument supplied to parallel with -a)
    SAMP=${1}

    #Use cutadapt to cut the adapters off
    #-e tells it the maximum error rate
    #-a, -g tells it the adapter to cut (we set the adapter above) (-a is 3' adapter -g is 5' adapter)
    #The input file goes at the end of the command
    #Specify the output file after the > for standard output stream of running a script
    cutadapt -e 0.01 --match-read-wildcards -a ${R1_adapter} -o ${HOME}CutAdapt_Out/${SAMP}_adtrim_R1.fastq.gz ${HOME}Raw_Lanes_Concatenated/${SAMP}_R1.fastq.gz
    cutadapt -e 0.01 --match-read-wildcards -a ${R2_adapter} -o ${HOME}CutAdapt_Out/${SAMP}_adtrim_R2.fastq.gz ${HOME}Raw_Lanes_Concatenated/${SAMP}_R2.fastq.gz
} 


# === Function for Cutadapt ===
parallelcut() {
    SAMP=${1}
    echo "=== Running Cutadapt on sample: ${SAMP} ==="

    cutadapt \
      -e 0.01 --match-read-wildcards \
      -a "$R1_adapter" \
      -A "$R2_adapter" \
      -q 20,20 \
      -m 40 \
      -o "${OUTDIR}/${SAMP}_adtrim_R1.fastq.gz" \
      -p "${OUTDIR}/${SAMP}_adtrim_R2.fastq.gz" \
      "${INDIR}/${SAMP}_R1.fastq.gz" \
      "${INDIR}/${SAMP}_R2.fastq.gz"

    echo "=== Finished ${SAMP} ==="
}

export -f parallelcut

cd "${WORKDIR}"

parallel --joblog cutadapt_parallel_logfile.txt -a SampleID.txt parallelcut
