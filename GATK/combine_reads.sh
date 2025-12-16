GAT#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=72:00:00
#SBATCH --mem=62gb
#SBATCH --output=combine.out
#SBATCH --error=combine.err

# ==== Directory setup ====
WORKDIR=/home/wrad07/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125
INDIR=${WORKDIR}/CutAdapt_Out
OUTDIR=${WORKDIR}/Concatenated_Reads

export WORKDIR INDIR OUTDIR

module load gcc/12.1.0 parallel/20240622

# ==== Function to combine all lanes per sample ====
parcomb() {
    LANEID=${1} 
    # strip off "_L00X" to get the base sample name
    SAMP=$(echo "${LANEID}" | sed 's/_L00[0-9]//')

    # ensure output directory exists
    mkdir -p "${OUTDIR}"

    echo "Processing sample: ${SAMP}"

    # locate all R1 and R2 files across lanes for this sample
    IN_R1="${INDIR}/${SAMP}_L*_adtrim_R1.fastq.gz"
    IN_R2="${INDIR}/${SAMP}_L*_adtrim_R2.fastq.gz"

    OUT_R1="${OUTDIR}/${SAMP}_adtrim_R1.fastq.gz"
    OUT_R2="${OUTDIR}/${SAMP}_adtrim_R2.fastq.gz"

    # check existence before concatenating
    if ls ${IN_R1} 1> /dev/null 2>&1 && ls ${IN_R2} 1> /dev/null 2>&1; then
        cat ${IN_R1} > "${OUT_R1}"
        cat ${IN_R2} > "${OUT_R2}"
        echo " → ${SAMP} concatenation complete."
    else
        echo " Missing input files for ${SAMP}"
    fi
}

# export function for parallel
export -f parcomb

# run in parallel
cd "${WORKDIR}"
parallel --joblog combine_reads_log.txt -a Sample_IDs.txt parcomb
