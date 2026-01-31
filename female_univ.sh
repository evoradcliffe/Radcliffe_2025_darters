#!/bin/bash
#SBATCH --job-name=female_univ_es
#SBATCH --nodes=1
#SBATCH --partition=tier1q
#SBATCH --ntasks-per-node=12
#SBATCH --time=168:00:00
#SBATCH --mem=60gb
#SBATCH --output=female_univ_es.out
#SBATCH --error=female_univ_es.err

module load gcc/12.1.0 jellyfish/2.3.0

### make female "k-verse", which acts as subtraction set, everything seen in females above certain threshold ###

set -euo pipefail

MINC=3
PAR="${SLURM_CPUS_PER_TASK:-12}"

# Choose a writable temp directory for sort
if [[ -n "${SLURM_TMPDIR:-}" && -d "$SLURM_TMPDIR" ]]; then
    SORTTMP="$SLURM_TMPDIR"
elif [[ -n "${SCRATCH:-}" && -d "$SCRATCH" ]]; then
    SORTTMP="$SCRATCH/sorttmp_${SLURM_JOB_ID}"
    mkdir -p "$SORTTMP"
else
    SORTTMP="$PWD/sorttmp_${SLURM_JOB_ID}"
    mkdir -p "$SORTTMP"
fi

echo "Using SORTTMP=$SORTTMP"
echo "Using PAR=$PAR"
echo "SLURM_CPUS_PER_TASK='${SLURM_CPUS_PER_TASK:-}'  SLURM_NTASKS='${SLURM_NTASKS:-}'"

cat *.F.kmers \
  | sort --temporary-directory="$SORTTMP" -S 8G --parallel="$PAR" \
  | uniq \
  > FEMALE_UNIVERSE.kmers

# Extract female kmers passing MINC
#for f in $(cat spectabile_females.list); do
#  awk -v m=$MINC '$2>=m {print $1}' ${f}.kmers.tsv | sort -S 4G > ${f}.F.kmers
#done
