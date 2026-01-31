#!/bin/bash
#SBATCH --job-name=Y_shared
#SBATCH --nodes=1
#SBATCH --partition=tier2q
#SBATCH --cpus-per-task=24
#SBATCH --time=168:00:00
#SBATCH --mem=124gb
#SBATCH --output=Y_shared.out
#SBATCH --error=Y_shared.err

module load gcc/12.1.0 jellyfish/2.3.0
set -euo pipefail

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

# Clean up temp dir if we created it (safe even if SLURM_TMPDIR)
trap '[[ "${SORTTMP}" == *sorttmp_* ]] && rm -rf "$SORTTMP" || true' EXIT
echo "Using SORTTMP=$SORTTMP"
echo "Using PAR=$PAR"

NMALES=$(wc -l < spectabile_males.list)
THRESH=$(( (8*NMALES + 9) / 10 ))
echo "NMALES=$NMALES  THRESH=$THRESH"

cat $(awk '{print $1".Ycand.kmers"}' spectabile_males.list) \
| LC_ALL=C sort --temporary-directory="$SORTTMP" -S 8G --parallel="$PAR" \
| uniq -c \
| awk -v t="$THRESH" '$1>=t {print $2}' \
> Y_SHARED_80pct.kmers
