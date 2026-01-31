#!/bin/bash
#SBATCH --job-name=candidateY_es
#SBATCH --nodes=1
#SBATCH --partition=tier1q
#SBATCH --ntasks-per-node=12
#SBATCH --time=168:00:00
#SBATCH --mem=60gb
#SBATCH --output=canY_es.out
#SBATCH --error=canY_es.err

### subtract universe from each male -> "candidate Y-mers per male" ###

module load gcc/12.1.0 jellyfish/2.3.0

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

for m in $(cat spectabile_males.list); do
  awk -v m=3 '$2>=m {print $1}' ${m}.kmers.tsv | sort --temporary-directory="$SORTTMP" -S 4G --parallel="$PAR" > ${m}.M.kmers

  comm -23 ${m}.M.kmers FEMALE_UNIVERSE.kmers > ${m}.Ycand.kmers
done
