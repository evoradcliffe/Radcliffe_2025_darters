#!/bin/bash
#SBATCH --job-name=combine_gvcf
#SBATCH --nodes=1
#SBATCH --partition=tier2q
#SBATCH --ntasks-per-node=24
#SBATCH --time=168:00:00
#SBATCH --mem=124gb
#SBATCH --output=combine_frag.out
#SBATCH --error=combine_frag.err

echo -n "Ran: "
date

# Load required module
module load gcc/12.1.0 gatk4/4.4.0

REF="/home/wrad07/rmoran-lab/projects/ps6/reference/GCF_008692095.1_UIUC_Espe_1.0_genomic.fna"
GVCF_DIR="/home/wrad07/rmoran-lab/projects/darter_WGS/gvcfs_2025/fragi"
OUTPUT="/home/wrad07/rmoran-lab/projects/darter_WGS/gvcfs_2025/fragi/efragi_combined.g.vcf"

gatk --java-options "-Xmx60g" CombineGVCFs \
     -R "$REF" \
     $(find "$GVCF_DIR" -name '*.g.vcf.gz' | sort | sed 's/^/--variant /') \
     -O "$OUTPUT"

### CHECK MALFORMED ###
#for f in *.g.vcf.gz; do   line=$(bcftools view -H -r NC_045733.1:29187 "$f" 2>/dev/null | head -n 1);   if [ -n "$line" ]; then     alt=$(echo "$line" | cut -f5);     if [[ "$alt" != *"<NON_REF>"* ]]; then       echo "BAD (no <NON_REF>): $f   ALT=$alt";     fi;   fi; done

#gatk CombineGVCFs \
#  -R "$REF" \
#  $(cat gvcf_inputs.list | sed 's/^/--variant /') \
#  -O fragi.combined.g.vcf.gz \
#  -XL bad_site.bed

echo -n "Done: "
date
