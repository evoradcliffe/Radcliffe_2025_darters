module load GCC/13.2.0 VCFtools/0.1.16

vcftools --vcf your_sympatric.vcf \
  --weir-fst-pop male_samples.txt \
  --weir-fst-pop female_samples.txt \
  --out MvF_fst

### Clean up data to remove any NaN, NA, etc.:
awk '!/nan/ {print $1, $2, $3}' MvF_fst.weir.fst > MvF_fst_clean.txt

### If you need to change scaffolds to their corresponding chromosomes (i.e., from the reference genome), run:
awk 'NR==FNR { map[$1]=$2; next }
     FNR==1 { print; next }
     {
       if ($1 in map) {
         $1 = map[$1]
       } else {
         print "WARNING: No mapping for", $1 > "/dev/stderr"
       }
       print
     }' scaffold_to_chromosome_map.txt OFS="\t" mvf.weir_clean.fst > mvf.chrom_fst
