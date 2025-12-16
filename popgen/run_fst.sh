module load GCC/13.2.0 VCFtools/0.1.16

vcftools --gzvcf ../merged_RB2RB_ChrOnly_filtered_RepsRem_Indels10bpRem_nohybrids.vcf.gz --weir-fst-pop male_sym_IL.txt --weir-fst-pop female_sym_IL.txt --out sym_IL_mvf
