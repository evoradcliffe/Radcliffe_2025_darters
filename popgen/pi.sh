module load GCC/13.2.0 VCFtools/0.1.16

vcftools --site-pi --gzvcf ../merged_RB2RB_ChrOnly_filtered_RepsRem_Indels10bpRem_nohybrids.vcf.gz --keep pooled_female.txt --out rbxrb_pooled_f_div
