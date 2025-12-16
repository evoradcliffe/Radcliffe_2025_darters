library("vcfR")

#Specify a population for each individual in your vcf, in the order they appear in the vcf
pop <- as.factor(c("F_allo_LA_Ecaeruleum",   "F_allo_LA_Ecaeruleum",     "F_allo_LA_Ecaeruleum",        "F_allo_LA_Ecaeruleum", "M_allo_LA_Ecaeruleum", "M_allo_LA_Ecaeruleum", "M_allo_LA_Ecaeruleum", "M_allo_LA_Ecaeruleum", "M_allo_LA_Ecaeruleum", "F_allo_MI_Ecaeruleum", "F_allo_MI_Ecaeruleum", "F_allo_MI_Ecaeruleum", "F_allo_MI_Ecaeruleum", "F_allo_MI_Ecaeruleum", "F_allo_MI_Ecaeruleum", "F_allo_MI_Ecaeruleum", "M_allo_MI_Ecaeruleum", "M_allo_MI_Ecaeruleum", "M_allo_MI_Ecaeruleum", "M_allo_MI_Ecaeruleum", "M_allo_MI_Ecaeruleum", "M_allo_MI_Ecaeruleum", "M_allo_MI_Ecaeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "F_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum", "M_sym_IL_caeruleum"))


#read in your vcf and write the output for each chromosome
vcf1 <- read.vcfR("merged_RB2RB_ChrOnly_filtered_RepsRem_Indels10bpRem_nohybrids.vcf.gz")
myDiff1 <- genetic_diff(vcf1, pops=pop, method='nei')
write.csv(myDiff1, 'sdr_Heterozygosity_bysex_nohybrids_RB2RB.csv')
