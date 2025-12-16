#!/bin/bash
#SBATCH --job-name=filter
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=12:00:00
#SBATCH --mem=102gb
#SBATCH --output=filter_snps_ind_gatk4.5.out
#SBATCH --error=filter_snps_ind_gatk4.5.err

module load GCCcore/12.3.0
module load GATK/4.5.0.0-Java-17

export _JAVA_OPTIONS="-Xmx100g" 
gatk VariantFiltration --tmp-dir /scratch/user/dkimfish --filter-name "QD2" --filter-expression "QD < 2.0" --filter-name "FS60" --filter-expression "FS > 60.0" --filter-name "MQ40" --filter-expression "MQ < 40.0" --filter-name "MQRankSum-12.5" --filter-expression "MQRankSum < -12.5" --filter-name "ReadPosRankSum-8" --filter-expression "ReadPosRankSum < -8.0" -R /scratch/group/moranlab/shared/darters/refs/Espec_sym_ref_NCBI/GCF_008692095.1_UIUC_Espe_1.0_genomic.fa -V /scratch/group/moranlab/shared/darters/RAW_WGS/Psychromaster_gvcf/genotyped/filtered/subset_snps/Psychromaster_cohort_subsetSNPs.vcf -O /scratch/group/moranlab/shared/darters/RAW_WGS/Psychromaster_gvcf/genotyped/filtered/filtered_snps/Psychromaster_cohort_filteredSNPs.vcf --verbosity INFO
gatk VariantFiltration --tmp-dir /scratch/user/dkimfish --filter-name "QD2" --filter-expression "QD < 2.0" --filter-name "FS60" --filter-expression "FS > 60.0" --filter-name "MQ40" --filter-expression "MQ < 40.0" -R /scratch/group/moranlab/shared/darters/RAW_WGS/GCA_026937815.1_UBuf_Eper_1.0_genomic.fna -V /scratch/group/moranlab/shared/darters/RAW_WGS/Psychromaster_gvcf/genotyped/filtered/subset_mono/Psychromaster_cohort_subsetMONO.vcf -O /scratch/group/moranlab/shared/darters/RAW_WGS/Psychromaster_gvcf/genotyped/filtered/filtered_mono/Psychromaster_cohort_filteredMONO.vcf --verbosity INFO
gatk VariantFiltration --tmp-dir /scratch/user/dkimfish --filter-name "QD2" --filter-expression "QD < 2.0" --filter-name "FS200" --filter-expression "FS > 200.0" --filter-name "ReadPosRankSum-20" --filter-expression "ReadPosRankSum < -20.0" -R /scratch/group/moranlab/shared/darters/RAW_WGS/GCA_026937815.1_UBuf_Eper_1.0_genomic.fna -V /scratch/group/moranlab/shared/darters/RAW_WGS/Psychromaster_gvcf/genotyped/filtered/subset_mixed/Psychromaster_cohort_subsetMixedIndels.vcf -O /scratch/group/moranlab/shared/darters/RAW_WGS/Psychromaster_gvcf/genotyped/filtered/filtered_mixed/Psychromaster_cohort_filteredMixedIndels.vcf --verbosity INFO
