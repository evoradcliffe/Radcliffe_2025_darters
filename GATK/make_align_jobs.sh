#!/bin/bash

#Originally copied from https://github.umn.edu/konox006/SEM_CaveFish/blob/master/Remapping/Scripts/Analysis/Make_Align_Jobs.sh by Tom Kono on 17 July 2018. Lightly modified by Adam Herman.

# Define software paths
module load GCC/12.2.0 seqtk/1.4

##SYMPATRIC E. caeruleum##
# Define output directory
#OUTDIR="/home/wrad07/moranlab/shared/darters/Ready_To_Align/symIL_Ecaer_reads/sym_Ecaer_10Xhap1"

# Define the reference index base      
#REFBASE="/home/wrad07/moranlab/shared/darters/10x_Data/scaffold_assemblies/ragtag10x/symEcaer_pseudohap1_output/symEcaer_pseudohap1.scaffold.fasta"

# Define the root of all reads                                                                                                                
#READS_ROOT="/home/wrad07/moranlab/shared/darters/Ready_To_Align/symIL_Ecaer_reads"

# How many threads do we want?                                                                                                                
#THREADS="16"

#### REPLACE ALL SAMPLE NAMES AND THEIR CORRESPONDING ENTRIES IN THE "declare -A READS"  SECTION WITH ALL E. SPECTABILE SAMPLES TO REPEAT FOR THAT SPECIES ####

#SAMPLES=(
#    "symEcaer_f1" "symEcaer_f2" "symEcaer_f3" "symEcaer_f4" "symEcaer_f5" "symEcaer_f6" "symEcaer_f7" "symEcaer_f8" "symEcaer_f9" "symEcaer_f10" "symEcaer_f11" "symEcaer_f12"
    
#    "symEcaer_m1" "symEcaer_m2" "symEcaer_m3" "symEcaer_m4" "symEcaer_m5" "symEcaer_m6" "symEcaer_m7" "symEcaer_m8" "symEcaer_m9" "symEcaer_m10" "symEcaer_m11" "symEcaer_m12"
#)

#declare -A READS
#READS["symEcaer_f1_R1"]=${READS_ROOT}/symIL-caeruleum-f1_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f1_R2"]=${READS_ROOT}/symIL-caeruleum-f1_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f1_U"]=${READS_ROOT}/symIL-caeruleum-f1_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f2_R1"]=${READS_ROOT}/symIL-caeruleum-f2_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f2_R2"]=${READS_ROOT}/symIL-caeruleum-f2_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f2_U"]=${READS_ROOT}/symIL-caeruleum-f2_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f3_R1"]=${READS_ROOT}/symIL-caeruleum-f3_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f3_R2"]=${READS_ROOT}/symIL-caeruleum-f3_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f3_U"]=${READS_ROOT}/symIL-caeruleum-f3_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f4_R1"]=${READS_ROOT}/symIL-caeruleum-f4_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f4_R2"]=${READS_ROOT}/symIL-caeruleum-f4_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f4_U"]=${READS_ROOT}/symIL-caeruleum-f4_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f5_R1"]=${READS_ROOT}/symIL-caeruleum-f5_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f5_R2"]=${READS_ROOT}/symIL-caeruleum-f5_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f5_U"]=${READS_ROOT}/symIL-caeruleum-f5_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f6_R1"]=${READS_ROOT}/symIL-caeruleum-f6_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f6_R2"]=${READS_ROOT}/symIL-caeruleum-f6_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f6_U"]=${READS_ROOT}/symIL-caeruleum-f6_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f7_R1"]=${READS_ROOT}/symIL-caeruleum-f7_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f7_R2"]=${READS_ROOT}/symIL-caeruleum-f7_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f7_U"]=${READS_ROOT}/symIL-caeruleum-f7_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f8_R1"]=${READS_ROOT}/symIL-caeruleum-f8_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f8_R2"]=${READS_ROOT}/symIL-caeruleum-f8_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f8_U"]=${READS_ROOT}/symIL-caeruleum-f8_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f9_R1"]=${READS_ROOT}/symIL-caeruleum-f9_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f9_R2"]=${READS_ROOT}/symIL-caeruleum-f9_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f9_U"]=${READS_ROOT}/symIL-caeruleum-f9_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f10_R1"]=${READS_ROOT}/symIL-caeruleum-f10_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f10_R2"]=${READS_ROOT}/symIL-caeruleum-f10_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f10_U"]=${READS_ROOT}/symIL-caeruleum-f10_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f11_R1"]=${READS_ROOT}/symIL-caeruleum-f11_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f11_R2"]=${READS_ROOT}/symIL-caeruleum-f11_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f11_U"]=${READS_ROOT}/symIL-caeruleum-f11_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_f12_R1"]=${READS_ROOT}/symIL-caeruleum-f12_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_f12_R2"]=${READS_ROOT}/symIL-caeruleum-f12_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_f12_U"]=${READS_ROOT}/symIL-caeruleum-f12_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m1_R1"]=${READS_ROOT}/symIL-caeruleum-m1_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m1_R2"]=${READS_ROOT}/symIL-caeruleum-m1_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m1_U"]=${READS_ROOT}/symIL-caeruleum-m1_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m2_R1"]=${READS_ROOT}/symIL-caeruleum-m2_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m2_R2"]=${READS_ROOT}/symIL-caeruleum-m2_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m2_U"]=${READS_ROOT}/symIL-caeruleum-m2_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m3_R1"]=${READS_ROOT}/symIL-caeruleum-m3_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m3_R2"]=${READS_ROOT}/symIL-caeruleum-m3_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m3_U"]=${READS_ROOT}/symIL-caeruleum-m3_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m4_R1"]=${READS_ROOT}/symIL-caeruleum-m4_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m4_R2"]=${READS_ROOT}/symIL-caeruleum-m4_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m4_U"]=${READS_ROOT}/symIL-caeruleum-m4_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m5_R1"]=${READS_ROOT}/symIL-caeruleum-m5_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m5_R2"]=${READS_ROOT}/symIL-caeruleum-m5_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m5_U"]=${READS_ROOT}/symIL-caeruleum-m5_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m6_R1"]=${READS_ROOT}/symIL-caeruleum-m6_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m6_R2"]=${READS_ROOT}/symIL-caeruleum-m6_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m6_U"]=${READS_ROOT}/symIL-caeruleum-m6_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m7_R1"]=${READS_ROOT}/symIL-caeruleum-m7_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m7_R2"]=${READS_ROOT}/symIL-caeruleum-m7_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m7_U"]=${READS_ROOT}/symIL-caeruleum-m7_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m8_R1"]=${READS_ROOT}/symIL-caeruleum-m8_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m8_R2"]=${READS_ROOT}/symIL-caeruleum-m8_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m8_U"]=${READS_ROOT}/symIL-caeruleum-m8_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m9_R1"]=${READS_ROOT}/symIL-caeruleum-m9_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m9_R2"]=${READS_ROOT}/symIL-caeruleum-m9_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m9_U"]=${READS_ROOT}/symIL-caeruleum-m9_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m10_R1"]=${READS_ROOT}/symIL-caeruleum-m10_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m10_R2"]=${READS_ROOT}/symIL-caeruleum-m10_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m10_U"]=${READS_ROOT}/symIL-caeruleum-m10_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m11_R1"]=${READS_ROOT}/symIL-caeruleum-m11_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m11_R2"]=${READS_ROOT}/symIL-caeruleum-m11_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m11_U"]=${READS_ROOT}/symIL-caeruleum-m11_adtrim_trim_unpair_all.fastq.gz
#READS["symEcaer_m12_R1"]=${READS_ROOT}/symIL-caeruleum-m12_adtrim_trim_pair_R1.fastq.gz
#READS["symEcaer_m12_R2"]=${READS_ROOT}/symIL-caeruleum-m12_adtrim_trim_pair_R2.fastq.gz
#READS["symEcaer_m12_U"]=${READS_ROOT}/symIL-caeruleum-m12_adtrim_trim_unpair_all.fastq.gz

#for s in ${SAMPLES[@]}
#do
    # Define output file                                                                                                                        
#    outfile=${OUTDIR}/${s}_raw.bam
    # Get the R1, R2, Un read paths                                                                                                             
#    r1=${READS["${s}_R1"]}
#    r2=${READS["${s}_R2"]}
#    un=${READS["${s}_U"]}
    # echo the commands. Note that we merge the R1 and R2 files, then append the                                                                
    # unpaired reads                                                                                                                            
#    echo "(seqtk mergepe ${r1} ${r2}; gzip -cd ${un}) | bwa mem -t ${THREADS} -k 12 -M -p ${REFBASE} - | samtools view -hbu - | samtools sort -o ${outfile} -@ ${THREADS} -"
#done

#ALLOPATRIC E. caeruleum##
OUTDIR="/home/wrad07/moranlab/shared/darters/Ready_To_Align/allo_Ecaer_reads/allo_Ecaer_10Xhap1/raw_bams"

# Define the reference index base                  
REFBASE="/home/wrad07/moranlab/shared/darters/10x_Data/Assemblies_UMGC/Allo_Ecaer_10X_pseudohap_assemblies/Allo_Ecaer_pseudohap2.1.fasta"

# Define the root of all reads
READS_ROOT="/home/wrad07/moranlab/shared/darters/Ready_To_Align/allo_Ecaer_reads/michigan"

# How many threads do we want?                                                            
THREADS="16"

SAMPLES=(
    "alloEcaer_MI_f1" "alloEcaer_MI_f2" "alloEcaer_MI_f3" "alloEcaer_MI_f4" "alloEcaer_MI_f5" "alloEcaer_MI_f6" "alloEcaer_MI_f7"                                                                                                                 
    "alloEcaer_MI_m1" "alloEcaer_MI_m2" "alloEcaer_MI_m3" "alloEcaer_MI_m4" "alloEcaer_MI_m5" "alloEcaer_MI_m6" "alloEcaer_MI_m7"
)

declare -A READS 
READS["alloEcaer_MI_f1_R1"]=${READS_ROOT}/alloMI-caeruleum-f1_adtrim_trim_pair_R1.fastq.gz
READS["alloEcaer_MI_f1_R2"]=${READS_ROOT}/alloMI-caeruleum-f1_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_f1_U"]=${READS_ROOT}/alloMI-caeruleum-f1_adtrim_trim_unpair_all.fastq.gz
READS["alloEcaer_MI_f2_R1"]=${READS_ROOT}/alloMI-caeruleum-f2_adtrim_trim_pair_R1.fastq.gz 
READS["alloEcaer_MI_f2_R2"]=${READS_ROOT}/alloMI-caeruleum-f2_adtrim_trim_pair_R2.fastq.gz
READS["alloEcaer_MI_f2_U"]=${READS_ROOT}/alloMI-caeruleum-f2_adtrim_trim_unpair_all.fastq.gz
READS["alloEcaer_MI_f3_R1"]=${READS_ROOT}/alloMI-caeruleum-f3_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_f3_R2"]=${READS_ROOT}/alloMI-caeruleum-f3_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_f3_U"]=${READS_ROOT}/alloMI-caeruleum-f3_adtrim_trim_unpair_all.fastq.gz
READS["alloEcaer_MI_f4_R1"]=${READS_ROOT}/alloMI-caeruleum-f4_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_f4_R2"]=${READS_ROOT}/alloMI-caeruleum-f4_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_f4_U"]=${READS_ROOT}/alloMI-caeruleum-f4_adtrim_trim_unpair_all.fastq.gz
READS["alloEcaer_MI_f5_R1"]=${READS_ROOT}/alloMI-caeruleum-f5_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_f5_R2"]=${READS_ROOT}/alloMI-caeruleum-f5_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_f5_U"]=${READS_ROOT}/alloMI-caeruleum-f5_adtrim_trim_unpair_all.fastq.gz
READS["alloEcaer_MI_f6_R1"]=${READS_ROOT}/alloMI-caeruleum-f6_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_f6_R2"]=${READS_ROOT}/alloMI-caeruleum-f6_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_f6_U"]=${READS_ROOT}/alloMI-caeruleum-f6_adtrim_trim_unpair_all.fastq.gz
READS["alloEcaer_MI_f7_R1"]=${READS_ROOT}/alloMI-caeruleum-f7_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_f7_R2"]=${READS_ROOT}/alloMI-caeruleum-f7_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_f7_U"]=${READS_ROOT}/alloMI-caeruleum-f7_adtrim_trim_unpair_all.fastq.gz
READS["alloEcaer_MI_m1_R1"]=${READS_ROOT}/alloMI-caeruleum-m1_adtrim_trim_pair_R1.fastq.gz
READS["alloEcaer_MI_m1_R2"]=${READS_ROOT}/alloMI-caeruleum-m1_adtrim_trim_pair_R2.fastq.gz 
READS["alloEcaer_MI_m1_U"]=${READS_ROOT}/alloMI-caeruleum-m1_adtrim_trim_unpair_all.fastq.gz                                                     
READS["alloEcaer_MI_m2_R1"]=${READS_ROOT}/alloMI-caeruleum-m2_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_m2_R2"]=${READS_ROOT}/alloMI-caeruleum-m2_adtrim_trim_pair_R2.fastq.gz 
READS["alloEcaer_MI_m2_U"]=${READS_ROOT}/alloMI-caeruleum-m2_adtrim_trim_unpair_all.fastq.gz
READS["alloEcaer_MI_m3_R1"]=${READS_ROOT}/alloMI-caeruleum-m3_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_m3_R2"]=${READS_ROOT}/alloMI-caeruleum-m3_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_m3_U"]=${READS_ROOT}/alloMI-caeruleum-m3_adtrim_trim_unpair_all.fastq.gz 
READS["alloEcaer_MI_m4_R1"]=${READS_ROOT}/alloMI-caeruleum-m4_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_m4_R2"]=${READS_ROOT}/alloMI-caeruleum-m4_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_m4_U"]=${READS_ROOT}/alloMI-caeruleum-m4_adtrim_trim_unpair_all.fastq.gz  
READS["alloEcaer_MI_m5_R1"]=${READS_ROOT}/alloMI-caeruleum-m5_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_m5_R2"]=${READS_ROOT}/alloMI-caeruleum-m5_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_m5_U"]=${READS_ROOT}/alloMI-caeruleum-m5_adtrim_trim_unpair_all.fastq.gz                   
READS["alloEcaer_MI_m6_R1"]=${READS_ROOT}/alloMI-caeruleum-m6_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_m6_R2"]=${READS_ROOT}/alloMI-caeruleum-m6_adtrim_trim_pair_R2.fastq.gz                                                       
READS["alloEcaer_MI_m6_U"]=${READS_ROOT}/alloMI-caeruleum-m6_adtrim_trim_unpair_all.fastq.gz                                                     
READS["alloEcaer_MI_m7_R1"]=${READS_ROOT}/alloMI-caeruleum-m7_adtrim_trim_pair_R1.fastq.gz                                                       
READS["alloEcaer_MI_m7_R2"]=${READS_ROOT}/alloMI-caeruleum-m7_adtrim_trim_pair_R2.fastq.g#z                                                       
READS["alloEcaer_MI_m7_U"]=${READS_ROOT}/alloMI-caeruleum-m7_adtrim_trim_unpair_all.fastq.gz   

for s in ${SAMPLES[@]}
do
    # Define output file                                                                                            
    outfile=${OUTDIR}/${s}_raw.bam 
    r1=${READS["${s}_R1"]}
    r2=${READS["${s}_R2"]}
    un=${READS["${s}_U"]}
 
    echo "(seqtk mergepe ${r1} ${r2}; gzip -cd ${un}) | bwa mem -t ${THREADS} -k 12 -M -p ${REFBASE} - | samtools view -hbu - | samtools sort -o ${outfile} -@ ${THREADS} -"
done
