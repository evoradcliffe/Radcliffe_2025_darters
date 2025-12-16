#!/bin/bash

#### AGAIN, REPLACE ALL SAMPLES WITH E. CAERULEUM AND E. SPECTABILE SAMPLES BEING ANALYZED / THEIR GIVEN DIRECTORIES ####  
#   These are the samples that we are analyzing. These will eventually be part of filenames, so we define them here, as a bash array.
SAMPLES=("Ecyanorum-BROK-F1" "Ecyanorum-BROK-F2" "Ecyanorum-BROK-F3" "Ecyanorum-BROK-F4" "Ecyanorum-BROK-F5" "Ecyanorum-BROK-M1" "Ecyanorum-BROK-M2" "Ecyanorum-BROK-M3" "Ecyanorum-BROK-M4" "Ecyanorum-BROK-M5"

         "Epulch-BROK-F1" "Epulch-BROK-F2" "Epulch-BROK-F3" "Epulch-BROK-F4" "Epulch-BROK-F5" "Epulch-BROK-M1" "Epulch-BROK-M2" "Epulch-BROK-M3" "Epulch-BROK-M4" "Epulch-BROK-M5"

         "Epulch-Brazos-F10" "Epulch-Brazos-F9" "Epulch-Brazos-M8" "Epulch-McCoy-F1" "Epulch-McCoy-F2" "Epulch-McCoy-F3" "Epulch-McCoy-F4" "Epulch-McCoy-M1" "Epulch-McCoy-M2" "Epulch-McCoy-M3" "Epulch-McCoy-M4"

         "Espectabile-Kankakee-F2" "Espectabile-Kankakee-F3" "Espectabile-Kankakee-M1" "Espectabile-Kankakee-M2"

        "Esquamosum-Neosho-F1" "Esquamosum-Neosho-F10" "Esquamosum-Neosho-F2" "Esquamosum-Neosho-F3" "Esquamosum-Neosho-F4" "Esquamosum-Neosho-F5" "Esquamosum-Neosho-F6" "Esquamosum-Neosho-F7" "Esquamosum-Neosho-F8" "Esquamosum-Neosho-F9" "Esquamosum-Neosho-M1" "Esquamosum-Neosho-M10" "Esquamosum-Neosho-M2" "Esquamosum-Neosho-M3" "Esquamosum-Neosho-M4" "Esquamosum-Neosho-M5" "Esquamosum-Neosho-M6" "Esquamosum-Neosho-M7" "Esquamosum-Neosho-M8" "Esquamosum-Neosho-M9"
)
#   Generate two commands for each sample to (1) mark duplicate and (2) add read groups (i.e. sample name)

module load openjdk/17.0.2 picard/3.4.0

for s in ${SAMPLES[@]}
do
    echo "java -Xmx4G -jar ${PICARD} AddOrReplaceReadGroups \
        TMP_DIR=/tmp \
        INPUT=/home/wrad07/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125/Aligned_Reads/${s}_raw.bam \
        OUTPUT=/home/wrad07/rmoran-lab/projects/darter_WGS/Raw_Reads_2025/25205Mra_N25125/Aligned_Reads/rgadd_out/${s}_rgadd.bam \
        RGID=${s} RGLB=${s} RGPL=Illumina RGPU=${s} RGSM=${s} \
        CREATE_INDEX=true"
done
