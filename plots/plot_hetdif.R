library(ggplot2)
library(dplyr)

#HetDifs all Chr (no unplaced scaffs)
#EspecPopGenStats<- read.csv("Espec_Chr_Het_Difs.csv")
EcaerPopGenStats<- read.table("FIXED_sdr_Heterozygosity_bysex_nohybrids_RB2RB_difs.txt", header=TRUE)

# Bin positions into windows of 5000 base pairs
window_size <- 1


#E caeruleum plots

# Loop through chromosomes 1 to 24
for (chr_num in 1:24) {
  
  # Generate chromosome name dynamically (e.g., Chr1 is "NC_045733.1", Chr24 is "NC_045756.1")
  chr_name <- paste0("NC_0457", 32 + chr_num, ".1")  # Adjust this to match chromosome names

  # Filter and summarize for Sym_Het_Sex_Dif
  EcaerSym_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(diff_sym_IL_caeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Sym_Het_Sex_Dif
  pdf(paste0("EcaerSym_HetDif_bySite_Chr", chr_num, ".pdf"), width = 11, height = 5)
  p1=ggplot(EcaerSym_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_point() +  scale_y_continuous(limits=c(-0.505,0.505)) + 
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer IL Value across Chromosome", chr_num, "in 1000 bp Windows")) +
    theme_minimal()
  print(p1)
  dev.off()

  # Filter and summarize for Allo_Het_Sex_Dif
  EcaerAllo_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(diff_allo_MI_Ecaeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Allo_Het_Sex_Dif
  pdf(paste0("EcaerAlloMI_HetDif_bySite_Chr", chr_num, ".pdf"), width = 11, height = 5)
  p2=ggplot(EcaerAllo_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_point() +  scale_y_continuous(limits=c(-0.505,0.505)) + 
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer MI Value across Chromosome", chr_num, "in 1000 bp Windows")) +
    theme_minimal()
  print(p2)
  dev.off()

 # Filter and summarize for Allo_Het_Sex_Dif
  EcaerAlloLA_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(diff_allo_LA_Ecaeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Allo_Het_Sex_Dif
  pdf(paste0("EcaerAlloLA_HetDif_bySite_Chr", chr_num, ".pdf"), width = 11, height = 5)
  p3=ggplot(EcaerAlloLA_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_point() +  scale_y_continuous(limits=c(-0.505,0.505)) + 
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer LA Value across Chromosome", chr_num, "in 1000 bp Windows")) +
    theme_minimal()
  print(p3)
  dev.off()
}
