library(ggplot2)
library(dplyr)

#setwd("~/Desktop/darters/SDR")

setwd("~/moranlab/shared/darters/fall2023darters_gvcf/clean_data/sdr_data")



#Het difs for only Chr 9
#EspecPopGenStats<- read.csv("Espec_Chr_Het_DifsONLY_Chr9.csv")
#EcaerPopGenStats<- read.csv("Espec_Chr_Het_DifsONLY_Chr9.csv")


#HetDifs all Chr (no unplaced scaffs)
EspecPopGenStats<- read.csv("Espec_Chr_Het_Difs.csv")
EcaerPopGenStats<- read.table("sdr_Heterozygosity_bysex_nohybrids_EcaerOnly_hetDifsOnly_Chr.txt", header=TRUE)

# Bin positions into windows of 5000 base pairs
window_size <- 5000


#E spectabile plots

# Loop through chromosomes 1 to 24
for (chr_num in 1:24) {
  
  # Generate chromosome name dynamically (e.g., Chr1 is "NC_045733.1", Chr24 is "NC_045756.1")
  chr_name <- paste0("NC_0457", 32 + chr_num, ".1")  # Adjust this to match chromosome names

  # Filter and summarize for Sym_Het_Sex_Dif
  EsSym_HetDif_binned <- EspecPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(Sym_Het_Sex_Dif, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Sym_Het_Sex_Dif
  pdf(paste0("EsSym_HetDif_binned_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EsSym_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_line() +
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Espec Sam Value across Chromosome", chr_num, "in 5000 bp Windows")) +
    theme_minimal()
  dev.off()

  # Filter and summarize for Allo_Het_Sex_Dif
  EsAllo_HetDif_binned <- EspecPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(Allo_Het_Sex_Dif, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Allo_Het_Sex_Dif
  pdf(paste0("EsAllo_HetDif_binned_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EsAllo_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_line() +
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Espec Allo Value across Chromosome", chr_num, "in 5000 bp Windows")) +
    theme_minimal()
  dev.off()
}


for (chr_num in 1:24) {
  
  # Generate chromosome name dynamically (e.g., Chr1 is "NC_045733.1", Chr24 is "NC_045756.1")
  chr_name <- paste0("NC_0457", 32 + chr_num, ".1")  # Adjust this to match chromosome names
  print(chr_name)
}

#E caeruleum plots

# Loop through chromosomes 1 to 24
for (chr_num in 1:24) {
  
  # Generate chromosome name (e.g., Chr1 is "NC_045733.1", Chr24 is "NC_045756.1")
  chr_name <- paste0("NC_0457", 32 + chr_num, ".1")  # Adjust this to match chromosome names

  # Filter and summarize for Sym_Het_Sex_Dif
  EcaerSym_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(M_F_HetDif_sym_IL_Ecaeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Sym_Het_Sex_Dif
  pdf(paste0("EcaerSym_HetDif_binned_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EcaerSym_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_line() +
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer IL Value across Chromosome", chr_num, "in 5000 bp Windows")) +
    theme_minimal()
  dev.off()

  # Filter and summarize for Allo_Het_Sex_Dif
  EcaerAllo_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(M_F_HetDif_allo_MI_Ecaeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Allo_Het_Sex_Dif
  pdf(paste0("EcaerAlloMI_HetDif_binned_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EcaerAllo_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_line() +
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer MI Value across Chromosome", chr_num, "in 5000 bp Windows")) +
    theme_minimal()
  dev.off()
 # Filter and summarize for Allo_Het_Sex_Dif
  EcaerAlloLA_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(M_F_HetDif_allo_LA_Ecaeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Allo_Het_Sex_Dif
  pdf(paste0("EcaerAlloLA_HetDif_binned_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EcaerAlloLA_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_line() +
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer LA Value across Chromosome", chr_num, "in 5000 bp Windows")) +
    theme_minimal()
  dev.off()
}



###################### chr 9 and 23 only
for (chr_num in 1:24) {
  
  # Generate chromosome name dynamically (e.g., Chr1 is "NC_045733.1", Chr24 is "NC_045756.1")
  chr_name <- paste0("NC_0457", 32 + chr_num, ".1")  # Adjust this to match chromosome names
  
  # Filter and summarize for Sym_Het_Sex_Dif
  EsSym_HetDif_binned <- EspecPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(Sym_Het_Sex_Dif, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Sym_Het_Sex_Dif
  pdf(paste0("EsSym_HetDif_bySite_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EsSym_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_point() + scale_y_continuous(limits=c(-0.501,0.501)) + 
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Espec Sym Value across Chromosome", chr_num, "by site")) +
    theme_minimal()
  dev.off()
  
  
  
  
  # Filter and summarize for Allo_Het_Sex_Dif
  EsAllo_HetDif_binned <- EspecPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(Allo_Het_Sex_Dif, na.rm = TRUE))  # Summarize with NA handling
  
  
  pdf(paste0("EsAllo_HetDif_bySite_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EsAllo_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_point() + scale_y_continuous(limits=c(-0.501,0.501)) + 
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Espec Allo Value across Chromosome", chr_num, "by site")) +
    theme_minimal()
  dev.off()
  
  
  
#E caeruleum #
  
  
  EcaerSym_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(M_F_HetDif_sym_IL_Ecaeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Sym_Het_Sex_Dif
  pdf(paste0("EcaerSym_HetDif_bySite_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EcaerSym_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_point() + scale_y_continuous(limits=c(-0.501,0.501)) + 
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer IL Value across Chromosome", chr_num, "by Site")) +
    theme_minimal()
  dev.off()
  
  # Filter and summarize for Allo_Het_Sex_Dif
  EcaerAllo_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(M_F_HetDif_allo_MI_Ecaeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Allo_Het_Sex_Dif
  pdf(paste0("EcaerAlloMI_HetDif_bySite_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EcaerAllo_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_point() + scale_y_continuous(limits=c(-0.501,0.501)) + 
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer MI Value across Chromosome", chr_num, "by site")) +
    theme_minimal()
  dev.off()
  
  # Filter and summarize for Allo_Het_Sex_Dif
  EcaerAlloLA_HetDif_binned <- EcaerPopGenStats %>%
    filter(CHROM == chr_name) %>%
    mutate(window = floor(POS / window_size) * window_size) %>%
    group_by(CHROM, window) %>%
    summarize(mean_value = mean(M_F_HetDif_allo_LA_Ecaeruleum, na.rm = TRUE))  # Summarize with NA handling
  
  # Save the plot for Allo_Het_Sex_Dif
  pdf(paste0("EcaerAlloLA_HetDif_bySite_Chr", chr_num, ".pdf"), width = 11, height = 5)
  ggplot(EcaerAlloLA_HetDif_binned, aes(x = window, y = mean_value)) +
    geom_point() + scale_y_continuous(limits=c(-0.501,0.501)) + 
    labs(x = "Position (bp)", y = "Mean Value", title = paste("Ecaer LA Value across Chromosome", chr_num, "by site")) +
    theme_minimal()
  dev.off()
}