
library(dplyr)
library(readr)
library(ggplot2)
require(gridExtra)
library(readxl)
library(writexl)

bcot_freq_grand=read_excel("/Users/rlmoran/Desktop/darters/SDR/backcross_data/MAF_bcot_freq_grand.xlsx", sheet="MAF_bcot_freq_grand")

bcrb_freq_grand=read_excel("/Users/rlmoran/Desktop/darters/SDR/backcross_data/MAF_bcrb_freq_grand.xlsx", sheet="MAF_bcrb_freq_grand")

maf_data=read_excel("/Users/rlmoran/Desktop/darters/SDR/backcross_data/MAF_bcrb_bcot_freq_grand.xlsx", sheet="Sheet1")

maf_delta <- maf_data %>%
  select(LG, POS, MAF_OBS, MT) %>%
  tidyr::pivot_wider(names_from = MT, values_from = MAF_OBS, names_prefix = "MAF_") %>%
  mutate(delta_MAF = MAF_OT - MAF_RB)

# Step to create a genomic position index for a Manhattan-style plot
maf_delta <- maf_delta %>%
  group_by(LG) %>%
  arrange(POS, .by_group = TRUE) %>%
  mutate(SNP_index = row_number()) %>%
  ungroup()

# Get cumulative position for each LG to space chromosomes on x-axis
lg_offsets <- maf_delta %>%
  group_by(LG) %>%
  summarize(max_index = max(SNP_index)) %>%
  mutate(offset = lag(cumsum(max_index), default = 0))

# Join offset back to data
maf_delta <- maf_delta %>%
  left_join(lg_offsets, by = "LG") %>%
  mutate(cum_pos = SNP_index + offset)

# Compute mean ΔMAF and x-axis range per LG
mean_lines <- maf_delta %>%
  group_by(LG) %>%
  summarize(mean_delta = mean(delta_MAF, na.rm = TRUE),
            xmin = min(cum_pos),
            xmax = max(cum_pos))

# Create alternating black/gray colors for LGs
lg_list <- sort(unique(maf_delta$LG))
lg_colors <- rep(c("deepskyblue3", "deepskyblue4"), length.out = length(lg_list))
names(lg_colors) <- lg_list  # match LG values to colors

### for pub ###
# Plot with alternating gray/black LG colors
ggplot(maf_delta, aes(x = cum_pos, y = delta_MAF)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.2682, ymax = 0.2736),
            fill = "#ededed", alpha = 0.1, inherit.aes = FALSE) +
  geom_point(aes(color = as.factor(LG)), size = 0.5, alpha = 0.6) +
  geom_segment(data = mean_lines, aes(x = xmin, xend = xmax, y = mean_delta, yend = mean_delta),
               color = "black", size = 0.8) +
  scale_x_continuous(
    breaks = mean_lines$xmin + (mean_lines$xmax - mean_lines$xmin) / 2,
    labels = mean_lines$LG
  ) +
  scale_color_manual(values = lg_colors) +
  labs(
    x = "Chromosome",
    y = expression(Delta~"MAF (orangethroat - rainbow)"),
    title = expression("ΔMAF by Mitochondrial Haplotype Across the Genome")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(vjust = 0.5),
    legend.position = "none",
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),           # Removes grid lines
    panel.background = element_blank(),     # Removes plot background
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.title.y = element_text(),
    axis.text.y = element_text(),
    axis.line = element_line(color = "black")  # Restores x and y axis lines
  )
  
  
ggplot(maf_delta %>% filter(LG == "9"), aes(x = POS, y = delta_MAF)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.2682, ymax = 0.2736),
            fill = "#ededed", alpha = 0.1, inherit.aes = FALSE) +
  geom_smooth(method = "loess", span = 0.02, se = FALSE, color = "black") +
  geom_point(aes(color = as.factor(LG)), size = 0.5, alpha = 0.6) +
  scale_color_manual(values = lg_colors) +
  scale_y_continuous(limits = c(-0.50,0.50)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  labs(
    x = "bp",
    y = expression(Delta~"MAF (orangethroat - rainbow)"),
    title = expression("Chromosome 9")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0.5),
    legend.position = "none",
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),           # Removes grid lines
    panel.background = element_blank(),     # Removes plot background
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.title.y = element_text(),
    axis.text.y = element_text(),
    axis.line = element_line(color = "black"))  # Restores x and y axis lines
  

#600x250
ggplot(maf_delta %>% filter(LG == "23"), aes(x = POS, y = delta_MAF)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.2682, ymax = 0.2736),
            fill = "#ededed", alpha = 0.1, inherit.aes = FALSE) +
  geom_smooth(method = "loess", span = 0.01, se = FALSE, color = "black") +
  #geom_point(aes(color = as.factor(LG)), size = 0.5, alpha = 0.6) +
  scale_color_manual(values = lg_colors) +
  labs(
    x = "bp",
    y = expression(Delta~"MAF (orangethroat - rainbow)"),
    title = expression("Chromosome 23")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(vjust = 0.5),
    legend.position = "none",
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),           # Removes grid lines
    panel.background = element_blank(),     # Removes plot background
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.title.y = element_text(),
    axis.text.y = element_text(),
    axis.line = element_line(color = "black")  # Restores x and y axis lines
  )

#600x250



delta_maf_data <- read_excel("/Users/rlmoran/Desktop/darters/SDR/backcross_data/delta_MAF_by_LG.xlsx")

# Order chromosomes
delta_maf_data <- delta_maf_data %>%
  mutate(LG = factor(LG, levels = unique(LG)))

# Alternate colors
delta_maf_data <- delta_maf_data %>%
  mutate(color_group = as.factor(as.numeric(LG) %% 2))

# Plot
ggplot(delta_maf_data, aes(x = POS, y = delta_MAF, color = color_group)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.2682, ymax = 0.2736),
            fill = "lightgray", alpha = 0.3, inherit.aes = FALSE) +
  geom_point(size = 0.5) +
  geom_hline(data = delta_maf_data %>%
               group_by(LG) %>%
               summarize(mean_delta_MAF = mean(delta_MAF, na.rm = TRUE)),
             aes(yintercept = mean_delta_MAF, group = LG), 
             color = "darkred", linetype = "dashed") +
  scale_color_manual(values = c("black", "gray")) +
  facet_wrap(~LG, scales = "free_x", nrow = 1) +
  labs(
    x = "Position",
    y = expression(Delta~MAF~(OT[MT] - RB[MT])),
    title = expression(Delta~MAF~Across~Genome)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2, "lines"),
    strip.text = element_text(size = 8)
  )


maf_wide <- maf_data %>%
  select(LG, POS, MAF_OBS, MT) %>%
  tidyr::pivot_wider(names_from = MT, values_from = MAF_OBS, names_prefix = "MAF_") %>%
  mutate(delta_MAF = MAF_OT - MAF_RB)

write_xlsx(maf_wide, "/Users/rlmoran/Desktop/darters/SDR/backcross_data/delta_MAF_by_LG.xlsx")



# Prepare the data for cumulative positions across chromosomes
bcrb_freq_grand <- bcrb_freq_grand %>%
  mutate(LG = as.factor(LG)) %>%
  group_by(LG) %>%
  arrange(POS) %>%
  mutate(chr_len = max(POS)) %>%
  ungroup() %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  mutate(BPcum = POS + tot[as.numeric(LG)]) 

# Get chromosome center positions for x-axis labels
axis_df <- bcrb_freq_grand %>%
  group_by(LG) %>%
  summarize(center = (min(BPcum) + max(BPcum)) / 2)

MAFdev = abs((bcrb_freq_grand$MAF_OBS)-0.25)
bcrb_freq_grand$MAFdev = MAFdev
# Now plot
gemmachr9=ggplot(bcrb_freq_grand, aes(x = BPcum, y = MAFdev)) +
  geom_point(aes(color = as.factor(LG)), alpha = 0.75, size = 0.9) +
  scale_color_manual(values = rep(c("gray10", "gray50"), 12)) +
  scale_x_continuous(label = axis_df$chr, breaks = axis_df$center) +
  scale_y_continuous(limits = c(0, 20)) +
  geom_hline(yintercept = -log10(5e-08), color = "red", linetype = "dashed") +
  labs(x = "Chromosome", y = "-log10(p-value)") +
  theme_minimal() + geom_vline(xintercept = 12051184, color = "red", linetype = "dashed", size = 1) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


FDR_P <- p.adjust(bcot_freq_grand$chisq.p, method="fdr")
sort_FDR = sort(FDR_P)
head(sort_FDR)
bcot_freq_grand$FDR_P = FDR_P

MAFdev = abs((bcot_freq_grand$MAF_OBS)-0.25)
bcot_freq_grand$MAFdev = MAFdev

# Line plot from twisst weights
p1 <- ggplot(df, aes(x = position, y = weight)) +
  geom_line(color = "black", size = 0.5) +
  coord_cartesian(xlim = c(0,	7000000)) +
  labs(x = NULL, y = "Weighting") +
  theme_minimal()

# Manhattan plot using ggplot
p2 <- ggplot(GEMMA_RB_allo_dat, aes(x = ps, y = -log10(p_wald))) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  coord_cartesian(xlim = c(0,	7000000)) +
  labs(x = "Position", y = "-log10(p-value)") +
  theme_minimal()


maf1 <- ggplot(bcot_freq_grand %>% filter(LG == "23"), aes(POS, MAF_OBS)) +
  geom_point(aes(color = FDR_P < 0.05), size = 0.5) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), name = "FDR < 0.05") +
  scale_y_continuous(limits = c(0, 0.51)) +
  coord_cartesian(xlim = c(0, 7000000)) +
  geom_hline(yintercept = 0.25, color = "gray", linetype = "dashed") +
  theme_minimal() + theme(legend.position = "none")

# Stack the plots
library(patchwork)
#p1 / p2


grid.arrange(p1,p2,maf1, nrow=3)


##################

# Line plot from twisst weights
p1 <- ggplot(df, aes(x = position, y = weight)) +
  geom_line(color = "black", size = 0.5) +
  coord_cartesian(xlim = c(0,	24000000)) +
  labs(x = NULL, y = "Weighting") +
  theme_minimal()


#GEMMA_RB<-read.table("chr23_processed_rb_pooled_assoc.txt",header=T)
#SNP <- c(1:(nrow(GEMMA_RB)))
#GEMMA_RB_allo_dat <- data.frame(SNP,GEMMA_RB)
# Manhattan plot using ggplot
p2 <- ggplot(GEMMA_RB_allo_dat, aes(x = ps, y = -log10(p_wald))) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  coord_cartesian(xlim = c(0,	24000000)) + 
  geom_vline(xintercept = 4519966, color = "green", linetype = "dashed", alpha = 0.8, size = 1) + #atp23
  
  labs(x = "Position", y = "-log10(p-value)") +
  theme_minimal()

chr23_MAF_subset=bcot_freq_grand %>% filter(LG == "23" & FDR_P < 0.059)

maf1 <- ggplot(bcot_freq_grand %>% filter(LG == "23"), aes(POS, MAF_OBS)) +
  geom_point(aes(color = FDR_P < 0.05), size = 0.5) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), name = "FDR < 0.05") +
  scale_y_continuous(limits = c(0, 0.51)) +
  coord_cartesian(xlim = c(0, 24000000)) +
  geom_hline(yintercept = 0.25, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 4519966, color = "green", linetype = "dashed", alpha = 0.8, size = 1) + #atp23
  geom_vline(xintercept = 12486820, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #chchd3a
  geom_vline(xintercept = 21618507, color = "purple", linetype = "dashed", alpha = 0.3, size = 1) + #ndufa5
  geom_vline(xintercept = 15641759.5, color = "purple", linetype = "dashed", alpha = 0.3, size = 1) + #ndufb2
  geom_vline(xintercept = 17601849.5, color = "purple", linetype = "dashed", alpha = 0.3, size = 1) + #ndufa12
  geom_vline(xintercept = 17780862, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #ssbp1
  geom_vline(xintercept = 10797990.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #mrpl42
  geom_vline(xintercept = 15833579.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #tigarb
  geom_vline(xintercept = 17719082.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #tmem243b
  geom_vline(xintercept = 2367628, color = "green", linetype = "dashed", alpha = 0.3, size = 1) + #apex1
  geom_vline(xintercept = 12775864, color = "orange", linetype = "dashed", alpha = 0.3, size = 1) + #lncRNA 
  geom_vline(xintercept = 4817904, color = "pink", linetype = "dashed", alpha = 0.3, size = 1) + #sox5
  geom_vline(xintercept = 12618509, color = "pink", linetype = "dashed", alpha = 0.3, size = 1) + #exoc4 
  geom_vline(xintercept = 1848888, color = "pink", linetype = "dashed", alpha = 0.3, size = 1) + #fish egg lectin-like LOC116673457 
  geom_vline(xintercept = 18183937, color = "yellow", linetype = "dashed", alpha = 0.5, size = 1) + #atp5f1c
  theme_minimal() + theme(legend.position = "none")



mean_dxy_10kb_chr23=ggplot(Dxy_dat_binned %>% filter(scaffold == "NC_045755.1"), aes(x = base, y = mean_dxy)) +
  geom_point(color = "black", alpha = 0.75, size = 0.8) +   coord_cartesian(xlim = c(0, 24000000)) +
  geom_vline(xintercept = 18183937, color = "pink", linetype = "dashed", alpha = 0.3, size = 1) + #atp5f1c
  theme_minimal() + theme(legend.position = "none")
  #geom_vline(xintercept = 12776198, color = "red", linetype = "dashed", alpha = 0.6, size = 1) + #lncRNA - gnrh1 and amhr2 like??







mean_dxy_10kb_chr23

grid.arrange(p1,p2,maf1,mean_dxy_10kb_chr23 , nrow=4)

#1600 x 600



###################
bin_size <- 100000  # 1 kb bins

# Add a bin column to each dataset
bcot_chr23 <- bcot_freq_grand %>%
  filter(LG == "23") %>%
  mutate(bin = floor(POS / bin_size))

GEMMA_chr23 <- GEMMA_RB_allo_dat %>%
  filter(chr == 23) %>%
  rename(POS = ps) %>%
  mutate(bin = floor(POS / bin_size))

# Average MAF per bin
bcot_summary <- bcot_chr23 %>%
  group_by(bin) %>%
  summarize(mean_maf = mean(chisq.p, na.rm = TRUE), .groups = "drop")

# Minimum p-value per bin (or mean/log-mean if preferred)
gemma_summary <- GEMMA_chr23 %>%
  group_by(bin) %>%
  summarize(min_p = min(p_wald, na.rm = TRUE), .groups = "drop") %>%
  mutate(log_p = -log10(min_p))

merged_bins <- inner_join(bcot_summary, gemma_summary, by = "bin")

cor.test(merged_bins$mean_maf, merged_bins$log_p, method = "pearson")

ggplot(merged_bins, aes(x = mean_maf, y = log_p)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "Mean MAF (per bin)", y = "-log10(min p-value per bin)",
       title = "Correlation between binned MAF and GWAS signal (Chr 23)") +
  theme_minimal()


################## chr 9

# Line plot from twisst weights
twisst_chr9 <- ggplot(df_chr9, aes(x = position, y = weight)) +
  geom_line(color = "black", size = 0.5) +
  coord_cartesian(xlim = c(0,	39000000)) +
  labs(x = NULL, y = "Weighting") +
  theme_minimal()




Espec_all<-read.table("/Users/rlmoran/Desktop/darters/SDR/no_unplaced_processed_ot_assoc.txt",header=T)

SNP <- c(1:(nrow(Espec_all)))
Espec_dat <- data.frame(SNP,Espec_all)

# Manhattan plot using ggplot
gemma_chr9 <- ggplot(Espec_dat, aes(x = ps, y = -log10(p_wald))) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  coord_cartesian(xlim = c(0,	39000000)) + 
  #geom_vline(xintercept = 15364664, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #ndufs7
  geom_vline(xintercept = 37376362, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #atp5f1d
  labs(x = "Position", y = "-log10(p-value)") +
  theme_minimal()

chr9_bcot_FDR=bcot_freq_grand %>% filter(LG == "9")
write.csv(chr9_bcot_FDR, "/Users/rlmoran/Desktop/darters/SDR/backcross_data/chr9_bcot_FDR.csv", row.names = FALSE)



SNP <- c(1:(nrow(delta_maf_data)))
delta_maf_data <- data.frame(SNP,delta_maf_data)


FDR_P <- p.adjust(delta_maf_data$Pvalue, method = "BH")
delta_maf_data$FDR_P = FDR_P


maf_all_bc<- ggplot(delta_maf_data%>% filter(LG == "23"), aes(SNP, -log10(Pvalue))) +
  geom_point(size = 0.5, alpha = 0.6) +
  facet_wrap(~LG, scales = "free_x", nrow = 1) +
  scale_x_continuous(
    breaks = mean_lines$xmin + (mean_lines$xmax - mean_lines$xmin) / 2,
    labels = mean_lines$LG
  ) +
  labs(
    x = "Position",
    y = "-log10(p)",
    title = "Mito-nuclear association backcrosses"
  ) + theme_minimal()+ theme(legend.position = "none")


ggplot(maf_delta, aes(x = cum_pos, y = delta_MAF)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.2682, ymax = 0.2736),
            fill = "#ededed", alpha = 0.1, inherit.aes = FALSE) +
  geom_point(aes(color = as.factor(LG)), size = 0.5, alpha = 0.6) +
  geom_segment(data = mean_lines, aes(x = xmin, xend = xmax, y = mean_delta, yend = mean_delta),
               color = "black", size = 0.8) +
  scale_x_continuous(
    breaks = mean_lines$xmin + (mean_lines$xmax - mean_lines$xmin) / 2,
    labels = mean_lines$LG
  ) +
  scale_color_manual(values = lg_colors) +
  labs(
    x = "Chromosome",
    y = expression(Delta~"MAF (orangethroat - rainbow)"),
    title = expression("ΔMAF by Mitochondrial Haplotype Across the Genome")
  ) +

  geom_point(aes(color = FDR_P < 0.059), size = 0.5) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), name = "FDR < 0.05") +
  scale_y_continuous(limits = c(0, 0.51)) +
  coord_cartesian(xlim = c(0, 39000000)) +
  geom_hline(yintercept = 0.25, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 15364664, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #ndufs7
  geom_vline(xintercept = 15719195, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #ndufa7
  geom_vline(xintercept = 30182131.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #lonp1
  geom_vline(xintercept = 37195911.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #polrmt
  geom_vline(xintercept = 20425783.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #mfn1
  geom_vline(xintercept = 38169326.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #gpx4
  geom_vline(xintercept = 18804830.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #opa1
  geom_vline(xintercept = 20476883, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #mrpl47
  geom_vline(xintercept = 37082900, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #micos13
  geom_vline(xintercept = 683500, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #cpt2
  geom_vline(xintercept = 20463981.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #aldh7a1
  geom_vline(xintercept = 20480532, color = "pink", linetype = "dashed", alpha = 0.7, size = 1) + #ndufb5
  geom_vline(xintercept = 30415936, color = "pink", linetype = "dashed", alpha = 0.7, size = 1) + #ndufa11
  geom_vline(xintercept = 36348608, color = "orange", linetype = "dashed", alpha = 0.3, size = 1) + #amh
  theme_minimal() + theme(legend.position = "none")





mean_dxy_10kb_chr9=ggplot(Dxy_dat_binned %>% filter(scaffold == "NC_045741.1"), aes(x = base, y = mean_dxy)) +
  geom_point(color = "black", alpha = 0.75, size = 0.8) +   coord_cartesian(xlim = c(0, 39000000)) +
  geom_vline(xintercept = 15364664, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #ndufs7
  geom_vline(xintercept = 15719195, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #ndufa7
  geom_vline(xintercept = 30182131.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #lonp1
  geom_vline(xintercept = 37195911.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #polrmt
  geom_vline(xintercept = 20425783.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #mfn1
  geom_vline(xintercept = 38169326.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #gpx4
  geom_vline(xintercept = 18804830.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #opa1
  geom_vline(xintercept = 20476883, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #mrpl47
  geom_vline(xintercept = 37082900, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #micos13
  geom_vline(xintercept = 683500, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #cpt2
  geom_vline(xintercept = 20463981.5, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #aldh7a1
  geom_vline(xintercept = 20480532, color = "blue", linetype = "dashed", alpha = 0.3, size = 1) + #ndufb5
  geom_vline(xintercept = 36348608, color = "orange", linetype = "dashed", alpha = 0.3, size = 1) + #amh
  theme_minimal() + theme(legend.position = "none")

36340165-36348608
lonp1
gpx4
mrpl47
ndufs7
polrmt
cpt2
mfn1
micos13
aldh7a1
opa1
ndufb5

30182131.5
37195911.5
20425783.5
38169326.5
18804830.5
20476883
37082900
683500
20463981.5
20480532


grid.arrange(twisst_chr9,gemma_chr9, maf1_chr9,mean_dxy_10kb_chr9, nrow=4)
#1500x600

