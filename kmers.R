library(tidyverse)
library(readr)
library(dplyr)

df <- read_tsv("Ykmer_ec_density.tsv", col_names = FALSE)
map <- read_tsv("scaffold_to_chromosome.tsv",col_types = cols()) %>%
  rename(chr_num = chr)
colnames(df) <- c("chr","start","end","n_hits","bases_cov","win_len","frac_cov")
subset_df_base <- subset(df, grepl("^NC", chr))

out <- subset_df_base %>%
  left_join(map, by = c("chr" = "scaffold")) %>%
  # warn if anything didn't map
  { 
    missing <- dplyr::filter(., is.na(chr_num)) %>% distinct(chr)
    if (nrow(missing) > 0) warning("Unmapped scaffolds: ", paste(missing$chr, collapse = ", "))
    .
  } %>%
  mutate(chr = chr_num) %>%              # replace scaffold IDs with numeric chr
  select(-chr_num)
head(out)

### 

out <- out %>%
  mutate(chr = as.integer(chr)) %>%
  arrange(chr, start)

chr_offsets <- out %>%
  group_by(chr) %>%
  summarise(chr_len = max(end)) %>%
  arrange(chr) %>%
  mutate(offset = lag(cumsum(chr_len), default = 0))

out2 <- out %>%
  mutate(mid = (start + end) / 2) %>%
  left_join(chr_offsets, by = "chr") %>%
  mutate(genome_pos = mid + offset)

# chromosome centers for labeling (needs one row per chr)
chr_labels <- chr_offsets %>%
  distinct(chr, chr_len, offset) %>%
  mutate(center = offset + chr_len/2) %>%
  arrange(as.integer(chr))

p_genome <- ggplot(out2, aes(x = genome_pos, y = n_hits)) +
  geom_line(color="darkviolet") +
  geom_vline(
    data = chr_labels,
    aes(xintercept = offset),
    linetype = "dashed", linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = chr_labels$center,
    labels = paste0(chr_labels$chr)
  ) +
  labs(title="E. caeruleum genome-wide Y-mer density", x = "Chromosome", y = "Y-kmer hits per 50 kb window") +
  theme_classic() +
  theme(axis.ticks.x = element_blank())

p_genome

chr23_df <- subset(out2, chr == 9)
chr9_2 <- chr23_df %>%
  mutate(mid = (start + end) / 2) %>%
  left_join(chr_offsets, by = "chr") %>%
  mutate(genome_pos = mid)

plot <- ggplot(chr9_2,aes(x=genome_pos,y=n_hits))+
  geom_line(color="darkorange") +   labs(title = "Chromosome 9 Y-mer density for E. spectabile",x = "Chromosome", y = "Y-kmer hits per 50 kb window") +
  theme_classic() +
  theme(axis.ticks.x = element_blank())

plot
### attempting to find if windows in high density are in amh ###
plot + 
  geom_vline(xintercept = 36340165, linetype = "dashed", linewidth = 0.6) + 
  geom_vline(xintercept = 36348608, linetype="dashed",linewidth = 1, col='green')
 
chr_zoom <- plot +  coord_cartesian(xlim = c(36300000, 36350000))
chr_zoom

### bar plot method of showing sex chrom peaks ###
kmer_chr_sum <- out2 %>%          # or use your df object
  mutate(CHROM = as.character(chr)) %>%     # match your naming convention
  group_by(CHROM) %>%
  summarise(
    mean = mean(n_hits, na.rm = TRUE),
    se   = sd(n_hits, na.rm = TRUE) / sqrt(sum(!is.na(n_hits))),
    n_windows = sum(!is.na(n_hits)),
    .groups = "drop"
  ) %>%
  mutate(
    CHROM = factor(CHROM, levels = as.character(1:24))  # ensure all 24 show up
  )
kmer_chr_sum_clean <- kmer_chr_sum %>%
  dplyr::filter(!is.na(CHROM))
ggplot(kmer_chr_sum_clean, aes(x = CHROM, y = mean)) +
  geom_col(fill = "darkviolet", width = 0.9) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "E. caeruleum Y-kmer density — mean per chromosome",
    y = "Mean kmer hits per 50 kb window",
    x = "Chromosome"
  ) + theme_minimal()
