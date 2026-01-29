library(readxl)
library(dplyr)
library(purrr)

setwd("~/Documents/Manuscripts/NatComms_SDR")


## 1. Read mito-nuclear gene coordinates (orthologs on OT genome)
mt <- read_excel("MTnucgenes_OTcoords.xlsx") %>%
  transmute(
    Chr   = as.integer(Chromosomes),
    start = `Annotation Genomic Range Start`,
    stop  = `Annotation Genomic Range Stop`,
    Symbol
  )



## 2. Read TRD results for both backcross directions
bc_ot <- read.csv("MAF_bcot_freq_grand.csv")  # OT-mom backcrosses
bc_rb <- read.csv("MAF_bcrb_freq_grand.csv")  # RB-mom backcrosses

# harmonize column names + chromosome field
bc_ot <- bc_ot %>%
  transmute(
    Chr = Chr, POS,
    MAF_OT = MAF_OBS,
    N_OT   = NCHROBS,
    kOT    = OBS_MIN,
    pTRD_OT = chisq.p
  )

bc_rb <- bc_rb %>%
  transmute(
    Chr = Chr, POS,
    MAF_RB = MAF_OBS,
    N_RB   = NCHROBS,
    kRB    = OBS_MIN,
    pTRD_RB = chisq.p
  )

# merge and compute summary stats
dat <- full_join(bc_ot, bc_rb, by = c("Chr","POS")) %>%
  mutate(
    sex_chr    = Chr %in% c(9, 23),
    
    # distortion (TRD magnitude) summarized across directions
    dist_score = (abs(MAF_OT - 0.25) + abs(MAF_RB - 0.25)) / 2,
    
    # cytoplasmic-background allele freq shift
    deltaMAF   = MAF_OT - MAF_RB
  )

# TRD significance (within-direction)
dat <- dat %>%
  mutate(
    trd_sig_OT  = p.adjust(pTRD_OT, method="BH") < 0.05,
    trd_sig_RB  = p.adjust(pTRD_RB, method="BH") < 0.05,
    trd_sig_any = trd_sig_OT | trd_sig_RB
  )

max(abs(dat$kOT - round(dat$kOT)), na.rm=TRUE)
max(abs(dat$kRB - round(dat$kRB)), na.rm=TRUE)

dat <- dat %>%
  mutate(
    kOT_i  = pmin(pmax(as.integer(round(kOT)), 0L), as.integer(round(N_OT))),
    kRB_i  = pmin(pmax(as.integer(round(kRB)), 0L), as.integer(round(N_RB))),
    N_OT_i = as.integer(round(N_OT)),
    N_RB_i = as.integer(round(N_RB))
  )

dat$cyto_p <- mapply(function(k1,n1,k2,n2){
  if (any(is.na(c(k1,n1,k2,n2))) || n1==0 || n2==0) return(NA_real_)
  mat <- matrix(c(k1, n1-k1, k2, n2-k2), nrow=2, byrow=TRUE)
  
  # Fisher if any expected small
  if (any(mat < 5)) {
    fisher.test(mat)$p.value
  } else {
    suppressWarnings(chisq.test(mat, correct = FALSE)$p.value)
  }
}, dat$kOT_i, dat$N_OT_i, dat$kRB_i, dat$N_RB_i)

dat$cyto_assoc <- p.adjust(dat$cyto_p, "BH") < 0.05



##Ask whether Chr9/23 have unusually high delta MAF
set.seed(1)
nperm <- 10000

auto <- dat %>% filter(!sex_chr, is.finite(deltaMAF))
chr_test <- function(chr_id) {
  chrdat <- dat %>% filter(Chr == chr_id, is.finite(deltaMAF))
  obs_n <- nrow(chrdat)
  obs_stat <- mean(abs(chrdat$deltaMAF))
  
  perm_stat <- replicate(nperm, {
    samp <- sample(auto$deltaMAF, size = obs_n, replace = FALSE)
    mean(abs(samp))
  })
  
  p <- (1 + sum(perm_stat >= obs_stat)) / (nperm + 1)
  list(chr=chr_id, obs=obs_stat, p=p)
}

chr_test(9)
chr_test(23)


##
dat2 <- dat %>% filter(is.finite(deltaMAF), is.finite(dist_score))
lm_abs <- lm(abs(deltaMAF) ~ sex_chr + dist_score, data=dat2)
summary(lm_abs)

dat2 <- dat %>% filter(is.finite(deltaMAF))
dat2 %>% group_by(sex_chr) %>% summarise(n=n(), mean_abs=mean(abs(deltaMAF)), med_abs=median(abs(deltaMAF)))


#Preserve distorion - stratify by distortion bins
nbins <- 10
auto2 <- auto %>% mutate(bin = cut(dist_score,
                                   breaks = quantile(dist_score, probs=seq(0,1,length.out=nbins+1),
                                                     na.rm=TRUE),
                                   include.lowest=TRUE))
chrdat <- dat %>% filter(Chr == 23, is.finite(deltaMAF)) %>%
  mutate(bin = cut(dist_score,
                   breaks = quantile(auto2$dist_score, probs=seq(0,1,length.out=nbins+1),
                                     na.rm=TRUE),
                   include.lowest=TRUE))
bin_counts <- table(chrdat$bin)

obs_stat <- mean(abs(chrdat$deltaMAF))

perm_stat <- replicate(nperm, {
  vals <- numeric(0)
  for (b in names(bin_counts)) {
    m <- as.integer(bin_counts[[b]])
    pool <- auto2$deltaMAF[auto2$bin == b]
    vals <- c(vals, sample(pool, size=m, replace=(length(pool)<m)))
  }
  mean(abs(vals))
})

p <- (1 + sum(perm_stat >= obs_stat)) / (nperm + 1)
p












#####
####
## 50 kb window analyses NDUF genes

library(readxl)
library(dplyr)
library(data.table)



## 1. Filter mito-nuclear table to complex I NDUF* genes
nduf <- mt %>%
  filter(grepl("^nduf", Symbol, ignore.case = TRUE)) %>%
  mutate(
    Chr   = as.integer(Chr),
    start = as.numeric(start),
    stop  = as.numeric(stop)
  )

nrow(nduf)  # sanity check: how many NDUF genes you have

## 2. Define ±50 kb windows around NDUF genes
window_size <- 50000

nduf_windows <- nduf %>%
  dplyr::mutate(
    win_start = pmax(start - window_size, 1),
    win_end   = stop + window_size
  ) %>%
  dplyr::select(Chr, win_start, win_end, Symbol)

## 3. Overlap SNPs with NDUF windows
dt_snps <- as.data.table(dat) %>%
  mutate(Chr = as.integer(Chr),
         POS = as.numeric(POS))

setkey(dt_snps, Chr, POS)

dt_nduf_win <- as.data.table(nduf_windows) %>%
  rename(start = win_start, end = win_end)
setkey(dt_nduf_win, Chr, start, end)

# Make SNPs into ranges [POS, POS]
dt_snps_range <- copy(dt_snps)
dt_snps_range[, `:=`(start = POS, end = POS)]
setkey(dt_snps_range, Chr, start, end)

# Find overlaps
nduf_overlaps <- foverlaps(
  dt_snps_range,
  dt_nduf_win,
  by.x = c("Chr", "start", "end"),
  by.y = c("Chr", "start", "end"),
  nomatch = 0L
)

nduf_snps <- unique(nduf_overlaps[, .(Chr, POS)])

# Flag NDUF-window SNPs
dt_snps_range[, nduf_window := FALSE]
dt_snps_range[nduf_snps, on = .(Chr, POS), nduf_window := TRUE]

dat <- as.data.frame(dt_snps_range)

###
###

## Sex chromosomes only
sex_dat <- dat %>%
  filter(sex_chr, !is.na(dist_score), !is.na(cyto_assoc), !is.na(nduf_window))

# Distortion within vs outside NDUF windows
wilcox.test(dist_score ~ nduf_window, data = sex_dat)

# Proportion mito-associated within vs outside NDUF windows
tab_sex_nduf <- table(sex_dat$nduf_window, sex_dat$cyto_assoc)
tab_sex_nduf
fisher.test(tab_sex_nduf)

## Autosomes only
auto_dat <- dat %>%
  filter(!sex_chr, !is.na(dist_score), !is.na(cyto_assoc), !is.na(nduf_window))

wilcox.test(dist_score ~ nduf_window, data = auto_dat)

tab_auto_nduf <- table(auto_dat$nduf_window, auto_dat$cyto_assoc)
tab_auto_nduf
fisher.test(tab_auto_nduf)






glm_nduf <- glm(
  cyto_assoc ~ sex_chr * nduf_window + dist_score,
  family = binomial,
  data = dat %>% filter(!is.na(cyto_assoc), !is.na(dist_score), !is.na(nduf_window))
)
summary(glm_nduf)


####
##make permutation figures
####

# Mean absolute deltaMAF for a chromosome
obs_mean_abs_delta <- function(dat, chr_id) {
  dat %>%
    filter(Chr == chr_id, is.finite(deltaMAF)) %>%
    summarise(m = mean(abs(deltaMAF)), n = n()) %>%
    as.list()
}

# Unstratified permutation null: sample same n from autosomes
perm_null_unstrat <- function(dat, chr_id, nperm = 10000, seed = 1) {
  set.seed(seed)
  
  chr_info <- obs_mean_abs_delta(dat, chr_id)
  obs_stat <- chr_info$m
  obs_n    <- chr_info$n
  
  auto <- dat %>%
    filter(!sex_chr, is.finite(deltaMAF)) %>%
    pull(deltaMAF)
  
  perm_stat <- replicate(nperm, {
    samp <- sample(auto, size = obs_n, replace = FALSE)
    mean(abs(samp))
  })
  
  p <- (1 + sum(perm_stat >= obs_stat)) / (nperm + 1)
  
  list(chr = chr_id, obs = obs_stat, n = obs_n, perm = perm_stat, p = p)
}

# Stratified permutation null: match distortion bins (DistScore) of the focal chr
perm_null_stratified <- function(dat, chr_id, nperm = 10000, nbins = 10, seed = 1) {
  set.seed(seed)
  
  auto2 <- dat %>%
    filter(!sex_chr, is.finite(deltaMAF), is.finite(dist_score)) %>%
    mutate(bin = cut(
      dist_score,
      breaks = quantile(dist_score, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE),
      include.lowest = TRUE
    ))
  
  chrdat <- dat %>%
    filter(Chr == chr_id, is.finite(deltaMAF), is.finite(dist_score)) %>%
    mutate(bin = cut(
      dist_score,
      breaks = quantile(auto2$dist_score, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE),
      include.lowest = TRUE
    ))
  
  obs_stat <- mean(abs(chrdat$deltaMAF))
  bin_counts <- table(chrdat$bin)
  
  perm_stat <- replicate(nperm, {
    vals <- numeric(0)
    for (b in names(bin_counts)) {
      m <- as.integer(bin_counts[[b]])
      pool <- auto2$deltaMAF[auto2$bin == b]
      vals <- c(vals, sample(pool, size = m, replace = (length(pool) < m)))
    }
    mean(abs(vals))
  })
  
  p <- (1 + sum(perm_stat >= obs_stat)) / (nperm + 1)
  
  list(chr = chr_id, obs = obs_stat, n = nrow(chrdat), perm = perm_stat, p = p,
       nbins = nbins)
}

#plot null & obs

plot_perm_hist <- function(res, title) {
  df <- data.frame(perm = res$perm)
  
  ggplot(df, aes(x = perm)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = res$obs, linewidth = 1.1) +
    labs(
      title = title,
      x = "Permutation mean(|ΔMAF|)",
      y = "Count"
    ) +
    annotate(
      "text",
      x = res$obs, y = Inf,
      label = paste0("obs = ", signif(res$obs, 4),
                     "\np = ", format(res$p, scientific = TRUE)),
      vjust = 1.2, hjust = -0.05
    ) +
    theme_classic()
}

# Run permutations
res9  <- perm_null_unstrat(dat, chr_id = 9,  nperm = 10000, seed = 1)
res23 <- perm_null_unstrat(dat, chr_id = 23, nperm = 10000, seed = 1)


res9_strat  <- perm_null_stratified(dat, chr_id = 9,  nperm = 10000, nbins = 10, seed = 1)
res23_strat <- perm_null_stratified(dat, chr_id = 23, nperm = 10000, nbins = 10, seed = 1)

# Plots
xlim_unstrat <- range(
  c(res9$perm,  res23$perm,  res9$obs,  res23$obs),
  na.rm = TRUE
)

xlim_strat <- range(
  c(res9_strat$perm, res23_strat$perm, res9_strat$obs, res23_strat$obs),
  na.rm = TRUE
)

# --- titles with statistic spelled out ---
ttlA <- "A. Chr 9: autosomal resampling null\nNull distribution of mean(|ΔMAF|)"
ttlB <- "B. Chr 23: autosomal resampling null\nNull distribution of mean(|ΔMAF|)"
ttlC <- "C. Chr 9: distortion-matched (binned) null\nNull distribution of mean(|ΔMAF|)"
ttlD <- "D. Chr 23: distortion-matched (binned) null\nNull distribution of mean(|ΔMAF|)"

# Plots (add consistent x limits)
pA <- plot_perm_hist(res9,  ttlA) + coord_cartesian(xlim = c(0.05,0.2))
pB <- plot_perm_hist(res23, ttlB) + coord_cartesian(xlim = c(0.05,0.2))

pC <- plot_perm_hist(res9_strat,  ttlC) + coord_cartesian(xlim = c(0.05,0.2))
pD <- plot_perm_hist(res23_strat, ttlD) + coord_cartesian(xlim = c(0.05,0.2))

figX <- (pA / pB / pC / pD)
figX

ggsave("SuppFigX_perm_nulls.pdf", figX, width = 7, height = 13)


