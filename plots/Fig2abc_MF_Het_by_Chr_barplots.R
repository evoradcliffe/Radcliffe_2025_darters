# Load libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(tidyverse)


setwd("/Users/rlmoran/Desktop/darters/SDR/heterozygosity/RB2RB/by_site_filtered")

### ggplot theme

MyTheme=theme_bw() +
  theme(strip.background = element_blank()) +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(hjust = 0.5,family="Arial", size=16, color="black", face = "bold"),
        panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
        panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
        #panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
        axis.line = element_line(colour = 'black', size = 0.5),#sets the axis line size
        axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
        axis.title.x = element_text(family="Arial", size=20, color="black", face = "bold"), #size of x-axis title
        axis.title.y = element_text(family="Arial", size=20, color="black"), #size of x-axis title
        axis.text.x = element_text(family="Arial", size=18, color="black"), #size of x-axis text and angle of labels (vjust = 0.5)
        axis.text.y = element_text(family="Arial", size=18, color="black"))#size of y-axis text

###


#read in OT data (mean & SE het dif by chr)

summary_OT=read.table("chrom_diff_summary_OT2OT_pooled.txt", header=T)


# Pivot longer
OT_long <- summary_OT %>%
  pivot_longer(
    cols = -CHROM,
    names_to = c("measure", "group"),
    names_sep = "_diff_"
  ) %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  )


#make sure all 24 chrs are listed as factors (so they all show up on the plot below)
OT_long$CHROM <- factor(OT_long$CHROM, levels = as.character(1:24))


#make the barplot
ggplot(OT_long, aes(x = CHROM, y = mean)) +
  geom_col(fill = "lightgray", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.9),
                width = 0.2) + ylim(-0.0001,0.0008) +
  labs(title = "Orangethroat darter - Mean M:F Heterozygosity Differences",
       y = "Mean M:F Heterozygosity",
       x = "Chromosome")  + MyTheme







#read in RB data (mean & SE het dif by chr)

summary_RB=read.table("chrom_diff_summary_RB2RB_pooled.txt", header=T)


# Pivot longer
RB_long <- summary_RB %>%
  pivot_longer(
    cols = -CHROM,
    names_to = c("measure", "group"),
    names_sep = "_diff_"
  ) %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  )


#make sure all 24 chrs are listed as factors (so they all show up on the plot below)
RB_long$CHROM <- factor(RB_long$CHROM, levels = as.character(1:24))

#make the barplot
ggplot(RB_long, aes(x = CHROM, y = mean)) +
  geom_col(fill = "lightgray", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.9),
                width = 0.2) + ylim(-0.0002,0.005) +
  labs(title = "Orangethroat darter - Mean M:F Heterozygosity Differences",
       y = "Mean M:F Heterozygosity",
       x = "Chromosome")  + MyTheme






