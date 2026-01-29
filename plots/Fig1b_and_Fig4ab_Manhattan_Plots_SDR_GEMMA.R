setwd("/home/rlmoran/moranlab/shared/darters/sdr_mysia/sex_association/output")

setwd("/Users/rlmoran/Desktop/darters/SDR")


library(qqman)
library(ggplot2)


Dxy_symEsEc=read.table("/Users/rlmoran/Desktop/darters/SDR/Dxy/symILcaeruleum_symILspectabile_Dxy_ChrOnly_noNAN_3cols_fixed.txt.gz", header = T)
sorted_Dxy_symEsEc<- Dxy_symEsEc %>%
  arrange(scaffold, base)

SNP <- c(1:(nrow(sorted_Dxy_symEsEc)))
Dxy_dat <- data.frame(SNP,sorted_Dxy_symEsEc)


Dxy_dat %>% filter(scaffold == "NC_045755.1")

mean_dxy_10kb_chr23=ggplot(Dxy_dat_binned %>% filter(scaffold == "NC_045755.1"), aes(x = base, y = mean_dxy)) +
  geom_point(aes(color = scaffold), alpha = 0.75, size = 0.8) + theme_minimal() + theme(legend.position = "none")

  
  
  
  scale_color_manual(values = rep(c("skyblue", "slateblue"), length(unique(Dxy_dat$scaffold)))) +
  scale_x_continuous(label = axis_df$scaffold, breaks = axis_df$center) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Scaffold", y = "dxy") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


##
Dxy_dat_binned <- Dxy_dat %>%
  mutate(window = floor(base / 1000) * 1000) %>%
  group_by(scaffold, window) %>%
  summarize(mean_dxy = mean(dxy, na.rm = TRUE),
            base = mean(base)) %>%
  ungroup()


scaffold_to_chr <- data.frame(
  scaffold = paste0("NC_0457", sprintf("%02d", 33:56), ".1"),
  chr = paste0("Chr_", 1:24),
  chr_num = 1:24
)

Dxy_dat_binned <- Dxy_dat_binned %>%
  left_join(scaffold_to_chr, by = "scaffold")


#manhattan_plot_Espec <- manhattan(Dxy_dat,chr="scaffold",bp="base",p="dxy",
                            #      snp="SNP",logp=TRUE,ylab="dxy",
                            #      suggestiveline = FALSE, ,ylim =c(0,1))

Dxy_dat <- Dxy_dat_binned %>%
  group_by(scaffold) %>%
  mutate(chr_len = max(base)) %>%
  ungroup() %>%
  mutate(scaffold = factor(scaffold, levels = unique(scaffold))) %>%
  group_by(scaffold) %>%
  summarise(chr_len = max(base)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  right_join(Dxy_dat, by = "scaffold") %>%
  mutate(pos_cum = base + tot)

axis_df <- Dxy_dat %>%
  group_by(scaffold) %>%
  summarize(center = (min(pos_cum) + max(pos_cum)) / 2)


pdf(paste0("Dxy_bySite_Chr.pdf"), width = 11, height = 5)

ggplot(Dxy_dat_binned, aes(x = pos_cum, y = dxy)) +
  geom_point(aes(color = scaffold), alpha = 0.75, size = 0.8) +
  scale_color_manual(values = rep(c("skyblue", "slateblue"), length(unique(Dxy_dat$scaffold)))) +
  scale_x_continuous(label = axis_df$scaffold, breaks = axis_df$center) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Scaffold", y = "dxy") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

##
Dxy_dat_binned <- Dxy_dat %>%
  mutate(window = floor(base / 1000) * 1000) %>%
  group_by(scaffold, window) %>%
  summarize(mean_dxy = mean(dxy, na.rm = TRUE),
            base = mean(base)) %>%
  ungroup()


scaffold_to_chr <- data.frame(
  scaffold = paste0("NC_0457", sprintf("%02d", 33:56), ".1"),
  chr = paste0("Chr_", 1:24),
  chr_num = 1:24
)

Dxy_dat_binned <- Dxy_dat_binned %>%
  left_join(scaffold_to_chr, by = "scaffold")

# Compute max base per chromosome
chr_offsets <- Dxy_dat_binned %>%
  group_by(chr_num) %>%
  summarise(chr_length = max(base, na.rm = TRUE)) %>%
  mutate(offset = lag(cumsum(chr_length), default = 0))

# Join offset and compute new position
Dxy_dat_binned <- Dxy_dat_binned %>%
  left_join(chr_offsets, by = "chr_num") %>%
  mutate(cum_pos = base + offset)

# Alternating colors
color_vec <- rep(c("skyblue3", "steelblue4"), length.out = 24)

axis_df <- Dxy_dat_binned %>%
  group_by(chr_num, chr) %>%
  summarise(center = (min(cum_pos) + max(cum_pos)) / 2)


ggplot(Dxy_dat_binned, aes(x = cum_pos, y = mean_dxy, color = as.factor(chr_num %% 2))) +
  geom_point(size = 0.5) +
  scale_color_manual(values = color_vec, guide = "none") +
  labs(x = "Genomic Position (across chromosomes)", y = "Mean Dxy", title = "Genome-wide Dxy") 
  + scale_x_continuous(
    label = axis_df$chr,
    breaks = axis_df$center
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




ggplot(Dxy_dat_binned, aes(x = base, y = mean_dxy)) +
  geom_line(size = 0.3, color = "steelblue") +
  facet_wrap(~ scaffold, scales = "free_x", ncol = 4) +
  labs(x = "Genomic Position", y = "Mean Dxy", title = "Genomic Landscape of Dxy by Chromosome") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




tiff(file="rainbow_LA_Rb.tiff",
     width=6, height=4, units="in", res=100, type="cairo")



no_unplaced_processed_ot_assoc.txt
no_unplaced_processed_rb_assoc.txt
## RB pooled ##
Ecaer_all<-read.table("no_unplaced_processed_rb_assoc.txt",header=T)

Ecaer_chr23 <- Ecaer_all %>%
  filter(chr == 23)


SNP <- c(1:(nrow(Ecaer_all)))
Ecaer_dat <- data.frame(SNP,Ecaer_all)

manhattan_plot_Ecaer <- manhattan(Ecaer_dat,chr="chr",bp="ps",p="p_wald",
                                  snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                                  suggestiveline = FALSE,   genomewideline = -log10(5e-08),ylim =c(0,15))


## OT pooled ##
Espec_all<-read.table("no_unplaced_processed_ot_assoc.txt",header=T)

SNP <- c(1:(nrow(Espec_all)))
Espec_dat <- data.frame(SNP,Espec_all)

manhattan_plot_Espec <- manhattan(Espec_dat,chr="chr",bp="ps",p="p_wald",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,   genomewideline = -log10(5e-08),ylim =c(0,15))

1100 x 400




## Epulch ##
Ep_all<-read.table("processed_Epulch.assoc.txt",header=T)

SNP <- c(1:(nrow(Ep_all)))
Ep_dat <- data.frame(SNP,Ep_all)

manhattan_plot_Ep <- manhattan(Ep_dat,chr="chr",bp="ps",p="p_wald",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,   genomewideline = -log10(5e-08))

Epulch_chr9 = Ep_dat %>% filter(chr == "9")
manhattan_plot_ob <- manhattan(Epulch_chr9,chr="chr",bp="ps",p="p_wald",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,   genomewideline = -log10(5e-08))


Epulch_chr23 = Ep_dat %>% filter(chr == "22")
manhattan_plot_ob <- manhattan(Epulch_chr23,chr="chr",bp="ps",p="p_wald",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,   genomewideline = -log10(5e-08))
#,

## OBs ##
ob_all<-read.table("processed_Orangebelly.assoc.txt",header=T)

SNP <- c(1:(nrow(ob_all)))
ob_dat <- data.frame(SNP,ob_all)

manhattan_plot_ob <- manhattan(ob_dat,chr="chr",bp="ps",p="p_wald",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,   genomewideline = -log10(5e-08),ylim =c(0,15))

ob_chr22 = ob_dat %>% filter(chr == "22")
manhattan_plot_ob <- manhattan(ob_chr22,chr="chr",bp="ps",p="p_wald",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,   genomewideline = -log10(5e-08),
                                 highlight =c(5504824:5504872))
manhattan_plot_ob_annotated = manhattan_plot_ob + annotate("text", x=14865758, y=9, label="ndufb4", color="red", size=5, hjust=0) 
  
ndufb4 = 14865758..14868522

########## patchwork 3 manhanttan plots ###########

library(patchwork)
class(manhattan_plot_Ecaer)
class(manhattan_plot_ob)
class(manhattan_plot_Espec)
# Suppress x-axis on first two plots
manhattan_plot_Ecaer_clean <- manhattan_plot_Ecaer +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

manhattan_plot_ob_clean <- manhattan_plot_ob +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

# Keep full x-axis for the bottom plot
manhattan_plot_Espec_clean <- manhattan_plot_Espec

# Stack vertically
combined_plot <- manhattan_plot_Ecaer_clean /
  manhattan_plot_ob_clean /
  manhattan_plot_Espec_clean

combined_plot
#######################################################

ob_chr9 = ob_dat %>% filter(chr == "9")
manhattan_plot_ob <- manhattan(ob_chr9,chr="chr",bp="ps",p="p_wald",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,   genomewideline = -log10(5e-08))

ob_chr23 = ob_dat %>% filter(chr == "23")
manhattan_plot_ob <- manhattan(ob_chr23,chr="chr",bp="ps",p="p_wald",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,   genomewideline = -log10(5e-08))

                    

## RBs LA## 
rb_all<-read.table("processed_allo_LA_rb.assoc.txt",header=T)

SNP <- c(1:(nrow(rb_all)))
rb_dat <- data.frame(SNP,rb_all)

manhattan_plot_rb <- manhattan(rb_dat,chr="chr",bp="ps",p="p_wald",
          snp="SNP",logp=TRUE,ylab="-log10(p-value)",
          suggestiveline = FALSE,   genomewideline = -log10(5e-08), ylim = c(0, 10))

#manhattan_plot_rb <- manhattan_plot_rb + 
#  annotate("text", x=10, y=5, label="My Label", color="red", size=5, hjust=0)

# Render the plot with the annotation
print(manhattan_plot_rb)

dev.off()

tiff(file="RB_LA_allo_only.tiff",
     width=8, height=4, units="in", res=100, type="cairo")



#### OT 2 OT pooled, chr 9 #####



ot_pooled_chr9<-read.table("Pooled_OT2OT_Chr9_assoc.txt",header=T)



SNP <- c(1:(nrow(ot_pooled_chr9)))
ot_dat_9 <- data.frame(SNP,ot_pooled_chr9)


manhattan_plot_ot <- manhattan(ot_dat_9,chr="chr",bp="ps",p="p_wald", 
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,
                               genomewideline = -log10(5e-08),
                               ylim = c(0, 20))


manhattan_plot_ot <- manhattan(ot_dat_9,chr="chr",bp="ps",p="p_wald", 
                              snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                              suggestiveline = FALSE,
                              genomewideline = -log10(5e-08),
                              ylim = c(0, 20) , xlim=c(12000000,39000000),
                              highlight =c(179089:179141, 113195:113203, 56093:56108))
dmrta2=113195:113203
amh=179089:179141
dmrt2b=12051184-12055409
manhattan_plot_ot <- manhattan_plot_ot +
  annotate("text", x=36340165, y=16, label="amh", color="red", size=5, hjust=0) +
  annotate("text", x=22382485, y=16, label="dmrta2", color="red", size=5, hjust=0) + 
  annotate("text", x=12051184, y=16, label="dmrt2b", color="red", size=5, hjust=0)


print(manhattan_plot_ot)

dev.off()



##### RB 2 RB pooled, chr 23 #####
chr23_only_RB2RB_assoc.txt




#ot_allo <- subset(ot_all, chr=="9")


SNP <- c(1:(nrow(ot_allo)))
ot_allo <- data.frame(SNP,ot_allo)


head(ot_allo)
manhattan_plot_ot <- manhattan(ot_allo,chr="chr",bp="ps",p="X1e.17",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,
                               genomewideline = -log10(5e-08))

,
                               highlight=c(8985:9005))
highlight =c(113190:113208,179080:179148))
dmrta2=113196:113203,
amh=179088:179140
manhattan_plot_ot <- manhattan_plot_ot +
  annotate("text", x=36340165, y=16, label="amh", color="red", size=5, hjust=0) +
  annotate("text", x=22382485, y=16, label="dmrta2", color="red", size=5, hjust=0) + 
  annotate("text", x=12051184, y=16, label="dmrt2b", color="red", size=5, hjust=0)


print(manhattan_plot_ot)





#rb 

rb_allo<-read.table("noINF_chr_processed_allo_rb_assoc.txt",header=T)


SNP <- c(1:(nrow(rb_allo)))
rb_allo <- data.frame(SNP,rb_allo)
rb_allo_23 <- subset(rb_allo, chr=="23")
rb_allo_4 <- subset(rb_allo, chr=="4")



head(rb_allo)
manhattan_plot_rb <- manhattan(rb_allo,chr="chr",bp="ps",p="X1e.17",
                               snp="SNP",logp=TRUE,ylab="-log10(p-value)",
                               suggestiveline = FALSE,
                               genomewideline = -log10(5e-08),
                               ylim = c(0, 20))

manhattan(rb_allo_23,chr="chr",bp="ps",p="X1e.17",
          snp="SNP",logp=TRUE,ylab="-log10(p-value)",
          suggestiveline = FALSE,
          ylim = c(0, 20),
          genomewideline = -log10(5e-08))

,
highlight=c(8985:9005))
highlight =c(113190:113208,179080:179148))
dmrta2=113196:113203,
amh=179088:179140
manhattan_plot_ot <- manhattan_plot_ot +
  annotate("text", x=36340165, y=16, label="amh", color="red", size=5, hjust=0) +
  annotate("text", x=22382485, y=16, label="dmrta2", color="red", size=5, hjust=0) + 
  annotate("text", x=12051184, y=16, label="dmrt2b", color="red", size=5, hjust=0)


