#### Manhattan plot of sex association test results - any mention of "allo" or "rb" below can be replaced with any other corresponding .assoc.txt file and naming system.####
library(qqman)
library(ggplot2)

# Load GEMMA output
ot_all <- read.table("chr_processed_allo_rb_assoc.txt", header = TRUE)

# Ensure chromosome is numeric
# Remove "chr" prefix if present, and convert to numeric
ot_all$chr <- as.numeric(gsub("chr", "", ot_all$chr))

# Filter out chr == 0
ot_all <- subset(ot_all, chr != 0)

# Create SNP index
SNP <- seq_len(nrow(ot_all))
ot_dat <- data.frame(SNP, ot_all)

# Clean data
ot_dat_clean <- ot_dat[is.finite(ot_dat$p_wald) & !is.na(ot_dat$p_wald), ]

# Manhattan plot with alternating colors
manhattan(ot_dat_clean,
          chr = "chr",
          bp = "ps",
          p = "p_wald",
          snp = "SNP",
          logp = TRUE,
          ylab = "-log10(p-value)",
          suggestiveline = FALSE,
          genomewideline = -log10(5e-08),
          ylim = c(0, 15),
          col = c("black", "grey"))  # Alternating colors

#### This is script for plotting location of genes of interest. ####
##CHROMOSOME 9 GENES##
text(x=36344386.5, y=13, labels="amh", col="seagreen", cex=1, pos=4,font = 2)
text(x=22384089, y=13, labels="dmrta2", col="seagreen", cex=1, pos=4,font = 2)
text(x=12053296.5,y=15,labels="dmrt2b",col="seagreen",cex=1,pos=4,font=2)
text(x=10246512,y=16,labels="sox2",col="seagreen",cex=1,pos=4,font=2)
text(x=2653441,y=16.5,labels="sox14",col="seagreen",cex=1,pos=4,font=2)
text(x=19529739,y=15.75,labels="LOC116695774",col="seagreen",cex=1,pos=4,font=2)
text(x=19556440.5,y=15,labels="vtg3",col="seagreen",cex=1,pos=4,font=2)
text(x=1534427.5,y=14.25,labels="brdt",col="seagreen",cex=1,pos=4,font=2)
text(x=1606810.5,y=15,labels="tgfbr3",col="seagreen",cex=1,pos=4,font=2)
text(x=30182131.5,y=16,labels="lonp1",col="darkorange",cex=1,pos=4,font = 2)
text(x=37195911.5,y=16.75,labels="polrmt",col="darkorange",cex=1,pos=4,font = 2)
text(x=20425783.5,y=16.75,labels="mfn1b",col="darkorange",cex=1,pos=4,font = 2)
text(x=18804830.5,y=13.5,labels="opa1",col="darkorange",cex=1,pos=4,font = 2)
text(x=20476883,y=17.5,labels="mrpl47",col="darkorange",cex=1,pos=4,font = 2)
text(x=37082900,y=17.5,labels="micos13",col="darkorange",cex=1,pos=4,font = 2)
text(x=683500,y=16,labels="cpt2",col="darkorange",cex=1,pos=4,font = 2)
text(x=20463981.5,y=18.25,labels="aldh7a1",col="darkorange",cex=1,pos=4,font = 2)
text(x=20480532,y=19,labels="ndufb5",col="darkorange",cex=1,pos=4,font = 2)
text(x=30418040.5,y=14.25,labels="ndufa11",col="darkorange",cex=1,pos=4,font=2)
text(x=38169326.5,y=18.25,labels="gpx4a",col="darkorange",cex=1,pos=4,font = 2)
text(x=15367103,y=15,labels="ndufs7",col="darkorange",cex=1,pos=4,font=2)
text(x=15725207,y=16, labels="ndufa7",col="darkorange",cex=1,pos=4,font=2)

##CHROMOSOME 23 GENES##
text(x=7575182.5,y=16,labels="asun",col="seagreen", cex=1, pos=4,font = 2)
text(x=3413640,y=14,labels="sox5",col="seagreen", cex=1, pos=4,font = 2)
text(x=3099637,y=16.75,labels="atp23",col="darkorange",cex=1,pos=4,font = 2)
text(x=17273146.5,y=15,labels="ndufa5",col="darkorange",cex=1,pos=4,font = 2)
text(x=13235413.5,y=16,labels="ndufa12",col="darkorange",cex=1,pos=4,font = 2)
text(x=15734044.5,y=16.75,labels="ssbp1",col="darkorange",cex=1,pos=4,font = 2)
text(x=13339154,y=17.5,labels="tmem243b",col="darkorange",cex=1,pos=4,font = 2)
text(x=7609552.5,y=15,labels="mrpl42",col="darkorange",cex=1,pos=4,font = 2)
text(x=6099088,y=15,labels="pnpla8",col="darkorange",cex=1,pos=4,font = 2)
text(x=1848888,y=15,labels="LOC116673457",col="seagreen", cex=1, pos=4,font = 2)
text(x=13810692,y=15,labels="gata3",col="seagreen",cex=1,pos=4,font=2)
text(x=9314724,y=16,labels="LOC116673454",col="magenta",cex=1,pos=4,font=2)
text(x=15500376,y=16,labels="wnt-7b",col="seagreen",cex=1,pos=4,font = 2)
text(x=9104324,y=14.75,labels="exoc4",col="seagreen",cex=1,pos=4,font = 2)
text(x=9018994,y=13.75,labels="chchd3",col="darkorange",cex=1,pos=4,font=2)
text(x=1372805.5,y=16,labels="apex1",col="darkorange",cex=1,pos=4,font = 2)
text(x=11815775,y=15,labels="ndufb2",col="darkorange",cex=1,pos=4,font=2)

print(manhattan_plot_ot)
