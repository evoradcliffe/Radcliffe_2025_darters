library(qqman)
fst <- read.table("intersex_rb_sym_cleaned.weir.fst", header=TRUE)
## similarly for allo rb, sym ot, and allo ot => file-naming system is same ##
fstsubset <- fst[complete.cases(fst),]
SNP <- c(1:(nrow(fstsubset)))
mydf <- data.frame(SNP,fstsubset)

tiff(file="sympatric_orangethroat.tiff",
     width=6, height=4, units="in", res=100, type="cairo")
plot <- manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weir and Cockerham Fst",genomewideline=0.5,ylim=c(0,1.5))
print(plot)
dev.off()

### Zooming in to only one chromosome of interest - this is for chromosome 23 ###
  ##### can be repeated with switching out of CHROM == value ####
# Subset to only chromosome 23
chr23_df <- subset(mydf, CHROM == 23)
plot <- manhattan(chr23_df,
                  chr = "CHROM",
                  bp = "POS",
                  p = "WEIR_AND_COCKERHAM_FST",
                  snp = "SNP",
                  logp = FALSE,
                  ylab = "Weir and Cockerham Fst",
                  genomewideline = 0.5,
                  ylim = c(0, 1.5))

# Sort by FST descending
top_snps <- chr23_df[order(-chr23_df$WEIR_AND_COCKERHAM_FST), ]

# View top 10 SNPs (or any number you want)
head(top_snps, 15)

gene <- read.table("symIL_by_gene_p1.out", header=TRUE)
top <- gene[order(-gene$Mean_FST),]
top <- gene[!grepl("^LOC", gene$Gene), ]
top <- top[order(-top$Mean_FST), ]

print <- head(top, 25)

### FST BY CHROM ###
## entire genome, pooled allopatric and sympatric samples ##
library(dplyr)
library(ggplot2)
library(readr)

# Load .fst files
ot_fst <- read_tsv("fixed_intersex_ot_POOLED.weir.fst")  # <- your actual file
rb_fst <- read_tsv("fixed_intersex_rb_POOLED.weir.fst")

# Convert CHROM to numeric if it's like "NC_045733.1" (optional step)
ot_fst$CHROM <- as.numeric(as.character(ot_fst$CHROM))
rb_fst$CHROM <- as.numeric(as.character(rb_fst$CHROM))

# Filter out rows with NA
ot_fst_clean <- ot_fst %>% filter(!is.na(CHROM), !is.na(FST))
rb_fst_clean <- rb_fst %>% filter(!is.na(CHROM), !is.na(FST))
