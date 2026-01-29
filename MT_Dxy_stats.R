

library(ape)
library(seqinr)


setwd("~/Documents/Manuscripts/NatComms_SDR")


aln <- read.dna("Oligo_MT_genome_alignment.fasta", format="fasta")  # alignment
D <- dist.dna(aln, model="raw", pairwise.deletion=T, as.matrix=TRUE)

D  # p-distance matrix (proportion differences)

L <- ncol(aln)
fixed_diffs <- round(D * L)
fixed_diffs

#

library(ape)

aln <- read.dna("Oligo_MT_genome_alignment.fasta", format="fasta")
L <- ncol(aln)

Ndiff <- dist.dna(aln, model="N", pairwise.deletion=FALSE, as.matrix=TRUE)
p     <- Ndiff / L

Ndiff
p

image(aln)                 # visual scan for messy regions
table(as.character(aln))   # see if any '-' or lots of 'n' remain

D_JC <- dist.dna(aln, model="JC69", pairwise.deletion=FALSE, as.matrix=TRUE)
D_K80 <- dist.dna(aln, model="K80", pairwise.deletion=FALSE, as.matrix=TRUE)

D_JC <- dist.dna(aln, model="JC69", pairwise.deletion=FALSE, as.matrix=TRUE)
D_K80 <- dist.dna(aln, model="K80",  pairwise.deletion=FALSE, as.matrix=TRUE)
D_JC
D_K80

aln_char <- as.character(aln)

# how much missing data per sequence?
miss <- t(apply(aln_char, 1, function(x) c(N=sum(x=="N"), gap=sum(x=="-"))))
miss

# confirm the 337 count directly from the alignment (no rounding)
sum(aln_char["OR552052.1EcaeruleumMitochondrion", ] !=
      aln_char["MW856838.1EcaeruleumMitochondrion", ])


aln_char <- as.character(aln)

x <- aln_char["OR552052.1EcaeruleumMitochondrion", ]
y <- aln_char["MW856838.1EcaeruleumMitochondrion", ]

# 1) What your manual count is doing (includes gap-vs-base as diffs)
sum(x != y)

# 2) Count only A/C/G/T vs A/C/G/T (this should match dist.dna model="N")
keep_bases <- x %in% c("A","C","G","T") & y %in% c("A","C","G","T")
sum(x[keep_bases] != y[keep_bases])

# 3) How many sites are excluded because of gaps/ambiguous (should be small)
sum(!keep_bases)

# 4) How many of the manual diffs are indel-coded (gap involved)
gap_involved <- (x == "-" | y == "-")
c(
  gap_sites_total = sum(gap_involved),
  gap_diffs       = sum(gap_involved & x != y),
  base_diffs      = sum(keep_bases & x != y)
)


sort(unique(x))
sort(unique(y))
xU <- toupper(x)
yU <- toupper(y)
bases <- c("A","C","G","T")

keep_bases <- xU %in% bases & yU %in% bases

c(
  L = length(xU),
  kept = sum(keep_bases),
  excluded = sum(!keep_bases),
  base_diffs = sum(xU[keep_bases] != yU[keep_bases]),
  gap_diffs  = sum((xU=="-" | yU=="-") & xU != yU)
)




# 1) Read clustal alignment
library(ape)

# reads CLUSTAL into DNAbin (works well after you removed spaces in names)
aln <- read.dna("Oligo_MT_genome_alignment.aln", format = "clustal")

D <- dist.dna(aln, model = "raw", pairwise.deletion = T, as.matrix = TRUE)
D

L <- ncol(aln)
diff_sites <- round(D * L)
diff_sites

# m is DNAbin (raw). Convert to character letters first:
mc <- as.character(aln)          # returns a character matrix (a/c/g/t/-/n)
table(as.vector(mc))
m[1, 1:60]             # first 60 sites

ape::base.freq(aln, freq = FALSE)   # counts
ape::base.freq(aln, freq = TRUE)    # frequencies

mean(as.character(aln) == "-")


D <- dist.dna(aln, model="raw", pairwise.deletion=FALSE, as.matrix=TRUE)

D  # p-distance matrix (proportion differences)

L <- ncol(aln)
fixed_diffs <- round(D * L)
fixed_diffs






library(ape)

read_clustal_dnabin2 <- function(path) {
  x <- readLines(path, warn = FALSE)
  
  # drop header, blank lines
  x <- x[!grepl("^\\s*$", x)]
  x <- x[!grepl("^CLUSTAL", x, ignore.case = TRUE)]
  
  # drop consensus lines like "  **:*.*"
  x <- x[!grepl("^\\s*[\\*\\.:]+\\s*$", x)]
  
  # keep lines that contain at least some plausible sequence letters/gaps
  dna_pat <- "[ACGTRYKMSWBDHVNacgtrykmswbdhvn\\-\\?\\.]{10,}"
  x <- x[grepl(dna_pat, x)]
  
  # name = first non-space token (accession), seq = first long DNA-like token
  nm <- sub("^\\s*(\\S+).*", "\\1", x)
  seqtok <- regmatches(x, regexpr(dna_pat, x))
  seqtok <- gsub("\\.", "-", seqtok)  # treat '.' as gap if present
  seqtok <- toupper(seqtok)
  
  u <- unique(nm)
  seqs <- vapply(u, function(z) paste0(seqtok[nm == z], collapse = ""), "")
  mat <- do.call(rbind, strsplit(seqs, "")); rownames(mat) <- u
  
  as.DNAbin(mat)
}

aln <- read_clustal_dnabin2("Oligo_MT_genome_alignment.aln")
D <- dist.dna(aln, model="raw", pairwise.deletion=FALSE, as.matrix=TRUE)
D

L <- ncol(aln)
fixed_diffs <- round(D * L)
fixed_diffs

