library(ggplot2)
library(RcppRoll)

rm(list=ls())
ls() 

print("SES Chr 23")
#snp data
Auto <- read.csv("SES_snpdensity_fc_autosomes_no23.txt", header=T)
Sex <- read.csv("SES_snpdensity_fc_sexchrom_23.txt", header=T)

#sort data
AutoSorted <-Auto[order(Auto$Chromosome,Auto$WindowStart),]
SexSorted <-Sex[order(Sex$Chromosome,Sex$WindowStart),]

#significance tests
median(AutoSorted$MFLogaverage)
median(SexSorted$MFLogaverage)
wilcox.test(AutoSorted$MFLogaverage, SexSorted$MFLogaverage)

print("AES Chr 23")
#snp data
Auto <- read.csv("AES_snpdensity_fc_autosomes_no23.txt", header=T)
Sex <- read.csv("AES_snpdensity_fc_sexchrom_23.txt", header=T)

#sort data
AutoSorted <-Auto[order(Auto$Chromosome,Auto$WindowStart),]
SexSorted <-Sex[order(Sex$Chromosome,Sex$WindowStart),]

#significance tests
median(AutoSorted$MFLogaverage)
median(SexSorted$MFLogaverage)
wilcox.test(AutoSorted$MFLogaverage, SexSorted$MFLogaverage)


print("AEP Chr 23")
#snp data
Auto <- read.csv("AEP_snpdensity_fc_autosomes_no23.txt", header=T)
Sex <- read.csv("AEP_snpdensity_fc_sexchrom_23.txt", header=T)

#sort data
AutoSorted <-Auto[order(Auto$Chromosome,Auto$WindowStart),]
SexSorted <-Sex[order(Sex$Chromosome,Sex$WindowStart),]

#significance tests
median(AutoSorted$MFLogaverage)
median(SexSorted$MFLogaverage)
wilcox.test(AutoSorted$MFLogaverage, SexSorted$MFLogaverage)
