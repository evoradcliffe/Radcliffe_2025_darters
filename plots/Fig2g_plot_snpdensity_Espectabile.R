#Load packages
library(ggplot2)
library(RcppRoll)

setwd("/Users/rlmoran/Desktop/darters/SDR/SNP_density/OT2OT")


rm(list=ls())
ls() 

#snp data
Auto <- read.csv("AES_snpdensity_fc_autosomes_ChrOnly.txt", header=T)
Sex <- read.csv("AES_snpdensity_fc_sexchromosomes.txt", header=T)

#sort data
AutoSorted <-Auto[order(Auto$Chromosome,Auto$WindowStart),]
SexSorted <-Sex[order(Sex$Chromosome,Sex$WindowStart),]

#significance tests
median(AutoSorted$MFLogaverage)
median(SexSorted$MFLogaverage)
wilcox.test(AutoSorted$MFLogaverage, SexSorted$MFLogaverage)


pdf_path="SNP_density_chromosome_results_AES_10kb_AutoSexDensityPlot.pdf"

pdf(file=pdf_path,onefile=TRUE)


#plot
AEs_snpdensity_plot <- ggplot(AutoSorted, aes(x=MFLogaverage)) + 
     geom_density(colour="#000000",fill="#5F808D",size=0.7,alpha=0.7,bw=0.05) +
     geom_density(data=SexSorted, aes(x=MFLogaverage),colour="#000000", fill="darkorange", alpha=0.7, size=0.7,bw=0.05)+
     geom_vline(aes(xintercept=mean(AutoSorted$MFLogaverage)),color="#5F808D",linetype="dashed",size=0.8) +
     geom_vline(aes(xintercept=mean(SexSorted$MFLogaverage)),color="darkorange",linetype="dashed",size=0.8) +
     coord_cartesian(ylim = c(0.0,4.0),xlim=c(-1,1)) + #Adjust accordingly for each dataset
     theme(panel.grid.minor = element_blank(), panel.background = element_blank())+
     theme(
     		axis.title.y=element_text(margin=margin(0,5,0,0)),
           text=element_text(size=12),
           plot.margin = unit(c(2.5,6,1,3),"lines"),
          axis.line.y = element_line(color="black", size = 0.5),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.text=element_text(size=10)
           ) +
       scale_y_continuous(expand = c(0,0)) + 
  xlab(expression('M:F log'[2]*' normalized SNP density '*italic("E. spectabile")*', Chr 9 vs. Autosomes')) +
  ylab("Density") +
       geom_text(x=0.5, y=2.5, label="W = 105624206 \n p-value < 2.2e-16 \n ***", size=4) #Remove if non-significant difference

print(snpdensity_plot)
