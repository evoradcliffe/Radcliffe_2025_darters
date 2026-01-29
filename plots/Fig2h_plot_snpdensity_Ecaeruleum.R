#Load packages
library(ggplot2)
library(RcppRoll)

rm(list=ls())
ls() 
setwd("/Users/rlmoran/Desktop/darters/SDR/SNP_density/RB2RB")

#snp data
Auto <- read.csv("snpdensity_All_RBs_fc_autosomes_ChrOnly.txt", header=T)
Sex <- read.csv("snpdensity_All_RBs_fc_sexchromosomes.txt", header=T)

# Convert  Inf to NA (only numeric columns)
Auto[] <- lapply(Auto, function(x) if (is.numeric(x)) replace(x, is.infinite(x), NaN) else x)
Sex[] <- lapply(Sex, function(x) if (is.numeric(x)) replace(x, is.infinite(x), NaN) else x)

#sort data
AutoSorted <-Auto[order(Auto$Chromosome,Auto$WindowStart),]
SexSorted <-Sex[order(Sex$Chromosome,Sex$WindowStart),]

#significance tests
median(AutoSorted$MFLogaverage, na.rm = TRUE)
median(SexSorted$MFLogaverage, na.rm = TRUE)
wilcox.test(AutoSorted$MFLogaverage, SexSorted$MFLogaverage)


pdf_path="SNP_density_chromosome_results_AllRB2RB_10kb_AutoSexDensityPlot.pdf"

pdf(file=pdf_path,onefile=TRUE)


#plot
snpdensity_plot <- ggplot(AutoSorted, aes(x=MFLogaverage)) + 
     geom_density(colour="#000000",fill="#5F808D",size=0.7,alpha=0.7) +
     geom_density(data=SexSorted, aes(x=MFLogaverage),colour="#000000", fill="#E1B91A", alpha=0.7, size=0.7)+
     geom_vline(aes(xintercept=median(AutoSorted$MFLogaverage, na.rm = TRUE)),color="#5F808D",linetype="dashed",size=0.8) +
     geom_vline(aes(xintercept=median(SexSorted$MFLogaverage, na.rm = TRUE)),color="#E1B91A",linetype="dashed",size=0.8) +
     coord_cartesian(ylim = c(0,10),xlim=c(-1,1)) + #Adjust accordingly for each dataset
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
       xlab(expression('M:F log'[2]*' normalised SNP density All IL and MI RB aligned to RB')) +
       ylab("Density") +
       geom_text(x=0.5, y=1.5, label="***", size=5) #Remove if non-significant difference

print(snpdensity_plot)
