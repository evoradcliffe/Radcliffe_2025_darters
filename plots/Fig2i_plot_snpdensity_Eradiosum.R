#Load packages
library(ggplot2)
library(RcppRoll)

setwd("/Users/rlmoran/Desktop/darters/SDR/SNP_density/OB")

rm(list=ls())
ls() 

#snp data
OBAuto <- read.csv("OB_snpdensity_fc_autosomes_1.txt", header=T)
OBSex <- read.csv("OB_snpdensity_fc_sexchromosomes_1.txt", header=T)
#OBAuto <- read.csv("Caddo_snpdensity_fc_autosomes.txt", header=T)
#OBSex <- read.csv("Caddo_snpdensity_fc_sexchromosomes.txt", header=T)


#sort data
OBAutoSorted <-OBAuto[order(OBAuto$Chromosome,OBAuto$WindowStart),]
OBSexSorted <-OBSex[order(OBSex$Chromosome,OBSex$WindowStart),]

#significance tests
median(OBAutoSorted$MFLogaverage)
median(OBSexSorted$MFLogaverage)
mean(OBAutoSorted$MFLogaverage)
mean(OBSexSorted$MFLogaverage)
wilcox.test(OBAutoSorted$MFLogaverage, OBSexSorted$MFLogaverage)


pdf_path="SNP_density_chromosome_results_OB_10kb_AutoSexDensityPlot.pdf"

pdf(file=pdf_path,onefile=TRUE)

grid.arrange(AEs_snpdensity_plot, MI_RB_snpdensity_plot, OB_snpdensity_plot, nrow=1)


#plot
OB_snpdensity_plot <- ggplot(OBAutoSorted, aes(x=MFLogaverage)) + 
  geom_density(colour="#000000",fill="#5F808D",size=0.7,alpha=0.7) +
  geom_density(data=OBSexSorted, aes(x=MFLogaverage),colour="#000000", fill="#2CA25F", alpha=0.7, size=0.7)+
  geom_vline(xintercept = mean(OBAutoSorted$MFLogaverage), color="#5F808D", linetype="dashed", size=0.8) +
  geom_vline(xintercept = mean(OBSexSorted$MFLogaverage), color="#2CA25F", linetype="dashed", size=0.8) +
  coord_cartesian(ylim = c(0,40),xlim=c(-0.2,0.2)) + #Adjust accordingly for each dataset
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
  xlab(expression('M:F log'[2]*' normalized SNP density '*italic("E. radiosum")*', Chr 22 vs. Autosomes')) +
  ylab("Density") +
  geom_text(x=0.08, y=32.5, label="W = 66404706 \n p-value < 2.2e-16 \n ***", size=4) #Remove if non-significant difference
#600x430

#
OB_snpdensity_plot <- ggplot(OBAutoSorted, aes(x=MFLogaverage)) + 
  geom_density(colour="#000000", fill="#5F808D", size=0.7, alpha=0.7, bw=0.03) +
  geom_density(data=OBSexSorted, aes(x=MFLogaverage),
               colour="#000000", fill="#2CA25F", alpha=0.7, size=0.7, bw=0.03) +
  geom_vline(xintercept = mean(OBAutoSorted$MFLogaverage, na.rm=TRUE),
             color="#5F808D", linetype="dashed", size=0.8) +
  geom_vline(xintercept = mean(OBSexSorted$MFLogaverage, na.rm=TRUE),
             color="#2CA25F", linetype="dashed", size=0.8) +
  coord_cartesian(ylim = c(0,8), xlim=c(-0.5,0.5)) +  # match the other figures
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.y=element_text(margin=margin(0,5,0,0)),
        text=element_text(size=12),
        plot.margin = unit(c(2.5,6,1,3),"lines"),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.text=element_text(size=10)) +
  scale_y_continuous(expand = c(0,0)) + 
  xlab(expression('M:F log'[2]*' normalized SNP density '*italic("E. radiosum")*', Chr 22 vs. Autosomes')) +
  ylab("Density") +
  geom_text(x=0.15, y=2.5,
            label="W = 66404706, p-value < 2.2e-16 \n ***", size=4)

print(snpdensity_plot)
