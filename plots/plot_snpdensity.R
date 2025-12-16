#Load packages
library(ggplot2)
library(RcppRoll)

#Define a moving average function
movingaverage <- function (x, window) {
  ma <- roll_mean(x, window, fill=NA)
}

windowsizeCI = 50
windowsizeMovAv = 90

pdf_path="SNP_density_chromosome_results_SES_10kb.pdf"

pdf(file=pdf_path,onefile=TRUE)

##### Orangethroat Allopatric #####
#Load data
datacov <- read.csv("snppdensity_fc_10kb_SES.txt", row.names=1,header=T)
datacovsorted <- datacov[order(datacov$Chromosome,datacov$WindowStart),]

#List of chromosomes
chromosomes <- c('NC_045733.1',
                 'NC_045734.1',
                 'NC_045735.1',
                 'NC_045736.1',
                 'NC_045737.1',
                 'NC_045738.1',
                 'NC_045739.1',
                 'NC_045740.1',
                 'NC_045741.1',
                 'NC_045742.1',
                 'NC_045743.1',
                 'NC_045744.1',
                 'NC_045745.1',
                 'NC_045746.1',
                 'NC_045747.1',
                 'NC_045748.1',
                 'NC_045749.1',
                 'NC_045750.1',
                 'NC_045751.1',
                 'NC_045752.1',
                 'NC_045753.1',
                 'NC_045754.1',
                 'NC_045755.1',
                 'NC_045756.1')

for (i in chromosomes)
{
    #Calculate confidence interval for moving average
    datacovminusSexChr <- datacovsorted[datacovsorted$Chromosome!=i,]
    MFpermutescov <- replicate(1000,mean(sample(datacovminusSexChr$MFLogaverage,windowsizeCI,replace = FALSE)))
    MFCI25cov <- quantile(MFpermutescov,c(.025,.5,.975))[[1]]
    MFCI975cov <- quantile(MFpermutescov,c(.025,.5,.975))[[3]]
    
    #plot sex chromosome coverage data
    datacovSexChr <- datacovsorted[datacovsorted$Chromosome==i,]
    smoothlinefccov = movingaverage(datacovSexChr$MFLogaverage,windowsizeMovAv)
    paneldatacov <- as.data.frame(smoothlinefccov)
    paneldatacov$startmb <- datacovSexChr$WindowStart/1000000
    paneldatacov <- na.omit(paneldatacov)

    # Dynamically calculate y-axis range
    y_range <- range(datacovSexChr$MFLogaverage, na.rm = TRUE)

    plotfcdens <- ggplot(datacovSexChr, aes(x = WindowStart / 1000000, y = MFLogaverage)) +
    geom_rect(
    xmax = max(datacovSexChr$WindowStart / 1000000, na.rm = TRUE),
    xmin = min(datacovSexChr$WindowStart / 1000000, na.rm = TRUE),
    ymax = quantile(MFpermutescov, c(.025, .5, .975))[[3]],
    ymin = quantile(MFpermutescov, c(.025, .5, .975))[[1]],
    fill = "grey", alpha = 0.08) +
    geom_point(colour = "#919191", fill = "#919191", alpha = 0.5, cex = 0.4) +
    geom_line(data = paneldatacov, aes(x = startmb, y = smoothlinefccov), size = 0.6, inherit.aes = FALSE) +
    theme(panel.grid.minor = element_blank(), panel.background = element_blank()) +
    theme(
    text = element_text(size = 12),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    axis.line.y = element_line(color = "black", size = 0.3),
    axis.line.x = element_line(color = "black", size = 0.3),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)) +
    scale_y_continuous(limits=c(-5.0,5.0)) +
    scale_x_continuous(breaks = seq(0, 40, 5)) +
    xlab("Start position (Mb)") +
    ylab(expression("M:F log"[2] * " SNP density")) +
    ggtitle("Orangethroat Sympatric 10kb windows", subtitle = i)
print(plotfcdens)
