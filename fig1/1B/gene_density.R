library(RIdeogram)
library(Cairo)
library(ggplot2)
library(rsvg)

setwd("E:\\a")

karyotype <- read.delim("E:\\a\\Dofficinale_karyotype2.txt",header = TRUE,stringsAsFactors = F)

density <- read.delim("E:\\a\\gene_density_500k.txt",header = FALSE,stringsAsFactors = F)
telomere <- read.delim("E:\\a\\telomere_plot.txt",header = TRUE,stringsAsFactors = F)
colnames(karyotype) <- c("Chr", "Start", "End")
colnames(density) <- c("Chr", "Start", "End", "Value")

ideogram(karyotype = karyotype, 
         overlaid = density,
         label= telomere,
         label_type = 'marker',
         output = 'chromosome.svg')

convertSVG("chromosome.svg", 
           file = 'chromosome2',
           device = "png")


rsvg_pdf("chromosome.svg", "chromosome2.pdf")

ideogram(karyotype = karyotype, 
         overlaid = density,
         label= telomere,
         label_type = 'marker',
         output = 'chromosome_new.svg')

convertSVG("chromosome_new.svg", 
           file = 'chromosome_new',
           device = "png")


rsvg_pdf("chromosome_new.svg", "chromosome_new.pdf")


