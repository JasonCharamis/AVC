## Principal components analysis based on SNPs identified through DNA-seq ##

suppressPackageStartupMessages(library("edgeR"))
library(ggfortify)
library(ggrepel)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(vcfR)
library(vegan)

## Get input file from standard input ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0, file.exists(args))
f_counts <- args

## open file from stdin ##
SNPs_cleaned <- read.csv(file=args[1], header=TRUE)

## make first column as rownames ##
rownames(SNPs_cleaned) <- SNPs_cleaned[,1]
SNPs_cleaned <- SNPs_cleaned[,-1]

## transform colnames to same syntax ##
colnames(SNPs_cleaned) <- str_replace(colnames(SNPs_cleaned), "results.|.s.bam|_", "")

SNPs_scaled <- scale(SNPs_cleaned)
pca_scaled <- prcomp(SNPs_scaled)

## generate sample names for grouping replicates ##
xt <- as.data.frame(vegan::scores(pca_scaled))

groups <- gsub("_\\d+$|d\\d+", "", rownames(xt))

## add column with sample name to group replicates ##
xtl <- xt %>% add_column(Sample = groups)

## draw auto-PCA plot with color mappings for groups ##
print ( autoplot(pca_scaled, data=xtl, colour='Sample', legend.size=5) +  labs(title = "Principal Component Analysis using log2TPM values" ) + theme_bw() + geom_text_repel(aes(label=rownames(xt),color=Sample), show.legend = FALSE  ) +
        theme(plot.title=element_text(face="bold",hjust=0.5), legend.title = element_text(size=12), legend.text = element_text(size=12), axis.title=element_text(size=12))   )

dev.off()
