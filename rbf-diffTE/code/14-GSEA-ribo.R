library(liger)
library(data.table)
library(dplyr)
library(BuenColors)

gs <- read.table("../ribo-geneset.txt", stringsAsFactors = FALSE)[,1]

x <- data.frame(fread("../output/24March2019-diffTE-cd34.tsv")) %>% arrange(log2FoldChange)
vals <- x$log2FoldChange ; names(vals) <- as.character(x$gene)
gsea(values=vals, LIN28B_targets, mc.cores=2, plot = TRUE) 
