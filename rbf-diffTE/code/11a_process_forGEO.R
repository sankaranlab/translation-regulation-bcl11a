library(dplyr)
library(data.table)

load("../output/24March2019-rbf_rna-counts-CD34.rda")

rbfnames <- c("Adult1-ribosomeprofiling","Adult2-ribosomeprofiling","Adult3-ribosomeprofiling","Adult4-ribosomeprofiling",
              "CordBlood1-ribosomeprofiling","CordBlood2-ribosomeprofiling","CordBlood3-ribosomeprofiling","CordBlood4-ribosomeprofiling")

rnanames <- c("Adult1-RNA","Adult2-RNA","Adult3-RNA","Adult4-RNA",
              "CordBlood1-RNA","CordBlood2-RNA","CordBlood3-RNA","CordBlood4-RNA")

df <- cbind(RBF.counts, RNA.counts)
df$gene <- rownames(RNA.counts)
colnames(df) <- c(rbfnames, rnanames, "gene")

write.table(df, file = "../output/RBF_RNAcounts_ADCB.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
