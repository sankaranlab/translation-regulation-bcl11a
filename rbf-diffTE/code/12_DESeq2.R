library(dplyr)
library(DESeq2)
library(BiocParallel)
library(IHW)
library(annotables)
register(MulticoreParam(2))

# Import raw expression values
raw <- read.table("../output/RNAseq_rawGeneCounts.tsv", header = TRUE)
RNA.counts <- raw[,1:12]
meta <- stringr::str_split_fixed(colnames(RNA.counts), "_", 4) %>% data.frame()
meta$population <- c(rep("FMR1-0507", 3), rep("FMR1-3731", 3), rep("FMR1-DMSO", 3), rep("Parent", 3))

# RNA DEseq2 setup
RNA.counts.df <- as.data.frame(RNA.counts)

# Establish column data
RNA.condition <-meta$population
colData <- as.data.frame(RNA.condition)
row.names(colData) <- colnames(RNA.counts.df)

# Run DEseq2
RNA.dds <- DESeqDataSetFromMatrix(countData = RNA.counts.df, colData = colData, design = ~ RNA.condition)
RNA.dds <- DESeq(RNA.dds, parallel = TRUE)

# All pairwise comparisons
comp <- expand.grid(unique(meta$population), unique(meta$population), stringsAsFactors = FALSE)
comp <- subset(comp, Var1!=Var2)
comp$cond <- "RNA.condition"
comp <- comp %>% 
  mutate(key = paste0(pmin(Var1, Var2), "_" ,pmax(Var1, Var2), sep = "")) %>%
  dplyr::distinct(key, .keep_all=TRUE)

# Function to do all pair-wise DEseq comparisons
all.pairwise <- function(x) {
  print(x)
  # Do the DESeq contrast
  resSig.RNA <- data.frame(results(RNA.dds, contrast=c(comp[x,3], comp[x,2], comp[x,1]), parallel = TRUE, filterFun=ihw, alpha = 0.05))
  resSig.RNA$gene <- raw[,13]
  resSig.RNA %>% filter(complete.cases(.)) %>% arrange(padj) -> out
  
  # Do some rounding
  out$baseMean <- round( out$baseMean, 1)
  out$log2FoldChange <- round( out$log2FoldChange, 1)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
   
  # Output table
  write.table(out[,c("gene", "baseMean", "log2FoldChange", "pvalue", "padj")],
              file = paste0("../output/RNA_DESeq2/RNAseq-", comp[x,4], "_all.tsv"), 
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  comp[x,4]
  
}
nothing <- lapply(1:dim(comp)[1], all.pairwise)
