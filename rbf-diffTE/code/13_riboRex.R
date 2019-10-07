library(riborex)
library(dplyr)
library(EnhancedVolcano)
load("../output/24March2019-rbf_rna-counts-CD34.rda")

condition1 <- "Adult"
condition2 <- "CordBlood"
# Subset out the specific sample conditions
condition <- c(rep(condition1, 4), rep(condition2, 4))
common <- intersect(rownames(RNA.counts), rownames(RBF.counts))
rbf <- RBF.counts[common,c(colnames(RBF.counts)[grep(condition1, colnames(RBF.counts))], colnames(RBF.counts)[grep(condition2, colnames(RBF.counts))])]
rna <- RNA.counts[common,c(colnames(RNA.counts)[grep(condition1, colnames(RNA.counts))], colnames(RNA.counts)[grep(condition2, colnames(RNA.counts))])]

res.deseq2 <- riborex(rna, rbf, condition, condition)
results.corrected <- correctNullDistribution(res.deseq2)

# Manual
results.corrected %>% data.frame() %>%
  mutate(gene = results.corrected@rownames) %>%
  arrange(padj) -> processed_df


title <- paste0(condition1, "-vs-", condition2)
EnhancedVolcano(data.frame(results.corrected) ,
                lab =  results.corrected@rownames,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-5, 5),
                ylab = bquote(~-log[10]~adjusted~italic(P)),
                title=title,
                pCutoff = 0.1,
                col=c('black', 'black', 'black', 'red3'),
                FCcutoff = 1.25,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) + theme(legend.position= "none")

write.table(processed_df, file = "../output/24March2019-diffTE-cd34.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
