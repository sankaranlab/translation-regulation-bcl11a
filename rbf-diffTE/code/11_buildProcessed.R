library(dplyr)
library(data.table)

# Now deal with expression
rnafiles <- list.files("../rnaseq", pattern=".tab", full.names = TRUE)
samplenames <- gsub("_1ReadsPerGene.out.tab", "", list.files("../rnaseq/", pattern=".tab"))
RNA.counts <- sapply(rnafiles, function(file) {
  data.frame(fread(paste0("cat < ", file), skip = 4))[,2]
}) %>% data.matrix()
colnames(RNA.counts) <- samplenames
RNA.counts <- data.frame(RNA.counts)

rownames(RNA.counts) <- data.frame(fread(paste0("cat < ", "../rnaseq/Adult_RNAseq_1_1ReadsPerGene.out.tab"), skip = 4))[,1]


# Now deal with rbf counts -- twice
rnafiles <- list.files("../rbf/runA_counts", pattern=".tab", full.names = TRUE)
samplenames <- gsub("_1ReadsPerGene.out.tab", "", list.files("../rbf/runA_counts", pattern=".tab"))
RBF.counts <- sapply(rnafiles, function(file) {
  data.frame(fread(paste0("cat < ", file), skip = 4))[,2]
}) %>% data.matrix()
colnames(RBF.counts) <- samplenames
RBF.countsA <- RBF.counts

rnafiles <- list.files("../rbf/runA_counts", pattern=".tab", full.names = TRUE)
samplenames <- gsub("_1ReadsPerGene.out.tab", "", list.files("../rbf/runA_counts", pattern=".tab"))
RBF.counts <- sapply(rnafiles, function(file) {
  data.frame(fread(paste0("cat < ", file), skip = 4))[,2]
}) %>% data.matrix()
colnames(RBF.counts) <- samplenames
RBF.countsB <- RBF.counts

RBF.counts <- data.frame(RBF.countsA + RBF.countsB)

rownames(RBF.counts)  <- data.frame(fread(paste0("cat < ", "../rnaseq/Adult_RNAseq_1_1ReadsPerGene.out.tab"), skip = 4))[,1]

save(RNA.counts, RBF.counts, file = "../output/24March2019-rbf_rna-counts-CD34.rda")

