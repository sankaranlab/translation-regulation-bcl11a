library(data.table)
library(dplyr)
library(GenomicRanges)
library(BuenColors)
library(diffloop)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SummarizedExperiment)

"%ni%" <- Negate("%in%")

# Annotations
x01 <- bedToGRanges("../annotations/UCSC_3primeUTR.bed")
x02 <- bedToGRanges("../annotations/UCSC_Exons.bed")
x03 <- bedToGRanges("../annotations/UCSC_5primeUTR.bed")
x04 <- bedToGRanges("../annotations/UCSC_Introns.bed")

peaks <- addchr(makeGRangesFromDataFrame(read.table("../output/LIN28B-IDR.bed",
                                                    header = TRUE), keep.extra.columns = TRUE))
#peaks<- resize(peaks, width = 60, fix="center")

SE <- SummarizedExperiment::SummarizedExperiment(
  rowRanges = peaks, 
  colData = data.frame(samples = ""), 
  assays = list())

# Compute deviation scores for 5-mers
kmer_ix <- matchKmers(4, SE, genome = BSgenome.Hsapiens.UCSC.hg19)
ov_1 <- findOverlaps(peaks, x01)
ov_2 <- findOverlaps(peaks, x02)
ov_3 <- findOverlaps(peaks, x03)
ov_4 <- findOverlaps(peaks, x04)

class <- ifelse(1:length(peaks) %in% queryHits(ov_1), "3UTR",
                ifelse(1:length(peaks) %in% queryHits(ov_2), "Exon",
                       ifelse(1:length(peaks) %in% queryHits(ov_3), "5UTR",
                              ifelse(1:length(peaks) %in% queryHits(ov_4), "Intron", "other"))))

dff <- data.frame(class, data.matrix(assays(kmer_ix)[["matches"]])*1)

dff %>% group_by(class) %>%  summarise_all(funs(mean)) %>%
  dplyr::select(-one_of("class")) %>% data.frame() %>% t() -> propDF

propDF["CCTC",]
propDF["CTCC",]
plot_df <- data.frame(x = propDF[,1]*100, y = propDF[,2]*100, motif = rownames(propDF))
plot_df$color <- plot_df$x > 60 & plot_df$y > 60

p1 <- ggplot(plot_df, aes(x = x, y = y, color = color)) + geom_point(size = 0.2) +
  pretty_plot(fontsize = 8) +
  L_border() + labs(x = "% of Peaks - 3'UTR", y = "% of Peaks - Coding Exon") +
  scale_color_manual(values = c("lightgrey", "black")) + theme(legend.position = "none") +
  geom_hline(yintercept = 60, linetype = 2) + geom_vline(xintercept = 60,  linetype = 2)

cowplot::ggsave(p1, file = "../output/plots/4mer.pdf", width = 2, height = 2)
