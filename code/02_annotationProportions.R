library(data.table)
library(dplyr)
library(GenomicRanges)
library(BuenColors)
library(diffloop)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

"%ni%" <- Negate("%in%")

# Import annotation
x01 <- bedToGRanges("../annotations/UCSC_3primeUTR.bed")
x02 <- bedToGRanges("../annotations/UCSC_Exons.bed")
x03 <- bedToGRanges("../annotations/UCSC_5primeUTR.bed")
x04 <- bedToGRanges("../annotations/UCSC_Introns.bed")

# Take a sample and run with it
getProportions <- function(file, sample){
  gr_t <- addchr(makeGRangesFromDataFrame(read.table(file, header = TRUE)))
  seq_t <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gr_t)

  # Do overlaps
  ov_1 <- findOverlaps(gr_t, x01)
  ov_2 <- findOverlaps(gr_t, x02)
  ov_3 <- findOverlaps(gr_t, x03)
  ov_4 <- findOverlaps(gr_t, x04)
  
  # Classify each variant
  class <- ifelse(1:length(gr_t) %in% queryHits(ov_1), "3UTR",
                  ifelse(1:length(gr_t) %in% queryHits(ov_2), "Exon",
                         ifelse(1:length(gr_t) %in% queryHits(ov_3), "5UTR",
                                ifelse(1:length(gr_t) %in% queryHits(ov_4), "Intron", "other"))))
  
  data.frame(class) %>%
    group_by(class) %>% 
    summarise (n = n()) %>%
    filter(class != "other") %>%
    mutate(freq = n / sum(n)) %>%
    reshape2::melt() %>% filter(variable == "freq") %>%
    mutate(sample = sample)
}

propdf <- rbind(getProportions("../output/LIN28B-IDR.bed", "Consensus"),
      getProportions("../output/LIN28B-Rep1.bed", "Rep1"),
      getProportions("../output/LIN28B-Rep2.bed", "Rep2"))

easter_colors <- c("3UTR" =	"#e0cdff", "Exon" = "#dcf9a8", "Intron" = "#c1f0fb")
ggplot(propdf, aes(x = sample, y = value*100, fill = class)) + 
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c( "gray30", "gray40", "lightgrey")) +
  labs(x = "", y = "% of peak annotations", fill = "")  + pretty_plot(fontsize = 8) +
  L_border() -> propPlot
cowplot::ggsave(propPlot, file = "../output/plots/proportion_plot.pdf", width = 2.2, height = 2)
