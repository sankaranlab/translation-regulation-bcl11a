library(data.table)
library(dplyr)
library(GenomicRanges)
library(BuenColors)
library(Biostrings)
library(diffloop)
library(BSgenome.Hsapiens.UCSC.hg19)

"%ni%" <- Negate("%in%")

LIN28B_binding <- makeGRangesFromDataFrame(data.frame(
  chr = "2", start = 60688540, end = 60688629
))

importPeaks_IDR <- function(file){
  df <- read.table(file)
  df$V2 <- pmin(df$start1, df$start2)
  df$V3 <- pmax(df$stop1, df$stop2)

  x <- unique(makeGRangesFromDataFrame(df[df$IDR < 0.01,], seqnames.field = "chr1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE))
  x <- x[seqnames(x) %in% c(as.character(1:22), "X")]
  return(x)
}

importPeaks <- function(file){
  df <- data.frame(fread(file))
  x <- unique(makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE))
  x <- x[seqnames(x) %in% c(as.character(1:22), "X")]
  return(x)
}

lin_rna <- importPeaks_IDR("../idrCode/LIN28B_IDR-overlapped-peaks.txt")
p1 <- importPeaks("../peaks/PURA-Abcam_noDup_noLambda_peaks.narrowPeak")
p2 <- importPeaks("../peaks/PURA-Bethyl_noDup_noLambda_peaks.narrowPeak")

lin28b_RNA_NoDupNoLambda_go <- lin_rna[1:length(lin_rna) %ni%  queryHits(findOverlaps(lin_rna, c(p1, p2)))]
seq_lin28b_RNA_NoDupNoLambda_go <- getSeq(BSgenome.Hsapiens.UCSC.hg19, addchr(lin28b_RNA_NoDupNoLambda_go))

motif_RNA <- as.numeric((vcountPattern("GGAGA", seq_lin28b_RNA_NoDupNoLambda_go) + vcountPattern("TCTCC", seq_lin28b_RNA_NoDupNoLambda_go)) > 0)
sumdf_RNA <- data.frame(lin28b_RNA_NoDupNoLambda_go, motif = motif_RNA, BCL11A = 1:length(lin28b_RNA_NoDupNoLambda_go) %in%
                          queryHits(findOverlaps(lin28b_RNA_NoDupNoLambda_go, LIN28B_binding)))

sumdf_RNA_noPURA <- data.frame(lin_rna, 
                               BCL11A = 1:length(lin_rna) %in%
                          queryHits(findOverlaps(lin_rna, LIN28B_binding)))

sumdf_RNA_noPURA %>% arrange(IDR) %>% 
  mutate(percentile_log10Q = 1:n()/n(), rank = 1:n(), total = n()) %>% filter(BCL11A) -> rankDF_RNA

sumdf_RNA %>% filter(IDR < 0.01) %>% dim()

sumdf_RNA %>% arrange(IDR) %>% 
  filter(IDR < 0.01) %>% 
  mutate(rank = 1:n()) %>%
  arrange(BCL11A) %>% 
  ggplot(aes(x = rank, y = -log10(IDR), color = BCL11A)) + geom_point(size = 0.2) +
  pretty_plot(fontsize = 8) + L_border() + labs(x = "Rank Sorted Peaks", y = "-log10 IDR") +
  scale_color_manual(values = c("grey", "black")) + theme(legend.position = "none") -> IDRP

cowplot::ggsave(cowplot::plot_grid(IDRP, nrow = 1), file = "../output/plots/rankSorted_IDR.pdf",  width = 2, height = 2)

out_idr <- sumdf_RNA[sumdf_RNA$IDR < 0.01 ,c(1,2,3,14,15)]; colnames(out_idr) <- c("chr", "start", "end", "IDR", "GGAGA_motif")
write.table(out_idr, file = "../output/LIN28B-IDR.bed", 
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

gro <- addchr(makeGRangesFromDataFrame(out_idr))
gro <- sortSeqlevels(gro); gro <- sort(gro)

write.table(data.frame(gro)[,c(1,2,3)], file = "../output/LIN28B-IGV_ready.bed", 
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
 
