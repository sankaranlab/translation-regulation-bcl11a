library(BuenColors)
library(ggrastr)
library(GenomicRanges)
library(rtracklayer)

df <- read.table("../output/24March2019-diffTE-cd34.tsv", header = TRUE)
df$BCL11A <- df$gene == "BCL11A"

gtf <- import("../other/hg19_10X.gtf.gz")
lin28b_bed <- diffloop::rmchr(import("../other/LIN28B-IGV_ready.bed"))
ov <- findOverlaps(gtf, lin28b_bed)
LIN28B_targets <- unique(mcols(gtf)$gene_name[queryHits(ov)])

df$LIN28B_target <-  df$gene %in% LIN28B_targets

p1 <- ggplot(df %>% arrange(BCL11A), aes(x = log2FoldChange, y = -log10(padj), color = BCL11A)) +
  geom_point_rast(size = 5, dpi = 500) +
  scale_color_manual(values = c("black", "firebrick")) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") 

p2 <- ggplot(shuf(df), aes(x = log2FoldChange, y = -log10(padj), color = LIN28B_target)) +
  geom_point_rast(size = 5, dpi = 500) +
  scale_color_manual(values = c("black", "firebrick")) +
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none") 

cowplot::ggsave(cowplot::plot_grid(p1, p2, nrow = 1), 
                width = 4, height = 2, filename = "../output/riboProfile_Lin28B.pdf")
