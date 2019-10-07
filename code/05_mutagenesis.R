library(data.table)
library(dplyr)
library(BuenColors)

rep1 <- read.table("../mutations/LIN28B_CB2_rRNAAligned.out.bam.all_positions.txt", header = TRUE)
rep2 <- read.table("../mutations/LIN28B_CB2Aligned.out.bam.all_positions.txt", header = TRUE)

total_df <- data.frame(
  coordinate = rep1$BP,
  TOTAL = rep1$TOTAL + rep2$TOTAL,
  REP1_altAll = apply((rep1[,c("A", "C", "G", "T")])/rep1$TOTAL, 1, function(X) sort(X, decreasing = TRUE)[2]),
  REP2_altAll = apply((rep2[,c("A", "C", "G", "T")])/rep2$TOTAL, 1, function(X) sort(X, decreasing = TRUE)[2])
)

mdf <- reshape2::melt(total_df[,c("coordinate", "REP1_altAll", "REP2_altAll")], id.vars = "coordinate")

p1 <- ggplot(total_df, aes(x = coordinate, y = TOTAL)) + geom_line() +
  pretty_plot() + L_border() + theme(legend.position = "none") + labs(x = "", y = "Read count")

p2 <- ggplot(mdf, aes(x = coordinate, y = -1*value*100, color = variable)) + geom_line() +
  pretty_plot() + L_border() + theme(legend.position = "none") + labs(x = "", y = "% alternate allele") +
  scale_color_manual(values = c("firebrick", "dodgerblue3"))

cowplot::ggsave(cowplot::plot_grid(p1, p2, nrow = 2), file = "../output/mutationRate.pdf", width = 4, height = 4)
