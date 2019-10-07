library(data.table)
library(dplyr)
library(BuenColors)
library(ggbeeswarm)

ad1 <- data.frame(fread("../rna-seqdata/kallisto/Ad1-S1/abundance.tsv"))
ad2 <- data.frame(fread("../rna-seqdata/kallisto/Ad2-S1/abundance.tsv"))
cb1 <- data.frame(fread("../rna-seqdata/kallisto/CB1-S1/abundance.tsv"))
cb2 <- data.frame(fread("../rna-seqdata/kallisto/CB2-S1/abundance.tsv"))

df <- data.frame(
  transcript = ad1$target_id,
  AD1 = ad1$tpm,
  AD2 = ad2$tpm,
  CB1 = cb1$tpm,
  CB2 = cb2$tpm
)

BCL11A <- c("ENST00000538214",
            "ENST00000537768", 
            "ENST00000492272",
            "ENST00000489516",
            "ENST00000489183",
            "ENST00000477659",
            "ENST00000409351",
            "ENST00000359629",
            "ENST00000358510",
            "ENST00000356842",
            "ENST00000335712")

df %>% filter(transcript %in% BCL11A) %>%
  reshape2::melt(id.vars = "transcript") %>%
  mutate(condition = ifelse(substring(variable, 1, 1) == "A", "Adult", "Newborn")) %>%
  group_by(transcript, condition) %>%
  summarize(mean = mean(log2(value + 1)), sd = sd(log2(value + 1)),
            p1 = min(log2(value + 1)), p2 = max(log2(value + 1))) -> summary_df
sdf <- reshape2::melt(summary_df[,c("transcript", "condition", "p1", "p2")], id.vars = c("transcript", "condition"))
colnames(sdf) <- c("transcript", "condition", "variable", "mean")

p1 <- ggplot(summary_df, aes(x = transcript, y = mean, fill = condition, shape = condition)) +
  geom_bar(position=position_dodge(), stat="identity", color = "black", width = 0.4) +
  scale_fill_manual(values = c("lightgrey", "gray60")) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd),
                width=.2, position=position_dodge(0.4)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "", y = paste0("log2 TPM")) +
  geom_point(data = sdf, position = position_dodge2(width = 0.4, padding = 0.1), size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.position = "bottom")
cowplot::ggsave(p1, file = "../output/TPM-Kallisto.pdf", width = 4, height = 3)


df %>% 
  mutate(meanNB = (CB1 + CB2)/2, 
         meanAD = (AD1 + AD2)/2) -> allDF

