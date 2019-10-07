library(sleuth)
library(BuenColors)
library(dplyr)

base_dir <- "../rna-seqdata/kallisto"
sample_id <- dir(file.path(base_dir))
kal_dirs <- paste0(base_dir, "/", sample_id)

s2c <- data.frame(
 sample = c("Ad1", "Ad2", "CB1", "CB2"),
 condition = c("Ad", "Ad", "CB", "CB"), 
 path = kal_dirs,
 stringsAsFactors = FALSE
)

# Do biomart things to collapse transcripts to gene level counts
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'grch37.ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g,
                  aggregation_column = 'ens_gene', gene_mode = TRUE)

so$bs_summary$obs_tpm["ENSG00000187772",] %>%
  reshape2::melt() %>%
  mutate(condition = c("Adult", "Adult", "Newborn", "Newborn")) %>%
  group_by(condition) %>% summarize(mean = mean(log2(value + 1)),
                                    sd = sd(log2(value + 1))) -> summary_df

summary_df$condition <- factor(as.character(summary_df$condition), c("Newborn", "Adult"))
p1 <- ggplot(summary_df, aes(x = condition, y = mean)) +
  geom_bar(position=position_dodge(), stat="identity",
           color = "black", width = 0.4, fill = "grey") +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd),
                width=.2, position=position_dodge(.4)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "", y = paste0("log2 TPM")) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.position = "bottom")
ggsave(p1, file = "../output/TPM-LIN28B.pdf", width = 2, height = 2)

so$bs_summary$obs_tpm["ENSG00000131914",] %>%
  reshape2::melt() %>%
  mutate(condition = c("Adult", "Adult", "Newborn", "Newborn")) %>%
  group_by(condition) %>%
  summarize(mean = mean(log2(value + 1)), sd = sd(log2(value + 1)),
            p1 = min(log2(value + 1)), p2 = max(log2(value+1))) -> summary_df

summary_df$condition <- factor(as.character(summary_df$condition), c("Newborn", "Adult"))
p1 <- ggplot(summary_df, aes(x = condition, y = mean)) +
  geom_bar(position=position_dodge(), stat="identity",
           color = "black", width = 0.4, fill = "grey") +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd),
                width=.2, position=position_dodge(.4)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "", y = paste0("log2 TPM")) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.position = "bottom")
ggsave(p1, file = "../output/TPM-LIN28A.pdf", width = 2, height = 2)

so$bs_summary$obs_tpm["ENSG00000119866",] %>%
  reshape2::melt() %>%
  mutate(condition = c("Adult", "Adult", "Newborn", "Newborn")) %>%
  group_by(condition) %>% summarize(mean = mean(log2(value + 1)),
                                    sd = sd(log2(value + 1)),
                                    p1 = min(log2(value + 1)), p2 = max(log2(value+1))) -> summary_df

sdf <- reshape2::melt(summary_df[,c("condition",  "p1", "p2")], id.vars = c("condition"))
colnames(sdf) <- c("condition", "variable", "mean")


p1 <- ggplot(summary_df, aes(x = condition, y = mean)) +
  geom_bar(position=position_dodge(), stat="identity",
           color = "black", width = 0.4, fill = "grey") +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd),
                width=.2, position=position_dodge(.4)) +
  geom_point(data = sdf, position = position_dodge2(width = 0.4, padding = 0.1), size = 0.8) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "", y = paste0("log2 TPM")) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.position = "bottom")
cowplot::ggsave(p1, file = "../output/TPM-BCL11A.pdf", width = 2, height = 2)


two_df <- data.frame(Ad = log2(rowSums(so$bs_summary$obs_tpm[,c(1,2)]) + 1),
                     New =  log2(rowSums(so$bs_summary$obs_tpm[,c(3,4)]) + 1))         

pC <- ggplot(two_df, aes(x = New, y = Ad)) +
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2 TPM Newborn", y = "log2 TPM Adult") +
  scale_y_continuous( expand = c(0, 0.2)) +
  scale_x_continuous( expand = c(0, 0.2)) +
  geom_abline(intercept = 0, slope = 1, color = "firebrick")
cowplot::ggsave(pC, file = "../output/plots/f2a_correlation.pdf", height = 2.7, width = 2.7)

write.table(round(so$bs_summary$obs_tpm,2), file = "../output/TPM-gene.tsv",
            sep = "\t", row.names = TRUE,
            col.names = TRUE, quote = FALSE)
