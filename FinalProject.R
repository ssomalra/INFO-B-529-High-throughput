library(DESeq2)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

colnames <- c("chr", "m6A_start",	"m6A_end", "strand", "annotation", "gene_start", "gene_end", "gene_name", "transcript_id", "description")
m6A_df <- read.csv('~/Documents/IUPUI/24-25 Courses/High-Throughput/FinalProject/m6A_annotations_clean.tsv', sep='\t', col.names=colnames)

# GENERATE COUNT TABLE OF # m6A SITES PER GENE =================================
# get unique m6A sites for each gene
m6A_unique_df <- m6A_df %>%
  distinct(gene_name, m6A_start, m6A_end, .keep_all = TRUE)  # Keep the first occurrence of each unique m6A position per gene

m6A_gene_counts <- m6A_unique_df %>%
  group_by(gene_name) %>%
  summarise(
    description = unique(description)[1],  # Getting the first unique description
    m6A_count = n(),
  ) %>%
  ungroup()

# sort the genes by the number of m6A sites
m6A_gene_counts <- m6A_gene_counts %>%
  arrange(desc(m6A_count))

write.table(m6A_gene_counts, file = "~/Documents/IUPUI/24-25 Courses/High-Throughput/FinalProject/m6A_gene_counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# M6A LOCATION DENSITY PLOT ====================================================
calculate_m6a_location <- function(df) {
  df <- df %>%
    mutate(
      gene_length = gene_end - gene_start,
      m6A_distance = m6A_start - gene_start,
      m6a_location = m6A_distance / gene_length,
      m6a_location = if_else(strand == "-", 1 - m6a_location, m6a_location)
    )
  return(df)
}

m6A_unique_df <- calculate_m6a_location(m6A_unique_df)

# calculate density
m6A_d <- density(m6A_unique_df$m6a_location)


# B_WT_IP
png(filename = '~/Documents/IUPUI/24-25 Courses/High-Throughput/FinalProject/density_m6A_location.png', res=500, width = 4000, height = 3200)
plot(m6A_d, main = "m6A Distribution Across Normalized Gene Lengths", 
     xlim = c(0, 1), ylim = c(0, 2.0),
     xlab='m6A Distribution Across Normalized Gene Lengths', ylab = "Density", col = "plum4", lwd = 2)
polygon(m6A_d, col = rgb(0.7, 0.5, 0.9, 0.4), border = NA)
dev.off()

# TOP 10 GENES WITH HIGHEST NUMBER OF M6A SITES =============================
m6A_gene_counts <- m6A_gene_counts %>%
  mutate(label = paste0(gene_name, " (", description, ")"))

# select top 10 genes
top10_genes <- m6A_gene_counts %>%
  arrange(desc(m6A_count)) %>%
  slice_head(n = 10)

png(filename = '~/Documents/IUPUI/24-25 Courses/High-Throughput/FinalProject/top10_genes.png', res=500, width = 7000, height = 3200)
ggplot(top10_genes, aes(x = reorder(label, m6A_count), y = m6A_count)) +
  geom_bar(stat = "identity", fill = "#9165a8") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 10 Genes with Highest m6A Site Counts",
    x = "Gene (Description)",
    y = "Number of m6A Sites"
  ) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
dev.off()

# CALCULATING TRANSCRIPTIONAL EFFICIENCY =======================================
# read in count tables
rna_counts <- read.table("~/Documents/IUPUI/24-25 Courses/High-Throughput/FinalProject/translational_efficiency/rnaseq_featureCounts.txt", header=TRUE, row.names=1)
ribo_counts_new <- read.table("~/Documents/IUPUI/24-25 Courses/High-Throughput/FinalProject/translational_efficiency/featureCounts_output_new.txt", header=TRUE, row.names=1)

# extract just the count columns (columns 6 to 9)
ribo_counts_new <- ribo_counts_new[, 6:9]
rna_counts  <- rna_counts[, 6:9]

# combine the count tables
combined_counts_new <- cbind(rna_counts, ribo_counts_new)

write.table(combined_counts_new, "~/Downloads/count_table.tsv", sep='\t', row.names=TRUE,col.names=TRUE, quote=FALSE)

# create meta data
condition_deseq_new <- data.frame(
  row.names = colnames(combined_counts_new),
  condition = c("KO","KO","WT","WT","KO","KO","WT","WT"),
  seqtype = rep(c("RNA", "Ribo"), each = 4)
)

# use DESeq Size Factor to normalize
dds.all.new <- DESeqDataSetFromMatrix(countData = combined_counts_new,
                                  colData = condition_deseq_new,
                                  design = ~ condition+seqtype+condition:seqtype)
keep.new <- rowSums(counts(dds.all.new)) >= 1
dds.all.new <- dds.all.new[keep.new,]
dds.all.new <- DESeq(dds.all.new)

## extract the normalized read counts for the translation efficiency calculation
count.df.new <- counts(dds.all.new, normalize = T) %>% as.data.frame()

## mean value for each condition
count.df.new$mean_WT_rna <- (count.df.new$alignment.rnaseq_WT_1Aligned.sortedByCoord.out.bam + count.df.new$alignment.rnaseq_WT_2Aligned.sortedByCoord.out.bam)/2
count.df.new$mean_KO_rna <- (count.df.new$alignment.rnaseq_KO_1Aligned.sortedByCoord.out.bam + count.df.new$alignment.rnaseq_KO_2Aligned.sortedByCoord.out.bam)/2
count.df.new$mean_WT_ribo <- (count.df.new$alignment_new.riboseq_WT_1Aligned.sortedByCoord.out.bam + count.df.new$alignment_new.riboseq_WT_2Aligned.sortedByCoord.out.bam)/2
count.df.new$mean_KO_ribo <- (count.df.new$alignment_new.riboseq_KO_1Aligned.sortedByCoord.out.bam + count.df.new$alignment_new.riboseq_KO_2Aligned.sortedByCoord.out.bam)/2

## TE calculation
count.df.new$TE_wt <- count.df.new$mean_WT_ribo/count.df.new$mean_WT_rna
count.df.new$TE_ko <- count.df.new$mean_KO_ribo/count.df.new$mean_KO_rna

# filter the dataframe
filtered_df <- count.df.new[
  !( (count.df.new$TE_wt == 0 & count.df.new$TE_ko == 0) |  # both zero
       is.nan(count.df.new$TE_wt) | is.nan(count.df.new$TE_ko) |  # either NaN
       is.infinite(count.df.new$TE_wt) | is.infinite(count.df.new$TE_ko)  # either Inf
  ), 
]

filtered_df$TE_log2FC <- log2((filtered_df$TE_ko)/(filtered_df$TE_wt))

# is there a significant difference between WT and KO?
mean(filtered_df$TE_wt)
mean(filtered_df$TE_ko)
wilcox.test(filtered_df$TE_wt, filtered_df$TE_ko)

filtered_df$log_TE_wt <- log(filtered_df$TE_wt)
filtered_df$log_TE_ko <- log(filtered_df$TE_ko)

filtered_df <- as_tibble(filtered_df)

te_long_new <- filtered_df %>%
  dplyr::select(TE_wt, TE_ko) %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "Condition",
                      values_to = "TE") %>%
  dplyr::mutate(Condition = ifelse(Condition == "TE_wt", "WT", "KD"))

png(filename = '~/Documents/IUPUI/24-25 Courses/High-Throughput/FinalProject/TE_WTvKD.png', res=500, width = 3200, height = 2800)
ggplot(te_long_new, aes(x = TE, color = Condition)) +
  stat_ecdf(geom = "line", size = 1) +
  scale_x_log10() +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  ) +
  labs(
    title = "Cumulative Translational Efficiency",
    x = "log(TE)",
    y = "Cumulative Proportion"
  ) +
  scale_color_manual(values = c("WT" = "#15a0bf", "KD" = "#f59c20")) +
  coord_cartesian(xlim = c(1, 30))
dev.off()

# CORRELATION BETWEEN M6A SITES AND TE =========================================
# add row names as a new column
count.df.new$gene_name <- rownames(count.df.new)

# merge with m6A data
merged_df <- dplyr::left_join(count.df.new, m6A_gene_counts, by = "gene_name")

# check for NAs and filter
merged_df <- merged_df %>%
  dplyr::filter(!is.na(m6A_count), !is.na(TE_wt), is.finite(TE_wt))

cor.test(merged_df$m6A_count, merged_df$TE_wt, method = "pearson")

png(filename = '~/Documents/IUPUI/24-25 Courses/High-Throughput/FinalProject/correlation.png', res=500, width = 3200, height = 2800)
ggplot(merged_df, aes(x = m6A_count, y = log10(TE_wt))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
dev.off()