# deseq2_strict_volcano.R

# 1) libraries
if (!requireNamespace("DESeq2",   quietly=TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ggplot2",  quietly=TRUE)) install.packages("ggplot2")
library(DESeq2)
library(ggplot2)

# 2) counts
fc        <- read.table("efaecium_featureCounts.txt",
                        header=TRUE, sep="\t", comment.char="#",
                        row.names=1, check.names=FALSE)
countData <- fc[, -(1:6)]
cleanCols <- sub("\\.sorted\\.bam$", "", basename(colnames(countData)))
colnames(countData) <- cleanCols

# 3) metadata
colData <- read.csv("sample_metadata.csv", row.names=1, stringsAsFactors=FALSE)
colData <- colData[colnames(countData), , drop=FALSE]
stopifnot(all(colnames(countData) == rownames(colData)))

# 4) DESeq2
dds  <- DESeqDataSetFromMatrix(countData, colData, design=~condition)
dds  <- dds[rowSums(counts(dds)) >= 10, ]
dds  <- DESeq(dds)
res  <- results(dds)

# 5) summary table (strict criteria)
res_strict   <- subset(res, padj < 0.001 & abs(log2FoldChange) > 1)
n_tested     <- sum(!is.na(res$padj))
n_strict     <- nrow(res_strict)
n_up_strict  <- sum(res_strict$log2FoldChange >  1)
n_down_strict<- sum(res_strict$log2FoldChange < -1)

summary_table <- data.frame(
  Category = c(
    "Tested after filtering",
    "Strict (padj<0.001 & |log₂FC|>1)",
    "Up in Serum (log₂FC>1)",
    "Down in Serum (log₂FC<−1)"
  ),
  Count = c(n_tested, n_strict, n_up_strict, n_down_strict)
)
print(summary_table)

# 6) volcano plot
plot_df <- as.data.frame(res)
plot_df$log10p <- -log10(plot_df$padj)
plot_df$strict <- with(plot_df, padj < 0.001 & abs(log2FoldChange) > 1)
plot_df$log10p[is.infinite(plot_df$log10p)] <- max(plot_df$log10p[is.finite(plot_df$log10p)])*1.1

p <- ggplot(plot_df, aes(log2FoldChange, log10p)) +
  geom_point(aes(color=strict), alpha=0.6, size=1.5) +
  scale_color_manual(values=c("FALSE"="grey70","TRUE"="firebrick")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.001), linetype="dashed") +
  labs(title="Volcano Plot (padj<0.001 & |log₂FC|>1)",
       x=expression(log[2]*" Fold Change"),
       y=expression(-log[10]*" adj. p-value")) +
  theme_minimal(base_size=14) +
  theme(legend.position="none")

print(p)
ggsave("volcano_strict.png", p, width=6, height=5)
