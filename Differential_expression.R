#!/usr/bin/env Rscript
# plot_top5_volcano.R

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
for (pkg in c("DESeq2","EnhancedVolcano","dplyr","tibble")) {
  if (!requireNamespace(pkg, quietly=TRUE))
    BiocManager::install(pkg, ask=FALSE)
}
suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(dplyr)
  library(tibble)
})

fc <- read.table("efaecium_featureCounts.txt", header=TRUE, sep="\t",
                 row.names=1, comment.char="#", check.names=FALSE)
countData <- fc[, -(1:6)]
colnames(countData) <- sub("\\.sorted\\.bam$", "", basename(colnames(countData)))

colData <- read.csv("sample_metadata.csv", row.names=1, stringsAsFactors=FALSE)
colData <- colData[colnames(countData), , drop=FALSE]
stopifnot(all(colnames(countData) == rownames(colData)))

dds <- DESeqDataSetFromMatrix(countData, colData, design=~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds)

res_df <- as.data.frame(res) %>%
  rownames_to_column("locus_tag") %>%
  mutate(
    pvalue = ifelse(pvalue == 0,
                    min(pvalue[pvalue > 0], na.rm=TRUE) * 1e-1,
                    pvalue)
  )

anno <- read.delim("efaecium.tsv", stringsAsFactors=FALSE) %>%
  select(locus_tag, gene) %>%
  mutate(label = ifelse(gene == "", locus_tag, gene))
res_df <- left_join(res_df, anno, by="locus_tag")

top5 <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%
  pull(label)

message("Top 5 genes: ", paste(top5, collapse=", "))

png("top5_volcano.png", width=6*300, height=6*300, res=300)
EnhancedVolcano(
  res_df,
  lab             = res_df$label,
  x               = 'log2FoldChange',
  y               = 'padj',
  selectLab       = top5,
  xlab            = bquote(~Log[2]~"(Serum/BH)"),
  ylab            = bquote(~-Log[10]~"(adj. p-value)"),
  title           = 'Volcano: Serum vs BH',
  pCutoff         = 0.001,
  FCcutoff        = 1,
  pointSize       = 2,
  labSize         = 4,
  col             = c('grey70','forestgreen','royalblue','red3'),
  colAlpha        = 0.8,
  drawConnectors  = TRUE,
  widthConnectors = 0.75,
  arrowheads      = TRUE,
  legendPosition  = 'right',
  legendLabSize   = 12,
  legendIconSize  = 4
)
dev.off()
