library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggrepel)
library(UpSetR)
library(clusterProfiler)
library(org.Dm.eg.db)

base_dir <- "/Volumes/My_Book/12926_Wang_Nanopore/bam"
results_dir <- file.path(base_dir, "deseq2_results")
figures_dir <- "figures"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Reconstruct deseq2 dds and vsd from raw counts and metadata
counts_raw <- read.csv("/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/Wang-Nanopore-Thorax/Head_Thorax_Comparison/13026-nanopore-Head_Thorax_Comparison-counts-all-samples.csv", row.names = 1, check.names = FALSE)
meta <- read.delim("/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/Wang-Nanopore-Thorax/Head_Thorax_Comparison/Higginson_metadata_batch_head_thorax.txt", row.names = 1)
meta$genotype <- factor(gsub(" (Head|Thorax)$", "", meta$Comparison), levels = c("WT", "G11A", "G146C"))
meta$tissue   <- factor(ifelse(grepl("Head", meta$Comparison), "head", "thorax"), levels = c("thorax", "head"))
meta$batch    <- factor(meta$batch)
counts_matrix <- as.matrix(counts_raw[, rownames(meta)])
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = meta, design = ~ genotype + tissue + genotype:tissue)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
save(dds, vsd, meta, file = file.path(results_dir, "deseq2_objects.RData"))

# PCA
pca_data <- plotPCA(vsd, intgroup = c("tissue", "genotype"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = genotype, shape = tissue)) +
  stat_ellipse(aes(fill = tissue, group = tissue), geom = "polygon", alpha = 0.12, color = NA, level = 0.95) +
  geom_point(size = 4, alpha = 0.85) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  scale_color_manual(values = c("WT" = "#2166AC", "G11A" = "#D6604D", "G146C" = "#4DAF4A")) +
  scale_fill_manual(values = c("head" = "#F46D43", "thorax" = "#74ADD1")) +
  scale_shape_manual(values = c("head" = 16, "thorax" = 17)) +
  theme_minimal(base_size = 14) +
  ggtitle("Principal Component Analysis (PCA)")

ggsave(file.path(figures_dir, "PCA_all_samples.pdf"), p_pca, width = 8, height = 6)
print(p_pca)

# Volcano
int_G11A_name  <- grep("G11A.*head|head.*G11A", resultsNames(dds), value = TRUE)
int_G146C_name <- grep("G146C.*head|head.*G146C", resultsNames(dds), value = TRUE)

res_G11A_head <- results(dds, contrast = list(c("genotype_G11A_vs_WT", int_G11A_name)), alpha = 0.05)
res_G146C_head <- results(dds, contrast = list(c("genotype_G146C_vs_WT", int_G146C_name)), alpha = 0.05)
res_G11A_thorax  <- results(dds, contrast = c("genotype", "G11A", "WT"), alpha = 0.05)
res_G146C_thorax <- results(dds, contrast = c("genotype", "G146C", "WT"), alpha = 0.05)

get_sig_genes <- function(res) {
  df <- as.data.frame(res)
  rownames(df)[!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1.5]
}

sig_G11A_head   <- get_sig_genes(res_G11A_head)
sig_G11A_thorax <- get_sig_genes(res_G11A_thorax)
sig_G146C_head   <- get_sig_genes(res_G146C_head)
sig_G146C_thorax <- get_sig_genes(res_G146C_thorax)

both_G11A  <- intersect(sig_G11A_head, sig_G11A_thorax)
both_G146C <- intersect(sig_G146C_head, sig_G146C_thorax)

make_volcano <- function(res, sig_this_tissue, both_genes, mutant_name, tissue_name) {
  df <- as.data.frame(res)
  df <- df[complete.cases(df), ]
  df$gene <- rownames(df)
  df$category <- "NS"
  passes_cutoff <- abs(df$log2FoldChange) > 1.5 & df$padj < 0.05
  df$category[passes_cutoff & df$gene %in% setdiff(sig_this_tissue, both_genes)] <- "Tissue-specific"
  df$category[passes_cutoff & df$gene %in% both_genes]      <- "Both tissues"
  df$category <- factor(df$category, levels = c("NS", "Tissue-specific", "Both tissues"))
  df <- df[order(df$category), ]
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = category, size = category)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("NS" = "grey80", "Tissue-specific" = "#D62728", "Both tissues" = "#2166AC"), name = "Category") +
    scale_size_manual(values = c("NS" = 0.6, "Tissue-specific" = 1.8, "Both tissues" = 1.8), guide = "none") +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    theme_minimal(base_size = 12) +
    labs(title = paste0(mutant_name, " vs WT - ", tissue_name), x = "log2FoldChange", y = "-log10(adjusted p-value)")
}

p_v_g11a_head   <- make_volcano(res_G11A_head,   sig_G11A_head,   both_G11A,  "G11A",  "Head")
p_v_g11a_thorax <- make_volcano(res_G11A_thorax, sig_G11A_thorax, both_G11A,  "G11A",  "Thorax")
p_v_g146c_head   <- make_volcano(res_G146C_head,   sig_G146C_head,   both_G146C, "G146C", "Head")
p_v_g146c_thorax <- make_volcano(res_G146C_thorax, sig_G146C_thorax, both_G146C, "G146C", "Thorax")

ggsave(file.path(figures_dir, "volcano_G11A_vs_WT_head.pdf"), p_v_g11a_head, width = 9, height = 7)
print(p_v_g11a_head)
ggsave(file.path(figures_dir, "volcano_G11A_vs_WT_thorax.pdf"), p_v_g11a_thorax, width = 9, height = 7)
print(p_v_g11a_thorax)
ggsave(file.path(figures_dir, "volcano_G146C_vs_WT_head.pdf"), p_v_g146c_head, width = 9, height = 7)
print(p_v_g146c_head)
ggsave(file.path(figures_dir, "volcano_G146C_vs_WT_thorax.pdf"), p_v_g146c_thorax, width = 9, height = 7)
print(p_v_g146c_thorax)

# Scatter
load_res_for_scatter <- function(res) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df %>% dplyr::select(gene_id, log2FoldChange, padj)
}

head_g11a_df <- load_res_for_scatter(res_G11A_head)
thorax_g11a_df <- load_res_for_scatter(res_G11A_thorax)
g11a_merged <- inner_join(head_g11a_df, thorax_g11a_df, by = "gene_id", suffix = c("_head", "_thorax")) %>% mutate(Genotype = "G11A")

head_g146c_df <- load_res_for_scatter(res_G146C_head)
thorax_g146c_df <- load_res_for_scatter(res_G146C_thorax)
g146c_merged <- inner_join(head_g146c_df, thorax_g146c_df, by = "gene_id", suffix = c("_head", "_thorax")) %>% mutate(Genotype = "G146C")

tissue_scatter_data <- bind_rows(g11a_merged, g146c_merged) %>%
  mutate(
    sig_head = !is.na(padj_head) & padj_head < 0.05 & abs(log2FoldChange_head) > 1.5,
    sig_thorax = !is.na(padj_thorax) & padj_thorax < 0.05 & abs(log2FoldChange_thorax) > 1.5,
    Significance = case_when(
      sig_head & sig_thorax ~ "Dysregulated in Both",
      sig_head & !sig_thorax ~ "Dysregulated in Head Only",
      !sig_head & sig_thorax ~ "Dysregulated in Thorax Only",
      TRUE ~ "Not Significant"
    )
  )

sig_counts_t <- tissue_scatter_data %>% count(Significance)
tissue_scatter_data <- tissue_scatter_data %>%
  left_join(sig_counts_t, by = "Significance") %>%
  mutate(Sig_Label = paste0(Significance, " (n = ", format(n, big.mark=","), ")"))

new_levels_t <- tissue_scatter_data %>% dplyr::select(Significance, Sig_Label) %>% distinct() %>% 
  arrange(match(Significance, c("Not Significant", "Dysregulated in Thorax Only", "Dysregulated in Both", "Dysregulated in Head Only")))
tissue_scatter_data$Significance <- factor(tissue_scatter_data$Significance, levels = c("Not Significant", "Dysregulated in Thorax Only", "Dysregulated in Both", "Dysregulated in Head Only"))

cor_data_t <- tissue_scatter_data %>%
  group_by(Genotype) %>%
  summarize(r2 = cor(log2FoldChange_thorax, log2FoldChange_head, use = "complete.obs")^2, .groups = 'drop') %>%
  mutate(label = sprintf("R² = %.2f", r2))

limit_t <- min(max(c(abs(tissue_scatter_data$log2FoldChange_head), abs(tissue_scatter_data$log2FoldChange_thorax)), na.rm=TRUE) * 1.05, 12)

p_scatter_tissue <- ggplot(tissue_scatter_data %>% arrange(Significance), aes(x = log2FoldChange_thorax, y = log2FoldChange_head)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", color = "black", alpha = 0.4) +
  geom_smooth(method = "lm", color = "#e66101", linewidth = 1, se = FALSE) +
  geom_point(aes(color = Significance), alpha = 0.7, size = 1.2) +
  geom_text(data = cor_data_t, aes(x = -limit_t * 0.9, y = limit_t * 0.9, label = label), hjust = 0, vjust = 1, fontface = "bold", size = 4.0, inherit.aes = FALSE) +
  facet_wrap(~ Genotype, ncol = 2) +
  scale_color_manual(
    values = c("Not Significant" = "grey85", "Dysregulated in Thorax Only" = "#4575b4", "Dysregulated in Both" = "#542788", "Dysregulated in Head Only" = "#d73027"),
    labels = setNames(new_levels_t$Sig_Label, new_levels_t$Significance)
  ) +
  coord_fixed(xlim = c(-limit_t, limit_t), ylim = c(-limit_t, limit_t)) +
  theme_bw(base_size = 12) +
  labs(title = "Tissue Susceptibility to Exosome Dysfunction", x = "Thorax LFC", y = "Head LFC", color = "Significance")

ggsave(file.path(figures_dir, "tissue_susceptibility_scatter_with_counts.pdf"), p_scatter_tissue, width = 10, height = 6.5)
print(p_scatter_tissue)

head_merged <- inner_join(head_g11a_df, head_g146c_df, by = "gene_id", suffix = c("_G11A", "_G146C")) %>% mutate(Tissue = "Head")
thorax_merged <- inner_join(thorax_g11a_df, thorax_g146c_df, by = "gene_id", suffix = c("_G11A", "_G146C")) %>% mutate(Tissue = "Thorax")

genotype_scatter_data <- bind_rows(head_merged, thorax_merged) %>%
  mutate(
    sig_g11a = !is.na(padj_G11A) & padj_G11A < 0.05 & abs(log2FoldChange_G11A) > 1.5,
    sig_g146c = !is.na(padj_G146C) & padj_G146C < 0.05 & abs(log2FoldChange_G146C) > 1.5,
    Significance = case_when(
      sig_g11a & sig_g146c ~ "Dysregulated in Both",
      sig_g11a & !sig_g11a ~ "Dysregulated in G11A Only",
      !sig_g11a & sig_g146c ~ "Dysregulated in G146C Only",
      TRUE ~ "Not Significant"
    )
  )

sig_counts_g <- genotype_scatter_data %>% count(Significance)
genotype_scatter_data <- genotype_scatter_data %>%
  left_join(sig_counts_g, by = "Significance") %>%
  mutate(Sig_Label = paste0(Significance, " (n = ", format(n, big.mark=","), ")"))

new_levels_g <- genotype_scatter_data %>% dplyr::select(Significance, Sig_Label) %>% distinct() %>% 
  arrange(match(Significance, c("Not Significant", "Dysregulated in G11A Only", "Dysregulated in G146C Only", "Dysregulated in Both")))
genotype_scatter_data$Significance <- factor(genotype_scatter_data$Significance, levels = c("Not Significant", "Dysregulated in G11A Only", "Dysregulated in G146C Only", "Dysregulated in Both"))

cor_data_g <- genotype_scatter_data %>%
  group_by(Tissue) %>%
  summarize(r2 = cor(log2FoldChange_G11A, log2FoldChange_G146C, use = "complete.obs")^2, .groups = 'drop') %>%
  mutate(label = sprintf("R² = %.2f", r2))

limit_g <- min(max(c(abs(genotype_scatter_data$log2FoldChange_G11A), abs(genotype_scatter_data$log2FoldChange_G146C)), na.rm=TRUE) * 1.05, 12)

p_scatter_genotype <- ggplot(genotype_scatter_data %>% arrange(Significance), aes(x = log2FoldChange_G11A, y = log2FoldChange_G146C)) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", color = "black", alpha = 0.4) +
  geom_smooth(method = "lm", color = "#e66101", linewidth = 1, se = FALSE) +
  geom_point(aes(color = Significance), alpha = 0.7, size = 1.2) +
  geom_text(data = cor_data_g, aes(x = -limit_g * 0.9, y = limit_g * 0.9, label = label), hjust = 0, vjust = 1, fontface = "bold", size = 4.0, inherit.aes = FALSE) +
  facet_wrap(~ Tissue, ncol = 2) +
  scale_color_manual(
    values = c("Not Significant" = "grey85", "Dysregulated in G11A Only" = "#4575b4", "Dysregulated in G146C Only" = "#d73027", "Dysregulated in Both" = "#542788"),
    labels = setNames(new_levels_g$Sig_Label, new_levels_g$Significance)
  ) +
  coord_fixed(xlim = c(-limit_g, limit_g), ylim = c(-limit_g, limit_g)) +
  theme_bw(base_size = 12) +
  labs(title = "Concordance Between Exosome Mutations", x = "G11A LFC", y = "G146C LFC", color = "Significance")

ggsave(file.path(figures_dir, "genotype_scatter_correlation_counts.pdf"), p_scatter_genotype, width = 10, height = 6.5)
print(p_scatter_genotype)

# UpSet
get_sig_for_upset <- function(res, direction = "up") {
  df <- as.data.frame(res)
  df <- df[!is.na(df$padj), ]
  if (direction == "up") rownames(df[df$padj < 0.05 & df$log2FoldChange > 1.5, ])
  else rownames(df[df$padj < 0.05 & df$log2FoldChange < -1.5, ])
}

sets_up <- list("G11A Head" = get_sig_for_upset(res_G11A_head, "up"), "G11A Thorax" = get_sig_for_upset(res_G11A_thorax, "up"), "G146C Head" = get_sig_for_upset(res_G146C_head, "up"), "G146C Thorax" = get_sig_for_upset(res_G146C_thorax, "up"))
sets_down <- list("G11A Head" = get_sig_for_upset(res_G11A_head, "down"), "G11A Thorax" = get_sig_for_upset(res_G11A_thorax, "down"), "G146C Head" = get_sig_for_upset(res_G146C_head, "down"), "G146C Thorax" = get_sig_for_upset(res_G146C_thorax, "down"))

pdf(file.path(figures_dir, "upset_upregulated.pdf"), width = 10, height = 6, onefile = FALSE)
upset(fromList(sets_up), order.by = "freq", keep.order = TRUE, nsets = 4, sets.bar.color = c("#0072B2", "#009E73", "#D55E00", "#56B4E9"))
dev.off()
upset(fromList(sets_up), order.by = "freq", keep.order = TRUE, nsets = 4, sets.bar.color = c("#0072B2", "#009E73", "#D55E00", "#56B4E9"))

pdf(file.path(figures_dir, "upset_downregulated.pdf"), width = 10, height = 6, onefile = FALSE)
upset(fromList(sets_down), order.by = "freq", keep.order = TRUE, nsets = 4, sets.bar.color = c("#0072B2", "#009E73", "#D55E00", "#56B4E9"))
dev.off()
upset(fromList(sets_down), order.by = "freq", keep.order = TRUE, nsets = 4, sets.bar.color = c("#0072B2", "#009E73", "#D55E00", "#56B4E9"))

# Antisense RNA Analysis
as_csv_path <- file.path(results_dir, "as_RNA_analysis", "asRNA_sense_overlap_DEG_results.csv")
overlap_base <- read.csv(as_csv_path, stringsAsFactors = FALSE)
class_df <- read.csv(file.path(results_dir, "as_RNA_analysis", "asRNA_TSS_classifications.csv"), stringsAsFactors = FALSE)

neuronal_mapping <- data.frame(
  asRNA_symbol = c("asRNA.CR45281", "asRNA.CR46354", "asRNA.CR45009", "asRNA.CR43959", "asRNA.CR45124", "asRNA.CR46484", "asRNA.CR44031"),
  label = c("asRNA.CR45281 / Alp6\n(Alkaline Phosphatase)", "asRNA.CR46354 / side-IV\n(Motor Axon Guidance)", "asRNA.CR45009 / FucTC\n(Fucosyltransferase C)", "asRNA.CR43959 / Hk\n(Neuronal Potassium Channel)", "asRNA.CR45124 / PlexB\n(Axon Guidance Receptor)", "asRNA.CR46484 / CG11396\n(Predicted Protein Kinase)", "asRNA.CR44031 / Mpv17\n(Mitochondrial Inner Membrane)")
)

panel_a_df <- overlap_base %>%
  inner_join(neuronal_mapping, by = "asRNA_symbol") %>%
  dplyr::select(label, G11A_Head_asRNA_LFC, G11A_Head_sense_LFC, G146C_Head_asRNA_LFC, G146C_Head_sense_LFC, G11A_Thorax_asRNA_LFC, G11A_Thorax_sense_LFC, G146C_Thorax_asRNA_LFC, G146C_Thorax_sense_LFC) %>%
  pivot_longer(cols = -label, names_to = c("Condition", "TranscriptType"), names_pattern = "([A-Z0-9a-z_]+)_([a-zA-Z]+)_LFC", values_to = "LFC")

panel_a_df$Condition <- factor(panel_a_df$Condition, levels = c("G11A_Head", "G146C_Head", "G11A_Thorax", "G146C_Thorax"), labels = c("G11A Head", "G146C Head", "G11A Thorax", "G146C Thorax"))
panel_a_df$TranscriptType <- factor(panel_a_df$TranscriptType, levels = c("asRNA", "sense"), labels = c("Antisense RNA", "Sense Host Gene"))
panel_a_df$label <- factor(panel_a_df$label, levels = rev(neuronal_mapping$label))

p_2b <- ggplot(panel_a_df, aes(x = LFC, y = label, fill = TranscriptType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Condition, ncol = 4) +
  theme_bw() +
  scale_fill_manual(values = c("Antisense RNA" = "#E66151", "Sense Host Gene" = "#4393C3")) +
  labs(title = "Coordinated Regulation of Key Upregulated Targets", x = "log2FoldChange", y = NULL, fill = "Transcript Type")

ggsave(file.path(figures_dir, "panel_A_coordinated_regulation.pdf"), p_2b, width = 12, height = 5)
print(p_2b)

sig_asRNAs <- overlap_base %>%
  dplyr::filter((G11A_Head_asRNA_pvalue < 0.05 & G11A_Head_asRNA_LFC > 0) | (G146C_Head_asRNA_pvalue < 0.05 & G146C_Head_asRNA_LFC > 0) | (G11A_Thorax_asRNA_pvalue < 0.05 & G11A_Thorax_asRNA_LFC > 0) | (G146C_Thorax_asRNA_pvalue < 0.05 & G146C_Thorax_asRNA_LFC > 0)) %>%
  dplyr::filter(!is.na(asRNA_symbol))

panel_b_list <- list()
conditions_c <- c("G11A_Head", "G146C_Head", "G11A_Thorax", "G146C_Thorax")
for (cond in conditions_c) {
  panel_b_list[[cond]] <- sig_asRNAs %>%
    mutate(Condition = cond, as_LFC = .data[[paste0(cond, "_asRNA_LFC")]], sense_LFC = .data[[paste0(cond, "_sense_LFC")]], sense_pval = .data[[paste0(cond, "_sense_pvalue")]]) %>%
    mutate(Response = case_when(is.na(sense_LFC) | is.na(sense_pval) ~ "Unchanged / Non-Significant", sense_pval < 0.05 & sense_LFC > 0 ~ "Coordinated Upregulation\n(Sense Up)", sense_pval < 0.05 & sense_LFC < 0 ~ "Opposite Regulation\n(Sense Down)", TRUE ~ "Unchanged / Non-Significant")) %>%
    dplyr::select(asRNA_symbol, Condition, Response)
}

panel_b_df <- do.call(rbind, panel_b_list)
panel_b_df$Condition <- factor(panel_b_df$Condition, levels = conditions_c, labels = c("G11A Head", "G146C Head", "G11A Thorax", "G146C Thorax"))
panel_b_df$Response <- factor(panel_b_df$Response, levels = c("Coordinated Upregulation\n(Sense Up)", "Unchanged / Non-Significant", "Opposite Regulation\n(Sense Down)"))

p_2c <- ggplot(panel_b_df, aes(x = Condition, fill = Response)) +
  geom_bar(position = "fill", width = 0.55) +
  theme_bw() +
  scale_fill_manual(values = c("Coordinated Upregulation\n(Sense Up)" = "#E66151", "Unchanged / Non-Significant" = "#E0E0E0", "Opposite Regulation\n(Sense Down)" = "#4393C3")) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportion of Overlapping Sense Gene Responses", x = "Condition", y = "Percentage of Targets", fill = "Sense Gene Response")

ggsave(file.path(figures_dir, "panel_C_overlapping_responses.pdf"), p_2c, width = 7, height = 6.5)
print(p_2c)

panel_c_df <- sig_asRNAs %>%
  dplyr::select(asRNA_symbol, G11A_Head_asRNA_LFC, G146C_Head_asRNA_LFC, G11A_Thorax_asRNA_LFC, G146C_Thorax_asRNA_LFC) %>%
  pivot_longer(cols = -asRNA_symbol, names_to = "Condition", values_to = "LFC") %>%
  mutate(Condition = gsub("_asRNA_LFC", "", Condition), AbsLFC = abs(LFC))

panel_c_df$Condition <- factor(panel_c_df$Condition, levels = conditions_c, labels = c("G11A Head", "G146C Head", "G11A Thorax", "G146C Thorax"))

p_2d <- ggplot(panel_c_df, aes(x = Condition, y = AbsLFC, fill = Condition)) +
  geom_violin(alpha = 0.4) + geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA) + geom_jitter(width = 0.12, size = 1.8, alpha = 0.85, aes(color = Condition)) +
  theme_bw() + scale_fill_manual(values = c("G11A Head" = "#4F46E5", "G146C Head" = "#06B6D4", "G11A Thorax" = "#10B981", "G146C Thorax" = "#F59E0B")) +
  scale_color_manual(values = c("G11A Head" = "#4F46E5", "G146C Head" = "#06B6D4", "G11A Thorax" = "#10B981", "G146C Thorax" = "#F59E0B")) +
  labs(title = "Dysregulation Severity (Significant asRNAs)", x = "Condition", y = "asRNA |log2FoldChange| vs. WT") + theme(legend.position = "none")

ggsave(file.path(figures_dir, "panel_D_dysregulation_severity.pdf"), p_2d, width = 7, height = 6.5)
print(p_2d)

# Target Category Upregulation (Fisher's exact barplots with brackets)
mito_entrez <- tryCatch(get("GO:0005739", org.Dm.egGO2ALLEGS), error = function(e) NULL)
mito_symbols <- if (!is.null(mito_entrez)) unique(bitr(mito_entrez, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Dm.eg.db)$SYMBOL) else c()
neuro_entrez <- unique(c(tryCatch(get("GO:0007399", org.Dm.egGO2ALLEGS), error = function(e) NULL), tryCatch(get("GO:0045202", org.Dm.egGO2ALLEGS), error = function(e) NULL)))
neuro_symbols <- if (length(neuro_entrez) > 0) unique(bitr(neuro_entrez, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Dm.eg.db)$SYMBOL) else c()

as_data <- overlap_base %>% mutate(is_mito = sense_symbol %in% mito_symbols, is_neuro = sense_symbol %in% neuro_symbols)

get_genotype_summary <- function(data, genotype) {
  df_flags <- data %>%
    mutate(up_Head = coalesce(.data[[paste0(genotype, "_Head_asRNA_pvalue")]] < 0.05 & .data[[paste0(genotype, "_Head_asRNA_LFC")]] > 0, FALSE), up_Thorax = coalesce(.data[[paste0(genotype, "_Thorax_asRNA_pvalue")]] < 0.05 & .data[[paste0(genotype, "_Thorax_asRNA_LFC")]] > 0, FALSE))
  bind_rows(list(
    data.frame(Category = "Neuronal", N = sum(df_flags$is_neuro), Tissue = c("Head (Neuro)", "Thorax (Muscle)"), Upregulated = c(sum(df_flags$up_Head[df_flags$is_neuro]), sum(df_flags$up_Thorax[df_flags$is_neuro]))),
    data.frame(Category = "Non-neuronal", N = sum(!df_flags$is_neuro), Tissue = c("Head (Neuro)", "Thorax (Muscle)"), Upregulated = c(sum(df_flags$up_Head[!df_flags$is_neuro]), sum(df_flags$up_Thorax[!df_flags$is_neuro]))),
    data.frame(Category = "Mitochondrial", N = sum(df_flags$is_mito), Tissue = c("Head (Neuro)", "Thorax (Muscle)"), Upregulated = c(sum(df_flags$up_Head[df_flags$is_mito]), sum(df_flags$up_Thorax[df_flags$is_mito]))),
    data.frame(Category = "Non-mitochondrial", N = sum(!df_flags$is_mito), Tissue = c("Head (Neuro)", "Thorax (Muscle)"), Upregulated = c(sum(df_flags$up_Head[!df_flags$is_mito]), sum(df_flags$up_Thorax[!df_flags$is_mito])))
  )) %>% mutate(Percentage = Upregulated / N * 100, Label = sprintf("%.1f%%\n(%d/%d)", Percentage, Upregulated, N), Category_Label = sprintf("%s\n(n = %d)", Category, N))
}

draw_brackets <- function(p, brackets_df, summary_df) {
  if (is.null(brackets_df) || nrow(brackets_df) == 0) return(p)
  tick_h <- max(summary_df$Percentage) * 0.025
  for (i in 1:nrow(brackets_df)) {
    b <- brackets_df[i, ]
    p <- p + annotate("segment", x = b$x1, xend = b$x2, y = b$y, yend = b$y, color = "black", linewidth = 0.4) +
      annotate("segment", x = b$x1, xend = b$x1, y = b$y, yend = b$y - tick_h, color = "black", linewidth = 0.4) +
      annotate("segment", x = b$x2, xend = b$x2, y = b$y, yend = b$y - tick_h, color = "black", linewidth = 0.4) +
      annotate("text", x = (b$x1 + b$x2)/2, y = b$y + tick_h * 0.4, label = b$label, size = 2.9, fontface = ifelse(grepl("\\*", b$label), "bold", "italic"), vjust = 0)
  }
  return(p)
}

create_enrichment_plot <- function(summary_df, genotype_name, subtitle_txt, brackets_df = NULL, fill_colors) {
  summary_df$Category_Label <- factor(summary_df$Category_Label, levels = unique(summary_df$Category_Label))
  summary_df$Tissue <- factor(summary_df$Tissue, levels = c("Head (Neuro)", "Thorax (Muscle)"))
  p <- ggplot(summary_df, aes(x = Category_Label, y = Percentage, fill = Tissue)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    geom_text(aes(label = Label), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.0, fontface = "bold") +
    theme_bw() + scale_fill_manual(values = fill_colors) +
    labs(title = genotype_name, subtitle = subtitle_txt, x = NULL, y = "Percentage Upregulated (%)") +
    scale_y_continuous(limits = c(0, max(summary_df$Percentage, if(!is.null(brackets_df)) max(brackets_df$y) else 0) * 1.15))
  draw_brackets(p, brackets_df, summary_df)
}

p_g11a <- create_enrichment_plot(
  get_genotype_summary(as_data, "G11A"), 
  "G11A Mutant (Severe)", 
  "Upregulated defined as pvalue < 0.05 & log2FoldChange > 0 compared to WT",
  brackets_df = data.frame(x1 = c(1 - 0.2, 1 + 0.2, 3 - 0.2, 3 + 0.2), x2 = c(2 - 0.2, 2 + 0.2, 4 - 0.2, 4 + 0.2), y = c(16.8, 9.8, 9.0, 12.0), label = c("Head: p = 0.0059**", "Thorax: p = 0.170 (n.s.)", "Head: p = 0.711 (n.s.)", "Thorax: p = 0.353 (n.s.)")),
  fill_colors = c("Head (Neuro)" = "#0072B2", "Thorax (Muscle)" = "#D55E00")
)

p_g146c <- create_enrichment_plot(
  get_genotype_summary(as_data, "G146C"), 
  "G146C Mutant (Mild)", 
  "Upregulated defined as pvalue < 0.05 & log2FoldChange > 0 compared to WT",
  brackets_df = data.frame(x1 = c(1 - 0.2, 1 + 0.2, 3 - 0.2, 3 + 0.2), x2 = c(2 - 0.2, 2 + 0.2, 4 - 0.2, 4 + 0.2), y = c(3.6, 1.2, 2.4, 1.6), label = c("Head: p = 0.296 (n.s.)", "Thorax: p = 1.000 (n.s.)", "Head: p = 1.000 (n.s.)", "Thorax: p = 1.000 (n.s.)")),
  fill_colors = c("Head (Neuro)" = "#0072B2", "Thorax (Muscle)" = "#D55E00")
)

ggsave(file.path(figures_dir, "asRNA_upregulation_G11A.pdf"), p_g11a, width = 8.5, height = 6)
print(p_g11a)
ggsave(file.path(figures_dir, "asRNA_upregulation_G146C.pdf"), p_g146c, width = 8.5, height = 6)
print(p_g146c)

# TSS Classification Pie Charts (Panel B)
plot_list_b <- list()
for (cond in conditions_c) {
  cond_sig <- overlap_base %>% dplyr::filter(.data[[paste0(cond, "_asRNA_pvalue")]] < 0.05 & .data[[paste0(cond, "_asRNA_LFC")]] > 0) %>% dplyr::select(asRNA_symbol)
  if (nrow(cond_sig) > 0) plot_list_b[[cond]] <- class_df %>% dplyr::filter(asRNA_Symbol %in% cond_sig$asRNA_symbol) %>% mutate(Condition = cond)
}

class_plot_df <- bind_rows(plot_list_b) %>% mutate(Genotype = ifelse(grepl("G11A", Condition), "G11A", "G146C"), Tissue = ifelse(grepl("Head", Condition), "Head", "Thorax"))
totals_df <- class_plot_df %>% group_by(Genotype, Tissue) %>% summarize(total_count = n(), .groups = 'drop')
class_counts <- class_plot_df %>% group_by(Genotype, Tissue, Classification) %>% summarize(count = n(), .groups = 'drop')

pie_data_complete <- expand.grid(Genotype = c("G11A", "G146C"), Tissue = c("Head", "Thorax"), Classification = c("PROMPT (Divergent Promoter)", "3' Overlapping / Read-through", "Gene-Body Overlapping"), stringsAsFactors = FALSE) %>%
  left_join(class_counts %>% left_join(totals_df, by = c("Genotype", "Tissue")) %>% mutate(percent = count / total_count, label = paste0(count, "\n(", round(percent * 100, 0), "%)")), by = c("Genotype", "Tissue", "Classification")) %>%
  mutate(count = coalesce(count, 0), percent = coalesce(percent, 0), label = coalesce(label, "")) %>%
  dplyr::select(-any_of("total_count")) %>% left_join(totals_df, by = c("Genotype", "Tissue")) %>% mutate(total_count = coalesce(total_count, 0), facet_label = paste0(Genotype, " ", Tissue, " (N = ", total_count, ")"))

facet_order <- sprintf("%s %s (N = %d)", c("G11A", "G146C", "G11A", "G146C"), c("Head", "Head", "Thorax", "Thorax"), sapply(list(c("G11A","Head"), c("G146C","Head"), c("G11A","Thorax"), c("G146C","Thorax")), function(x) coalesce(totals_df$total_count[totals_df$Genotype == x[1] & totals_df$Tissue == x[2]], 0L)))
pie_data_complete$facet_label <- factor(pie_data_complete$facet_label, levels = facet_order)
pie_data_complete$Classification <- factor(pie_data_complete$Classification, levels = c("PROMPT (Divergent Promoter)", "3' Overlapping / Read-through", "Gene-Body Overlapping"))

p_pie <- ggplot(pie_data_complete %>% dplyr::filter(total_count > 0), aes(x = 0.5, y = percent, fill = Classification)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) + coord_polar(theta = "y") + facet_wrap(~ facet_label, ncol = 2) + theme_void() +
  scale_fill_manual(values = c("PROMPT (Divergent Promoter)" = "#2CA02C", "3' Overlapping / Read-through" = "#D62728", "Gene-Body Overlapping" = "#1F77B4")) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold", color = "white") +
  labs(title = "Spatial Classification of Significant asRNAs", fill = "TSS Classification")

ggsave(file.path(figures_dir, "panel_B_spatial_classification.pdf"), p_pie, width = 8, height = 8.5)
print(p_pie)

# Composite summary layout
p_a_comb <- p_2b + theme(plot.title = element_text(size = 10, face = "bold"))
p_b_comb <- p_pie + theme(
  plot.title = element_text(size = 10, face = "bold"),
  strip.text = element_text(size = 8, face = "bold")
)
p_c_comb <- p_2c + theme(plot.title = element_text(size = 10, face = "bold"))
p_d_comb <- p_2d + theme(plot.title = element_text(size = 10, face = "bold"))

layout <- (p_a_comb) / (p_b_comb + ((p_c_comb / p_d_comb) + plot_layout(heights = c(1, 1)))) + 
  plot_layout(heights = c(1, 2.2), widths = c(1, 1.2)) +
  plot_annotation(
    title = "Characterization of Differentially Expressed Antisense RNAs (asRNAs)",
    subtitle = "Integrative analysis of genomic coordinates, TSS mapping, sense gene coordination, and exosome subunit pathology",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5, face = "italic")
    )
  )

ggsave(file.path(figures_dir, "asRNA_findings_summary_composite.pdf"), layout, width = 15, height = 13)
print(layout)

# Mitochondrial Dysregulation
# GSEA BP
mito_entrez <- AnnotationDbi::get("GO:0005739", org.Dm.egGO2ALLEGS)
nemg_symbols <- AnnotationDbi::select(org.Dm.eg.db, keys = mito_entrez, columns = c("SYMBOL", "CHR"), keytype = "ENTREZID") %>% dplyr::filter(!CHR %in% c("MT", "M", "mitochondrion_genome", "Un", NA)) %>% pull(SYMBOL) %>% unique()

conditions_gsea <- list("G11A_Head" = res_G11A_head, "G146C_Head" = res_G146C_head, "G11A_Thorax" = res_G11A_thorax, "G146C_Thorax" = res_G146C_thorax)
gsea_results_list <- list()

for (cond_name in names(conditions_gsea)) {
  df <- as.data.frame(conditions_gsea[[cond_name]])
  df$SYMBOL <- rownames(df)
  df <- df %>% dplyr::filter(!is.na(stat) & SYMBOL %in% nemg_symbols)
  if (nrow(df) == 0) next
  df_mapped <- merge(df, suppressMessages(bitr(df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dm.eg.db)), by = "SYMBOL") %>% group_by(ENTREZID) %>% slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  gene_list <- sort(setNames(df_mapped$stat, df_mapped$ENTREZID), decreasing = TRUE)
  
  gsea_res <- tryCatch(gseGO(geneList = gene_list, OrgDb = org.Dm.eg.db, ont = "BP", minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE), error = function(e) NULL)
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) gsea_results_list[[cond_name]] <- as.data.frame(setReadable(gsea_res, OrgDb = org.Dm.eg.db, keyType = "ENTREZID")) %>% mutate(Condition = cond_name)
}

if (length(gsea_results_list) > 0) {
  all_gsea <- bind_rows(gsea_results_list)
  plot_gsea_df <- all_gsea %>% dplyr::filter(Description %in% unique(c(all_gsea %>% dplyr::filter(NES > 0) %>% group_by(Condition) %>% slice_max(order_by = NES, n = 3, with_ties = FALSE) %>% pull(Description), all_gsea %>% dplyr::filter(NES < 0) %>% group_by(Condition) %>% slice_min(order_by = NES, n = 3, with_ties = FALSE) %>% pull(Description))))
  plot_gsea_df$Condition <- factor(plot_gsea_df$Condition, levels = c("G11A_Head", "G146C_Head", "G11A_Thorax", "G146C_Thorax"), labels = c("G11A Head", "G146C Head", "G11A Thorax", "G146C Thorax"))
  plot_gsea_df$Direction <- factor(ifelse(plot_gsea_df$NES > 0, "Upregulated Pathways (Positive NES)", "Downregulated Pathways (Negative NES)"))
  plot_gsea_df$padj_log <- pmin(-log10(plot_gsea_df$p.adjust), 10)
  
  p_gsea <- ggplot(plot_gsea_df, aes(x = Condition, y = Description)) +
    geom_point(aes(size = padj_log, color = NES)) +
    scale_color_gradient2(low = "#2166AC", mid = "white", high = "#D62728", midpoint = 0, name = "Normalized Enrichment\nScore (NES)") +
    scale_size_continuous(range = c(3, 9), name = "Significance\n-log10(padj)") +
    facet_grid(Direction ~ ., scales = "free_y", space = "free_y") + theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    labs(title = "NEMGs GSEA: Biological Process Dysregulation", x = NULL, y = NULL)
  
  ggsave(file.path(figures_dir, "gsea_combined_dotplot_GSEA_BP_NEMGs.pdf"), p_gsea, width = 12, height = 7)
  print(p_gsea)
}

dynamics_genes <- c("Drp1", "Fis1", "Marf", "Opa1", "Pink1", "park")
nemg_no_dynamics <- setdiff(nemg_symbols, dynamics_genes)

extract_res_for_heatmap <- function(res, genotype, tissue) {
  as.data.frame(res) %>% mutate(gene = rownames(.), Condition = paste(genotype, tissue))
}

combined_res_h <- bind_rows(extract_res_for_heatmap(res_G11A_head, "G11A", "Head"), extract_res_for_heatmap(res_G146C_head, "G146C", "Head"), extract_res_for_heatmap(res_G11A_thorax, "G11A", "Thorax"), extract_res_for_heatmap(res_G146C_thorax, "G146C", "Thorax")) %>% dplyr::filter(!is.na(log2FoldChange))
top_10_nemgs <- combined_res_h %>% dplyr::filter(gene %in% nemg_no_dynamics & !grepl("^CG[0-9]", gene, ignore.case = TRUE) & padj < 0.05 & log2FoldChange > 0) %>% group_by(gene) %>% summarize(max_lfc = max(log2FoldChange), .groups = 'drop') %>% slice_max(order_by = max_lfc, n = 10) %>% pull(gene)

plot_genes_df <- combined_res_h %>% dplyr::filter(gene %in% c(top_10_nemgs, dynamics_genes)) %>% mutate(Category = ifelse(gene %in% dynamics_genes, "Mitochondrial Dynamics", "Top 10 Upregulated Structural & Metabolic NEMGs"))
plot_genes_df$Condition <- factor(plot_genes_df$Condition, levels = c("G11A Head", "G146C Head", "G11A Thorax", "G146C Thorax"))
gene_order <- plot_genes_df %>% group_by(gene) %>% summarize(mean_lfc = mean(log2FoldChange), Category = first(Category), .groups = 'drop') %>% arrange(Category, mean_lfc) %>% pull(gene)
plot_genes_df$gene <- factor(plot_genes_df$gene, levels = gene_order)

p_heatmap <- ggplot(plot_genes_df, aes(x = Condition, y = gene, fill = log2FoldChange)) +
  geom_tile(color = "white") + facet_wrap(~ Category, nrow = 1, scales = "free_y") +
  scale_fill_gradient2(low = "#3B82F6", mid = "white", high = "#EF4444", midpoint = 0, limit = c(-2.0, 3.5), name = expression("log"[2]*"FC")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9), axis.text.y = element_text(size = 9, face = "italic"), strip.text = element_text(face = "bold", color = "white"), strip.background = element_rect(fill = "#374151"))

ggsave(file.path(figures_dir, "NEMG_and_Dynamics_Heatmap.pdf"), p_heatmap, width = 10, height = 5)
print(p_heatmap)

# Individual NEMG volcano plots
make_NEMG_volcano <- function(res, genotype, tissue) {
  df <- as.data.frame(res) %>% mutate(gene = rownames(.)) %>% dplyr::filter(gene %in% c(nemg_no_dynamics, dynamics_genes) & !is.na(log2FoldChange) & !is.na(pvalue)) %>%
    mutate(is_sig = pvalue < 0.05 & abs(log2FoldChange) > 1.0, category = case_when(!is_sig ~ "NS", gene %in% dynamics_genes ~ "Mitochondrial Dynamics", TRUE ~ "Mitochondrial Structural & Metabolic NEMGs"))
  df$category <- factor(df$category, levels = c("NS", "Mitochondrial Structural & Metabolic NEMGs", "Mitochondrial Dynamics"))
  
  labeled_genes <- unique(c(df %>% dplyr::filter(category != "NS" & !grepl("^CG[0-9]", gene)) %>% slice_min(order_by = pvalue, n = 5, with_ties = FALSE) %>% pull(gene), intersect(df$gene, dynamics_genes)))
  limit_val <- max(2.5, ceiling(max(abs(df$log2FoldChange)) * 10) / 10)
  
  p <- ggplot(df[order(df$category), ], aes(x = log2FoldChange, y = -log10(pvalue), color = category, size = category)) +
    geom_point(alpha = 0.6) + scale_color_manual(values = c("NS" = "grey80", "Mitochondrial Structural & Metabolic NEMGs" = "#D62728", "Mitochondrial Dynamics" = "#2166AC"), name = "Category") +
    scale_size_manual(values = c("NS" = 0.6, "Mitochondrial Structural & Metabolic NEMGs" = 1.8, "Mitochondrial Dynamics" = 1.8), guide = "none") +
    geom_point(data = df %>% dplyr::filter(gene %in% labeled_genes), aes(x = log2FoldChange, y = -log10(pvalue)), shape = 1, color = "black", size = 2.8, stroke = 1.0, inherit.aes = FALSE) +
    geom_text_repel(data = df %>% dplyr::filter(gene %in% labeled_genes), aes(label = gene), size = 3.2, fontface = "bold", color = "black", max.overlaps = 15) +
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) + xlim(-limit_val, limit_val) +
    theme_minimal(base_size = 11) + theme(legend.position = "right") + labs(title = paste(genotype, "vs WT -", tissue), x = "log2FC", y = "-log10(p-value)")
  
  fname <- paste0("volcano_NEMGs_", genotype, "_vs_WT_", tolower(tissue), ".pdf")
  ggsave(file.path(figures_dir, fname), p, width = 9, height = 7)
  print(p)
}

make_NEMG_volcano(res_G11A_head, "G11A", "Head")
make_NEMG_volcano(res_G11A_thorax, "G11A", "Thorax")
make_NEMG_volcano(res_G146C_head, "G146C", "Head")
make_NEMG_volcano(res_G146C_thorax, "G146C", "Thorax")

# Mito Dynamics (Separate Head and Thorax plots)
make_dynamics_plot <- function(res_g11a, res_g146c, tissue_name) {
  df_g11a <- as.data.frame(res_g11a) %>% mutate(gene = rownames(.), Mutant = "G11A")
  df_g146c <- as.data.frame(res_g146c) %>% mutate(gene = rownames(.), Mutant = "G146C")
  df_all <- bind_rows(df_g11a, df_g146c) %>% dplyr::filter(gene %in% dynamics_genes) %>%
    mutate(Category = case_when(
      gene %in% c("Drp1", "Fis1") ~ "Mitochondrial Fission",
      gene %in% c("Marf", "Opa1") ~ "Mitochondrial Fusion",
      TRUE ~ "Mitophagy / Quality Control"
    ))
  df_all$gene <- factor(df_all$gene, levels = dynamics_genes)
  df_all$Category <- factor(df_all$Category, levels = c("Mitochondrial Fission", "Mitochondrial Fusion", "Mitophagy / Quality Control"))
  
  ggplot(df_all, aes(x = gene, y = log2FoldChange, fill = Mutant)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text(aes(label = ifelse(!is.na(padj) & padj < 0.05, "*", ""), y = log2FoldChange + sign(log2FoldChange) * 0.1), position = position_dodge(width = 0.8), size = 6, vjust = 0.5) +
    facet_grid(~ Category, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c("G11A" = "#009E73", "G146C" = "#56B4E9")) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 12, color = "white"),
      strip.background = element_rect(fill = "grey40", color = NA),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = paste("Mitochondrial Dynamics in Mutant", tissue_name),
      x = "Key Mitochondrial Gene",
      y = expression("Fold Change (" ~ log[2] ~ "FC) vs WT"),
      fill = "Genotype"
    )
}

p_dyn_head <- make_dynamics_plot(res_G11A_head, res_G146C_head, "Head")
p_dyn_thorax <- make_dynamics_plot(res_G11A_thorax, res_G146C_thorax, "Thorax")

ggsave(file.path(figures_dir, "mitochondria_dynamics_head.pdf"), p_dyn_head, width = 10, height = 6)
print(p_dyn_head)
ggsave(file.path(figures_dir, "mitochondria_dynamics_thorax.pdf"), p_dyn_thorax, width = 10, height = 6)
print(p_dyn_thorax)

# Barplots degradosome & thorax enriched
extract_gene_data_for_barplot <- function(res_obj, cond_name, genes) {
  present_genes <- intersect(genes, rownames(res_obj))
  as.data.frame(res_obj[present_genes, ]) %>% mutate(Gene = rownames(.), Condition = cond_name) %>% dplyr::select(Gene, Condition, log2FoldChange, lfcSE, pvalue, padj)
}

degradosome_genes <- c("PNPase", "Suv3", "RNaseZ", "Trm1", "mtTFB2", "Nfs1", "bsf")
plot_df_g <- bind_rows(extract_gene_data_for_barplot(res_G11A_head, "G11A Head", degradosome_genes), extract_gene_data_for_barplot(res_G146C_head, "G146C Head", degradosome_genes), extract_gene_data_for_barplot(res_G11A_thorax, "G11A Thorax", degradosome_genes), extract_gene_data_for_barplot(res_G146C_thorax, "G146C Thorax", degradosome_genes))
plot_df_g$Condition <- factor(plot_df_g$Condition, levels = c("G11A Head", "G146C Head", "G11A Thorax", "G146C Thorax"))
plot_df_g$Gene <- factor(plot_df_g$Gene, levels = degradosome_genes)
plot_df_g <- plot_df_g %>% mutate(Group = ifelse(Gene %in% c("PNPase", "Suv3"), "Mito RNA Degradosome", "tRNA/mRNA Processing & Maturation"), SigLabel = case_when(is.na(padj) ~ "", padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""))

fill_colors_bp <- c("G11A Head" = "#4F46E5", "G146C Head" = "#06B6D4", "G11A Thorax" = "#10B981", "G146C Thorax" = "#F59E0B")

p_degradosome <- ggplot(plot_df_g, aes(x = Gene, y = log2FoldChange, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), position = position_dodge(width = 0.8), width = 0.25, color = "black", linewidth = 0.4) +
  geom_text(aes(y = ifelse(log2FoldChange >= 0, log2FoldChange + lfcSE + 0.15, log2FoldChange - lfcSE - 0.25), label = SigLabel), position = position_dodge(width = 0.8), size = 3.5, fontface = "bold", vjust = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x", labeller = label_wrap_gen(width = 20)) + theme_bw() +
  scale_fill_manual(values = fill_colors_bp) + labs(title = "Dysregulation of Mitochondrial RNA Degradosome", x = NULL, y = "log2FC") + theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "degradosome_processing_barplot.pdf"), p_degradosome, width = 11, height = 7)
print(p_degradosome)

thorax_enriched_genes <- c("Gasz", "Tcs4", "Nsun2", "angel", "mRpL44", "Parp1", "THG")
plot_df_h <- bind_rows(extract_gene_data_for_barplot(res_G11A_head, "G11A Head", thorax_enriched_genes), extract_gene_data_for_barplot(res_G146C_head, "G146C Head", thorax_enriched_genes), extract_gene_data_for_barplot(res_G11A_thorax, "G11A Thorax", thorax_enriched_genes), extract_gene_data_for_barplot(res_G146C_thorax, "G146C Thorax", thorax_enriched_genes))
plot_df_h$Condition <- factor(plot_df_h$Condition, levels = c("G11A Head", "G146C Head", "G11A Thorax", "G146C Thorax"))
plot_df_h$Gene <- factor(plot_df_h$Gene, levels = thorax_enriched_genes)
plot_df_h <- plot_df_h %>% mutate(Group = case_when(Gene %in% c("Gasz", "angel") ~ "RNA Decay & piRNA Processing", Gene %in% c("Tcs4", "Nsun2", "THG") ~ "tRNA Modification & Guanylylation", Gene %in% c("mRpL44") ~ "Mito Translation & Ribosome", TRUE ~ "Chromatin & RNA Regulation"), SigLabel = case_when(is.na(padj) ~ "", padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*", TRUE ~ ""))

p_thorax_enriched <- ggplot(plot_df_h, aes(x = Gene, y = log2FoldChange, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), position = position_dodge(width = 0.8), width = 0.25, color = "black", linewidth = 0.4) +
  geom_text(aes(y = ifelse(log2FoldChange >= 0, log2FoldChange + lfcSE + 0.15, log2FoldChange - lfcSE - 0.25), label = SigLabel), position = position_dodge(width = 0.8), size = 3.5, fontface = "bold", vjust = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x", labeller = label_wrap_gen(width = 20)) + theme_bw() +
  scale_fill_manual(values = fill_colors_bp) + labs(title = "Thorax-Enriched Mitochondrial RNA Processing", x = NULL, y = "log2FC") + theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "degradosome_processing_thorax_enriched_barplot.pdf"), p_thorax_enriched, width = 12, height = 7)
print(p_thorax_enriched)
