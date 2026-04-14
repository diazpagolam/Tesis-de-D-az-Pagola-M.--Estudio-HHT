# PROYECTO DE TESIS: ANÁLISIS DE MICROBIOMA
# SCRIPT DE PROCESAMIENTO Y ESTADÍSTICA

# LIBRERÍAS
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("phyloseq", "ggplot2", "readxl", "dplyr", "tibble", "vegan", 
              "patchwork", "gridExtra", "MASS", "DESeq2", "mixOmics", 
              "microbiome", "scales", "Polychrome", "lefser", 
              "SummarizedExperiment", "ggpubr", "apeglm", "ANCOMBC")  # ← saco vegan duplicado

lapply(packages, function(x) {
  if (!require(x, character.only = TRUE)) BiocManager::install(x, force = TRUE)
})

setwd("C:/Users/Usuario/Desktop/tesis/")

# IMPORTACIÓN Y PREPARACIÓN DE DATOS
COUNTS <- read.csv("counts_polished - COUNTS_polished.csv") %>% column_to_rownames("otu")
TAXCLA <- read.csv("taxla_polished.csv") %>% column_to_rownames("otu")
samples_df <- read.csv("Sample_datar - Hoja 1.csv") %>% column_to_rownames("Individuo")

hht <- phyloseq(
  otu_table(as.matrix(COUNTS), taxa_are_rows = TRUE),
  tax_table(as.matrix(TAXCLA)),
  sample_data(samples_df)
)

hht_filtered <- subset_samples(hht, sample_names(hht) != "HHT6")
hht_filtered <- prune_taxa(taxa_sums(hht_filtered) > 0, hht_filtered)

# RARECURVE
otu.rare <- otu_table(hht_filtered)
otu.rare <- as.data.frame(t(otu.rare))

otu.rarecurve <- rarecurve(otu.rare, step = 10000, label = TRUE)

# DIVERSIDAD ALFA
alpha_df <- estimate_richness(hht_filtered, measures = c("Observed", "Chao1", "Shannon"))
alpha_df$Grupo <- sample_data(hht_filtered)$Grupo

alpha_long <- alpha_df %>%
  tidyr::pivot_longer(cols = c("Observed", "Chao1", "Shannon"),
                      names_to = "Metric", values_to = "Value")

ggplot(alpha_long, aes(x = Grupo, y = Value, fill = Grupo)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f", "#fb9a99")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

kruskal.test(Shannon ~ Grupo, data = alpha_df)

# DIVERSIDAD BETA (HHT VS CONTROL)
phy_main <- subset_samples(hht_filtered, Grupo %in% c("HHT", "Control"))
phy_main_rel <- transform_sample_counts(phy_main, function(x) x / sum(x))
dist_bc <- phyloseq::distance(phy_main_rel, method = "bray")
ord_bc <- ordinate(phy_main_rel, method = "PCoA", distance = dist_bc)

meta_df <- as(sample_data(phy_main_rel), "data.frame") %>% rownames_to_column("SampleID")

# PERMANOVA y Betadisper
set.seed(123)
adonis_res <- adonis2(dist_bc ~ Grupo, data = meta_df, permutations = 999)
pval_permanova <- adonis_res$`Pr(>F)`[1]

disper <- betadisper(dist_bc, meta_df$Grupo)
pval_disper <- permutest(disper, permutations = 999)$tab$`Pr(>F)`[1]

scores_df <- as.data.frame(ord_bc$vectors) %>% 
  rownames_to_column("SampleID") %>%
  left_join(meta_df, by = "SampleID")

eig_vals <- ord_bc$values$Eigenvalues
var_exp <- (eig_vals / sum(eig_vals)) * 100

ggplot(scores_df, aes(x = Axis.1, y = Axis.2, color = Grupo)) +
  geom_point(size = 3, alpha = 0.85) +
  stat_ellipse(aes(group = Grupo), type = "t", level = 0.95, linewidth = 1, alpha = 0.6) +
  labs(title = "PCoA - Bray-Curtis",
       subtitle = paste0("PERMANOVA p = ", signif(pval_permanova, 3), " | Betadisper p = ", signif(pval_disper, 3)),
       x = paste0("PCoA 1 (", round(var_exp[1],1), "%)"),
       y = paste0("PCoA 2 (", round(var_exp[2],1), "%)")) +
  scale_color_manual(values = c("HHT" = "#a6cee3", "Control" = "#b2df8a")) +
  theme_bw() + coord_fixed(sqrt(eig_vals[2] / eig_vals[1]))

# FUNCIONES Y GRÁFICOS DE TAXONOMÍA
crear_paleta <- function(niveles) {
  paleta <- createPalette(N = length(niveles), 
                          seedcolors = c("#000000", "#E41A1C", "#377EB8", "#4DAF4A"), M = 5000)
  idx <- c(seq(1, length(paleta), by = 2), seq(2, length(paleta), by = 2))
  paleta <- paleta[idx]; names(paleta) <- niveles
  if("Otros" %in% niveles) paleta["Otros"] <- "#BDBDBD"
  return(paleta)
}

crear_plot_taxonomia <- function(df, rango, titulo, n_top = 10) {
  df <- df %>% filter(!is.na(!!sym(rango)))
  top_taxa <- df %>% group_by(!!sym(rango)) %>% 
    summarise(Total = sum(Abundance), .groups = "drop") %>%
    arrange(desc(Total)) %>% slice_head(n = n_top) %>% pull(!!sym(rango))
  
  df_plot <- df %>%
    mutate(Taxon = ifelse(!!sym(rango) %in% top_taxa, as.character(!!sym(rango)), "Otros")) %>%
    group_by(Sample, Taxon) %>% summarise(Abundance = sum(Abundance), .groups = "drop")
  
  paleta <- crear_paleta(unique(df_plot$Taxon))
  ggplot(df_plot, aes(x = Sample, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = paleta, name = rango) +
    scale_y_continuous(labels = percent) +
    labs(title = titulo, x = "Muestras", y = "Abundancia relativa") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
}

hht_rel <- transform_sample_counts(hht_filtered, function(x) x / sum(x))
df_rel <- psmelt(hht_rel) %>% mutate(Grupo = recode(Grupo, "Positive control" = "PC"))
df_main <- df_rel %>% filter(Grupo %in% c("HHT", "Control"))

print(crear_plot_taxonomia(df_main, "Phylum", "Abundancia por Filo"))
print(crear_plot_taxonomia(df_main, "Genus", "Abundancia por Género"))
print(crear_plot_taxonomia(df_main, "Species", "Abundancia por Especie"))

# ANÁLISIS LEFSE 
phy_genus <- tax_glom(phy_main, taxrank = "Genus")
phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))
sample_data(phy_genus_rel)$Grupo <- factor(sample_data(phy_genus_rel)$Grupo, levels = c("Control", "HHT"))

otu_mat_lefse <- as(otu_table(phy_genus_rel), "matrix")
if (!taxa_are_rows(phy_genus_rel)) otu_mat_lefse <- t(otu_mat_lefse)
tax_mat_lefse <- as.data.frame(tax_table(phy_genus_rel))
rownames(otu_mat_lefse) <- paste0("p_", tax_mat_lefse$Phylum, "|g_", 
                                  ifelse(is.na(tax_mat_lefse$Genus), "Unclassified", tax_mat_lefse$Genus))

se_lefse <- SummarizedExperiment(assays = list(counts = otu_mat_lefse), 
                                 colData = as.data.frame(sample_data(phy_genus_rel)))

set.seed(123)
lefse_result <- lefser(
  se_lefse,
  classCol = "Grupo",
  lda.cutoff = 2,
  wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05,
  bootstrap_n = 100
)

if (nrow(lefse_result) > 0) {
  print(lefserPlot(lefse_result) + ggtitle("Biomarcadores LEfSe") + theme_bw())
}


# DESEQ2 Y DISTRIBUCIÓN DE EFECTOS
phy_deseq <- subset_samples(hht, Grupo %in% c("HHT", "Control"))
phy_deseq <- prune_taxa(taxa_sums(phy_deseq) > 50, phy_deseq)

otu_mat_deseq <- as(otu_table(phy_deseq), "matrix")
if (!taxa_are_rows(phy_deseq)) otu_mat_deseq <- t(otu_mat_deseq)
keep_taxa <- apply(otu_mat_deseq, 1, function(x) sum(x > 0) >= 3)
phy_deseq <- prune_taxa(keep_taxa, phy_deseq)

sample_data(phy_deseq)$condition <- factor(sample_data(phy_deseq)$Grupo, levels = c("Control", "HHT"))

dds <- phyloseq_to_deseq2(phy_deseq, ~ condition)
dds <- DESeq(dds, sfType = "poscounts")
res_shrink <- lfcShrink(dds, coef = "condition_HHT_vs_Control", type = "apeglm")
res_df <- as.data.frame(res_shrink) %>% 
  mutate(padj = ifelse(is.na(padj), 1, padj), log10padj = -log10(padj))

# Volcano Plot
ggplot(res_df, aes(x = log2FoldChange, y = log10padj)) +
  geom_point(color = "grey70", alpha = 0.6) +
  geom_point(data = filter(res_df, padj < 0.05 & abs(log2FoldChange) > 1), color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw()

# Distribución de efectos
ggplot(res_df, aes(x = log2FoldChange)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Distribución de log2 Fold Change", y = "Número de taxones") +
  theme_bw()

# ANCOM-BC2 (PARAMETRIZADO)
phy_ancom <- subset_samples(hht, Grupo %in% c("HHT", "Control"))
sample_data(phy_ancom)$Grupo <- factor(sample_data(phy_ancom)$Grupo, levels = c("Control", "HHT"))

set.seed(123)
ancom_res <- ancombc2(
  data = phy_ancom,
  fix_formula = "Grupo",
  p_adj_method = "BH",
  prv_cut = 0.1,
  lib_cut = 1000,
  group = "Grupo",
  struc_zero = TRUE,
  alpha = 0.05
)

res_ancom <- data.frame(taxa = rownames(ancom_res$res$beta), 
                        lfc = ancom_res$res$beta[, 2], 
                        qval = ancom_res$res$q_val[, 2], 
                        diff = ancom_res$res$diff_robust[, 2])

res_ancom$qval[is.na(res_ancom$qval)] <- 1

# Merge y Top 5(Si lo anterior hubiera dado resultado)
tax_table_df <- as.data.frame(as(tax_table(phy_ancom), "matrix"))
res_completo <- merge(res_ancom, tax_table_df, by.x = "taxa", by.y = "row.names")
print(res_completo %>% arrange(qval) %>% head(5))
