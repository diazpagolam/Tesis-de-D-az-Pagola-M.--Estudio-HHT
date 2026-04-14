# 1. Librerías
library(tidyverse)
library(ggpubr)
library(janitor)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)


# 2. Cargar y limpiar datos
alfa_data <- read_excel("alpha-diversity(1).xlsx") %>% 
  clean_names() %>%
  dplyr::rename(grupo  set_name)

# 3. Transformación a formato largo
alfa_long <- alfa_data %>%
  dplyr::select(grupo, shannon, chao, simpson) %>%
  tidyr::pivot_longer(
    cols  c(shannon, chao, simpson), 
    names_to  "indice", 
    values_to  "valor"
  )

# 4. Ordenar grupos 
alfa_long$grupo <- factor(alfa_long$grupo, levels  c("Control", "HHT"))

# 5. Gráfico
ggplot(alfa_long, aes(x  grupo, y  valor, fill  grupo)) +
  geom_boxplot(alpha  0.7, outlier.shape  NA) +
  geom_jitter(width  0.1, alpha  0.5) +
  facet_wrap(~indice, scales  "free_y") + 
  theme_bw() +
  scale_fill_manual(values  c("Control"  "#3498db", "HHT"  "#e74c3c")) +
  labs(
    x  "Grupo", 
    y  "Valor", 
    title  "Diversidad Alfa mediante EZBioCloud"
  ) +
  stat_compare_means(method  "wilcox.test", label  "p.signif")

#Diversidad beta 

# 1. IMPORTAR EXCEL
setwd("C:/Users/Usuario/Desktop/tesis/")

df_ez <- read_excel("MTPSET_BetaDiversity_PCoA[C,m].xlsx")

# ver nombres reales de columnas
colnames(df_ez)

# 2. LIMPIEZA 

colnames(df_ez) <- c("SampleID", "Grupo", "PC1", "PC2", "PC3")

# Limpiar datos
df_ez <- df_ez %>%
  
  # eliminar fila de varianza explicada
  filter(SampleID ! "% variation explained") %>%
  
  # corregir formato
  mutate(
    Grupo  ifelse(Grupo  "HD", "Control", Grupo),
    
    PC1  as.numeric(gsub(",", ".", PC1)),
    PC2  as.numeric(gsub(",", ".", PC2))
  )

# 3. % VARIANZA EXPLICADA


var_exp <- c(43.005, 19.688)


# 4. PLOT 

p_ez_beta <- ggplot(df_ez, aes(x  PC1, y  PC2, color  Grupo)) +
  
  geom_point(size  3, alpha  0.85) +
  
  stat_ellipse(
    aes(group  Grupo),
    type  "t",
    level  0.95,
    linewidth  1.1,
    alpha  0.6,
    show.legend  FALSE
  ) +
  
  labs(
    title  "PCoA – EzBioCloud",
    subtitle  "Basado en datos procesados en EzBioCloud",
    x  paste0("PC1 (", var_exp[1], "%)"),
    y  paste0("PC2 (", var_exp[2], "%)"),
    color  "Grupo"
  ) +
  
  scale_color_manual(values  c(
    "HHT"  "#a6cee3",
    "Control"  "#b2df8a"
  )) +
  
  theme_bw() +
  theme(
    legend.position  "right",
    panel.grid  element_blank(),
    plot.title  element_text(size  14, face  "bold"),
    plot.subtitle  element_text(size  11)
  )

print(p_ez_beta)

#Grafico de abundancia
#Genero

# 2. CARGAR DATOS
 
df <- read_excel("MTPSET_Composition_Genus.xlsx") %>%
  clean_names()


# 3. SELECCIONAR Y RENOMBRAR

df <- df[, 1:3]
colnames(df) <- c("Taxon", "HHT", "Control")

# 4. LIMPIEZA

df <- df %>%
  mutate(
    Taxon  as.character(Taxon),
    HHT  as.numeric(as.character(HHT)),
    Control  as.numeric(as.character(Control))
  )

df$HHT[is.na(df$HHT)] <- 0
df$Control[is.na(df$Control)] <- 0

df <- as.data.frame(df)

# 5. TOP TAXA
 
df <- df %>%
  mutate(Total  HHT + Control)

top_n <- 10

df_top <- df %>%
  arrange(desc(Total)) %>%
  head(top_n)

# resto como "Otros" 
df_rest <- df %>%
  arrange(desc(Total)) %>%
  tail(nrow(df) - top_n) %>%
  summarise(
    Taxon  "Otros",
    HHT  sum(HHT),
    Control  sum(Control)
  )

df_plot <- bind_rows(df_top, df_rest)

# 6. FORMATO LARGO

df_long <- df_plot %>%
  pivot_longer(
    cols  c(HHT, Control),
    names_to  "Grupo",
    values_to  "Abundancia"
  )

# 7. NORMALIZAR
 
df_long <- df_long %>%
  group_by(Grupo) %>%
  mutate(Abundancia_relativa  Abundancia / sum(Abundancia) * 100) %>%
  ungroup()

# 8. ORDEN TAXA
 
df_long$Taxon <- factor(df_long$Taxon,
                        levels  df_plot$Taxon)

# 9. PALETA 

n_taxa <- length(unique(df_long$Taxon))

# Paleta base fuerte
colores_base <- brewer.pal(12, "Set3")

# Si hay más de 12 taxa, expandimos
colores <- colorRampPalette(colores_base)(n_taxa)

# 
# 10. GRÁFICO 
# 
ggplot(df_long, aes(x  Grupo,
                    y  Abundancia_relativa,
                    fill  Taxon)) +
  geom_bar(stat  "identity",
           position  "stack",
           color  "black",
           linewidth  0.2) +
  
  scale_fill_manual(values  colores) +
  
  scale_y_continuous(expand  c(0,0)) +
  
  labs(
    title  "Gráfico de abundancia relativa a nivel de Género a partir de EzBioCloud",
    x  "",
    y  "Abundancia relativa (%)",
    fill  "Género"
  ) +
  
  theme_classic() +
  theme(
    text  element_text(size  12),
    legend.position  "right",
    legend.key.size  unit(0.5, "cm")
  )


# 1. IMPORTAR DATOS
df <- read_excel("MTPSET_Composition_Species.xlsx")

colnames(df)[1] <- "Taxon"


# 2. LIMPIEZA
df_clean <- df %>%
  mutate(
    Taxon  as.character(Taxon),
    HHT  as.numeric(as.character(HHT)),
    Control  as.numeric(as.character(Control))
  )

df_clean$HHT[is.na(df_clean$HHT)] <- 0
df_clean$Control[is.na(df_clean$Control)] <- 0

# 3. TOP TAXA 
top_n <- 10

# Asegurar data.frame puro 
df_clean <- as.data.frame(df_clean)

# Calcular total
df_clean$total <- df_clean$HHT + df_clean$Control

# Ordenar
df_clean <- df_clean[order(-df_clean$total), ]

# Obtener top taxa
top_taxa <- head(df_clean$Taxon, top_n)

# Reasignar "Otros"
df_top <- df_clean
df_top$Taxon[!df_top$Taxon %in% top_taxa] <- "Otros"

# 4. AGRUPAR

df_grouped <- df_top %>%
  group_by(Taxon) %>%
  summarise(
    HHT  sum(HHT),
    Control  sum(Control),
    .groups  "drop"
  )


# 5. FORMATO LARGO + NORMALIZACIÓN
 
df_long <- df_grouped %>%
  pivot_longer(cols  c(HHT, Control),
               names_to  "Grupo",
               values_to  "Abundancia") %>%
  group_by(Grupo) %>%
  mutate(Abundancia_relativa  Abundancia / sum(Abundancia) * 100) %>%
  ungroup()

# Orden de taxa
df_long$Taxon <- factor(df_long$Taxon,
                        levels  df_grouped$Taxon)


# 6. PALETA
 
n_taxa <- length(unique(df_long$Taxon))
colores <- colorRampPalette(brewer.pal(12, "Set3"))(n_taxa)

# 7. GRÁFICO
 
ggplot(df_long, aes(x  Grupo,
                    y  Abundancia_relativa,
                    fill  Taxon)) +
  geom_bar(stat  "identity",
           color  "black",
           linewidth  0.2) +
  scale_fill_manual(values  colores) +
  scale_y_continuous(expand  c(0,0)) +
  labs(
    title  "Gráfico de abundancia relativa a nivel de Especie a partir de EzBioCloud",
    x  "",
    y  "Abundancia relativa (%)",
    fill  "Especie"
  ) +
  theme_classic() +
  theme(
    text  element_text(size  12),
    legend.position  "right"
  )

# 1. IMPORTAR DATOS

df <- read_excel("MTPSET_Composition_Phylum.xlsx") 

colnames(df)[1] <- "Taxon"

# 2. LIMPIEZA

df_clean <- df %>%
  mutate(
    Taxon  as.character(Taxon),
    HHT  as.numeric(as.character(HHT)),
    Control  as.numeric(as.character(Control))
  )

df_clean$HHT[is.na(df_clean$HHT)] <- 0
df_clean$Control[is.na(df_clean$Control)] <- 0


# 3. TOP TAXA 


top_n <- 10

# Asegurar data.frame puro 
df_clean <- as.data.frame(df_clean)

# Calcular total
df_clean$total <- df_clean$HHT + df_clean$Control

# Ordenar
df_clean <- df_clean[order(-df_clean$total), ]

# Obtener top taxa
top_taxa <- head(df_clean$Taxon, top_n)

# Reasignar "Otros"
df_top <- df_clean
df_top$Taxon[!df_top$Taxon %in% top_taxa] <- "Otros"
# 
# 4. AGRUPAR
# 
df_grouped <- df_top %>%
  group_by(Taxon) %>%
  summarise(
    HHT  sum(HHT),
    Control  sum(Control),
    .groups  "drop"
  )

# 
# 5. FORMATO LARGO + NORMALIZACIÓN
# 
df_long <- df_grouped %>%
  pivot_longer(cols  c(HHT, Control),
               names_to  "Grupo",
               values_to  "Abundancia") %>%
  group_by(Grupo) %>%
  mutate(Abundancia_relativa  Abundancia / sum(Abundancia) * 100) %>%
  ungroup()

# Orden de taxa
df_long$Taxon <- factor(df_long$Taxon,
                        levels  df_grouped$Taxon)

# 
# 6. PALETA
# 
n_taxa <- length(unique(df_long$Taxon))
colores <- colorRampPalette(brewer.pal(12, "Set3"))(n_taxa)

# 
# 7. GRÁFICO
# 
ggplot(df_long, aes(x  Grupo,
                    y  Abundancia_relativa,
                    fill  Taxon)) +
  geom_bar(stat  "identity",
           color  "black",
           linewidth  0.1) +
  scale_fill_manual(values  colores) +
  scale_y_continuous(expand  c(0,0)) +
  labs(
    title  "Gráfico de abundancia relativa a nivel de Filo a partir de EzBioCloud",
    x  "",
    y  "Abundancia relativa (%)",
    fill  "Filo"
  ) +
  theme_classic() +
  theme(
    text  element_text(size  12),
    legend.position  "right"
  )


# 
# 1. LIBRERÍAS
# 
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(janitor)

# 
# 2. LEER DATOS
# 
df <- read_excel("Taxonomic_biomarker_discoverykruskal(1).xlsx")

# 
# 3. LIMPIAR NOMBRES
# 
df <- df %>% clean_names()

# 
# 4. FORZAR NOMBRES CORRECTOS (CLAVE)
# 
colnames(df) <- c("taxon", "rank", "taxonomy", "p_value", "fdr", "control", "hht")

# 
# 5. CONVERTIR A NUMÉRICO
# 
df <- df %>%
  mutate(
    p_value  as.numeric(gsub(",", ".", p_value)),
    fdr      as.numeric(gsub(",", ".", fdr)),
    control  as.numeric(gsub(",", ".", control)),
    hht      as.numeric(gsub(",", ".", hht))
  )

# 
# 6. CREAR MÉTRICAS
# 
df <- df %>%
  mutate(
    diff  hht - control,
    abs_diff  abs(diff),
    mean_abundance  (hht + control) / 2,
    lda_like  abs_diff * log10(mean_abundance + 1)
  )

# 
# 7. FILTRADO
# 
df_filtrado <- df %>%
  dplyr::filter(
    p_value < 0.05,
    mean_abundance > 0.1
  ) %>%
  arrange(desc(abs_diff))

# 
# 8. GUARDAR TABLA
# 
write.csv(df_filtrado, "taxones_filtrados.csv", row.names  FALSE)

# 
# 9. TOP TAXONES
# 
top_taxa <- df_filtrado %>%
  slice_max(abs_diff, n  8)

# 
# 10. PREPARAR DATOS PARA PLOT
# 
df_plot <- top_taxa %>%
  dplyr::select(taxon, control, hht) %>%
  pivot_longer(
    cols  c(control, hht),
    names_to  "grupo",
    values_to  "abundancia"
  )

# 11. BARPLOT FINAL

ggplot(df_plot, aes(x  reorder(taxon, abundancia),
                    y  abundancia,
                    fill  grupo)) +
  geom_bar(stat  "identity", position  "dodge") +
  coord_flip() +
  labs(
    title  "Taxones con mayor diferencia entre HHT y controles",
    x  "",
    y  "Abundancia relativa"
  ) +
  theme_minimal()

