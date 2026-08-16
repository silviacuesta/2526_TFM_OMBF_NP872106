# ------------------------------------------------------------------------------
# 02_expresion_diferencial.R
# TFM - Análisis de expresión génica y su asociación con la supervivencia 
# en cáncer de mama mediante datos públicos RNA-seq
# Fase 2: Expresión diferencial (DESeq2), PCA y volcano plot
# Corresponde a los apartados 3.4 (PCA) y 4.1-4.2 (DESeq2) de la memoria
# Autora: Silvia Cuesta Cordón
# ------------------------------------------------------------------------------

library(DESeq2)
library(ggplot2)


# 0. Recuperamos los datos generados en 01_descarga_y_preprocesado.R------------

datos_completos <- readRDS(file.path("resultados", "rds", "datos_completos.rds"))
counts_filtrado <- readRDS(file.path("resultados", "rds", "counts_filtrado.rds"))


# 1. Preparación del objeto DESeq2----------------------------------------------

# Nos aseguramos de que las columnas de counts_filtrado y las filas de
# datos_completos estén en el mismo orden (imprescindible para DESeq2)

datos_completos <- datos_completos[match(colnames(counts_filtrado), datos_completos$barcode), ]
all(colnames(counts_filtrado) == datos_completos$barcode)  # debe dar TRUE

# Se emplea sample_type como única variable de diseño, ya que el objetivo de
# esta fase es comparar la expresión entre tejido tumoral y tejido normal
# (variable de interés principal del apartado 4.2).

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtrado,
  colData = datos_completos,
  design = ~ sample_type
)


# 2. PCA de control de calidad (apartado 3.4 de la memoria)---------------------

# Se utiliza blind = TRUE porque esta transformación se emplea únicamente
# con fines exploratorios y de control de calidad (PCA), sin tener en cuenta
# todavía el diseño experimental definido para el análisis diferencial.

vsd <- vst(dds, blind = TRUE)

pca_plot <- plotPCA(vsd, intgroup = "sample_type") + theme_minimal()
pca_plot
ggsave(file.path("resultados", "figuras", "figura_pca_muestras.png"),
       plot = pca_plot, width = 6, height = 4.5, dpi = 300)


# 3. Expresión diferencial: tumor vs. tejido normal-----------------------------

dds <- DESeq(dds)

# Se especifica el contraste en este orden (tumor frente a normal) para que
# un log2FoldChange positivo se interprete como sobreexpresión en tejido
# tumoral, facilitando la lectura de los resultados.

res <- results(dds, contrast = c("sample_type", "Primary Tumor", "Solid Tissue Normal"))
res <- res[order(res$padj), ]

summary(res)
head(res)

# Filtro estándar: significativo (padj < 0.05) Y cambio relevante (|log2FC| > 1)

sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
nrow(sig_genes)


# 4. Volcano plot (figura para la memoria)--------------------------------------

res_df <- as.data.frame(res)

res_df$categoria <- "No significativo"
res_df$categoria[res_df$padj < 0.05 & res_df$log2FoldChange > 1]  <- "Sobreexpresado"
res_df$categoria[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Infraexpresado"
res_df <- res_df[!is.na(res_df$padj), ]

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = categoria)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c(
    "Sobreexpresado" = "#E63946",
    "Infraexpresado" = "#457B9D",
    "No significativo" = "grey80"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    title = "Expresión diferencial: tumor vs. tejido normal",
    x = "log2(Fold Change)",
    y = "-log10(p-valor ajustado)",
    color = "Categoría"
  ) +
  theme_minimal()

volcano
ggsave(file.path("resultados", "figuras", "figura_volcano_plot.png"), 
       plot = volcano, width = 7, height = 5, dpi = 300)


# 5. Guardado de objetos para las siguientes fases-------------------------------

saveRDS(dds, file.path("resultados", "rds", "dds.rds"))
saveRDS(res, file.path("resultados", "rds", "res_deseq2.rds"))
saveRDS(sig_genes, file.path("resultados", "rds", "sig_genes.rds"))

# Para recuperar en la siguiente sesión, sin repetir DESeq2 (varios minutos):
# dds <- readRDS(file.path("resultados", "rds", "dds.rds"))
# res <- readRDS(file.path("resultados", "rds", "res_deseq2.rds"))
# sig_genes <- readRDS(file.path("resultados", "rds", "sig_genes.rds"))