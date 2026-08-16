# ------------------------------------------------------------------------------
# 06_interpretacion_biologica.R
# TFM - Análisis de expresión génica y su asociación con la supervivencia 
# en cáncer de mama mediante datos públicos RNA-seq
# Fase 6: Enriquecimiento funcional (GO / KEGG) de los genes candidatos
# Corresponde al apartado 5 de la memoria (objetivo específico 5)
# Autora: Silvia Cuesta Cordón
# ------------------------------------------------------------------------------

# Instalamos los paquetes necesarios (solo la primera vez)
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), 
# update = FALSE, ask = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)


# 0. Recuperamos los genes candidatos (top 50) de la fase 3---------------------

top_genes <- readRDS(file.path("resultados", "rds", "top_genes.rds"))

# Quitamos el sufijo de versión de los identificadores Ensembl

genes_sin_version <- sub("\\..*", "", top_genes$ensembl_id)


# 1. Conversión a identificadores ENTREZ (formato que requieren
#    los análisis de enriquecimiento de clusterProfiler)------------------------

genes_entrez <- bitr(
  genes_sin_version,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

nrow(genes_entrez)  # algunos genes pueden no tener ENTREZID asociado, es normal


# 2. Enriquecimiento funcional: Gene Ontology (procesos biológicos)-------------

# Se emplea la categoría BP (Biological Process) por ser la más relevante
# para interpretar mecanismos de progresión tumoral. Los umbrales de
# pvalueCutoff y qvalueCutoff siguen los valores por defecto recomendados
# por clusterProfiler para este tipo de análisis.

enriquecimiento_GO <- enrichGO(
  gene          = genes_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",              # BP = Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE                # convierte de nuevo ENTREZ a símbolo en el resultado
)

# Vemos los términos GO más significativos

head(as.data.frame(enriquecimiento_GO), 10)

# Gráfico de barras con los términos más relevantes

barplot_go <- barplot(enriquecimiento_GO, showCategory = 15) +
  labs(title = "Procesos biológicos enriquecidos (Gene Ontology)")

barplot_go
ggsave(file.path("resultados", "figuras", "figura_enriquecimiento_GO.png"), 
       plot = barplot_go, width = 8, height = 6, dpi = 300)


# 3. Enriquecimiento funcional: rutas KEGG--------------------------------------

enriquecimiento_KEGG <- enrichKEGG(
  gene          = genes_entrez$ENTREZID,
  organism      = "hsa",              # hsa = Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

head(as.data.frame(enriquecimiento_KEGG), 10)

barplot_kegg <- barplot(enriquecimiento_KEGG, showCategory = 15) +
  labs(title = "Rutas biológicas enriquecidas (KEGG)")

barplot_kegg
ggsave(file.path("resultados", "figuras", "figura_enriquecimiento_KEGG.png"), 
       plot = barplot_kegg, width = 8, height = 6, dpi = 300)

# El análisis de la sección 3 no devolvió ninguna ruta significativa tras la
# corrección por comparaciones múltiples (ver apartado 5.4 de la memoria).
# Se repite el análisis sin filtro de significancia (pvalueCutoff = 1) con
# fines descriptivos, para identificar qué rutas quedaron mejor posicionadas.

enriquecimiento_KEGG_exploratorio <- enrichKEGG(
  gene          = genes_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1
)

head(as.data.frame(enriquecimiento_KEGG_exploratorio), 10)


# 4. Guardado de resultados-----------------------------------------------------

saveRDS(enriquecimiento_GO, file.path("resultados", "rds", "enriquecimiento_GO.rds"))
saveRDS(enriquecimiento_KEGG, file.path("resultados", "rds", "enriquecimiento_KEGG.rds"))
saveRDS(enriquecimiento_KEGG_exploratorio, file.path("resultados", "rds", "enriquecimiento_KEGG_exploratorio.rds"))

# Para recuperar en la siguiente sesión:
# enriquecimiento_GO <- readRDS(file.path("resultados", "rds", "enriquecimiento_GO.rds"))
# enriquecimiento_KEGG <- readRDS(file.path("resultados", "rds", "enriquecimiento_KEGG.rds"))
# enriquecimiento_KEGG_exploratorio <- readRDS(file.path("resultados", "rds", "enriquecimiento_KEGG_exploratorio.rds"))