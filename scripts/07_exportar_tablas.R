# ------------------------------------------------------------------------------
# 07_exportar_tablas.R
# TFM - Análisis de expresión génica y su asociación con la supervivencia 
# en cáncer de mama mediante datos públicos RNA-seq
# Fase 7: Exportación a Excel de las tablas resumidas del apartado 5
# Autora: Silvia Cuesta Cordón
# ------------------------------------------------------------------------------

# Instalamos el paquete necesario (solo la primera vez)
# install.packages("openxlsx")

library(openxlsx)



# 0. Recuperamos los objetos necesarios------------------------------------------

top_genes <- readRDS(file.path("resultados", "rds", "top_genes.rds"))
resultados_supervivencia <- readRDS(file.path("resultados", "rds", "resultados_supervivencia.rds"))
modelo_cox <- readRDS(file.path("resultados", "rds", "modelo_cox.rds"))
enriquecimiento_KEGG_exploratorio <- readRDS(file.path("resultados", "rds", "enriquecimiento_KEGG_exploratorio.rds"))

# Los p-valores muy pequeños se formatean explícitamente como texto en
# notación científica, ya que Excel reinterpreta erróneamente estos valores
# si se introducen o se pegan como número (deformándolos a un exponente
# incorrecto, del tipo "E+14").

formato_cientifico <- function(x, decimales = 2) {
  formatC(x, format = "e", digits = decimales)
}


# 1. Tabla 1: los 4 genes pronósticos (expresión diferencial + supervivencia)---

genes_pronosticos_ids <- c("ENSG00000230838.1", "ENSG00000090889.12",
                            "ENSG00000095637.22", "ENSG00000175063.17")

tabla1 <- as.data.frame(top_genes)
tabla1 <- tabla1[tabla1$ensembl_id %in% genes_pronosticos_ids,
                  c("symbol", "log2FoldChange", "padj")]
tabla1 <- merge(tabla1, resultados_supervivencia[, c("symbol", "p_valor")], by = "symbol")
tabla1 <- tabla1[order(tabla1$padj), ]

tabla1_final <- data.frame(
  Gen = tabla1$symbol,
  `log2FC (tumor vs. normal)` = round(tabla1$log2FoldChange, 2),
  `padj (DESeq2)` = formato_cientifico(tabla1$padj),
  `p-valor (log-rank)` = round(tabla1$p_valor, 5),
  check.names = FALSE
)


# 2. Tabla 2: coeficientes del modelo de Cox multivariante----------------------

resumen_cox <- summary(modelo_cox)
tabla2_final <- data.frame(
  Gen = rownames(resumen_cox$coefficients),
  `HR (exp(coef))` = round(resumen_cox$coefficients[, "exp(coef)"], 4),
  `p-valor` = round(resumen_cox$coefficients[, "Pr(>|z|)"], 4),
  `IC 95% inf` = round(resumen_cox$conf.int[, "lower .95"], 4),
  `IC 95% sup` = round(resumen_cox$conf.int[, "upper .95"], 4),
  check.names = FALSE
)


# 3. Tabla 3: los 15 genes candidatos con menor padj + su p-valor de
#    supervivencia--------------------------------------------------------------

tabla3 <- as.data.frame(top_genes)
tabla3 <- tabla3[order(tabla3$padj), ][1:15, c("symbol", "log2FoldChange", "padj")]
tabla3 <- merge(tabla3, resultados_supervivencia[, c("symbol", "p_valor")],
                 by = "symbol", all.x = TRUE)
tabla3 <- tabla3[order(tabla3$padj), ]

tabla3_final <- data.frame(
  Gen = tabla3$symbol,
  `log2FC` = round(tabla3$log2FoldChange, 2),
  `padj (DESeq2)` = formato_cientifico(tabla3$padj),
  `p-valor supervivencia` = round(tabla3$p_valor, 3),
  check.names = FALSE
)


# 4. Tabla 4: rutas KEGG exploratorias (top 5)----------------------------------

tabla4 <- as.data.frame(enriquecimiento_KEGG_exploratorio)
tabla4 <- tabla4[order(tabla4$pvalue), ][1:5, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")]

tabla4_final <- data.frame(
  ID = tabla4$ID,
  `Ruta KEGG` = tabla4$Description,
  `Genes implicados` = tabla4$GeneRatio,
  `p (sin ajustar)` = round(tabla4$pvalue, 4),
  `p.adjust` = round(tabla4$p.adjust, 3),
  check.names = FALSE
)


# 5. Guardamos las 4 tablas en un único Excel, cada una en su pestaña-----------

wb <- createWorkbook()

addWorksheet(wb, "Tabla1_genes_pronosticos")
writeData(wb, "Tabla1_genes_pronosticos", tabla1_final)

addWorksheet(wb, "Tabla2_modelo_cox")
writeData(wb, "Tabla2_modelo_cox", tabla2_final)

addWorksheet(wb, "Tabla3_top15_genes")
writeData(wb, "Tabla3_top15_genes", tabla3_final)

addWorksheet(wb, "Tabla4_KEGG_exploratorio")
writeData(wb, "Tabla4_KEGG_exploratorio", tabla4_final)

# Se aplica el mismo estilo de cabecera (fondo azul, texto blanco) que el
# resto de tablas de la memoria, para mantener una identidad visual coherente.

estilo_cabecera <- createStyle(fgFill = "#2A78D6", fontColour = "white",
                                textDecoration = "bold", halign = "center")
for (hoja in c("Tabla1_genes_pronosticos", "Tabla2_modelo_cox",
               "Tabla3_top15_genes", "Tabla4_KEGG_exploratorio")) {
  addStyle(wb, hoja, estilo_cabecera, rows = 1, cols = 1:10, gridExpand = TRUE)
  setColWidths(wb, hoja, cols = 1:10, widths = "auto")
}

saveWorkbook(wb, file.path("resultados", "tablas", "tablas_apartado5_TFM.xlsx"), overwrite = TRUE)

cat("Archivo generado: resultados/tablas/tablas_apartado5_TFM.xlsx (4 pestañas)\n")