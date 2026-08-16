# ------------------------------------------------------------------------------
# 08_exportar_tablas_anexoA.R
# TFM - Análisis de expresión génica y su asociación con la supervivencia 
# en cáncer de mama mediante datos públicos RNA-seq
# Fase 8: Exportación de las tablas completas de resultados a CSV (Anexo A)
# Autora: Silvia Cuesta Cordón
# ------------------------------------------------------------------------------

# 0. Recuperamos los objetos generados en las fases anteriores------------------

top_genes <- readRDS(file.path("resultados", "rds", "top_genes.rds"))
resultados_supervivencia <- readRDS(file.path("resultados", "rds", "resultados_supervivencia.rds"))
enriquecimiento_GO <- readRDS(file.path("resultados", "rds", "enriquecimiento_GO.rds"))
enriquecimiento_KEGG_exploratorio <- readRDS(file.path("resultados", "rds", "enriquecimiento_KEGG_exploratorio.rds"))
modelo_cox <- readRDS(file.path("resultados", "rds", "modelo_cox.rds"))


# 1. Tabla A1: los 50 genes candidatos (expresión diferencial completa)---------

# Se exporta el conjunto completo de estadísticos de DESeq2 para los 50 genes
# candidatos, ordenados de menor a mayor p-valor ajustado (padj), como
# referencia exhaustiva del análisis resumido en el apartado 5.1.

tabla_top50 <- as.data.frame(top_genes)[, c("ensembl_id", "symbol", "baseMean",
                                            "log2FoldChange", "pvalue", "padj")]
tabla_top50 <- tabla_top50[order(tabla_top50$padj), ]
write.csv(tabla_top50,
          file.path("resultados", "tablas", "anexo_tabla1_top50_genes_expresion_diferencial.csv"),
          row.names = FALSE)


# 2. Tabla A2: resultados de supervivencia para los 50 genes candidatos---------

write.csv(resultados_supervivencia,
          file.path("resultados", "tablas", "anexo_tabla2_resultados_supervivencia.csv"),
          row.names = FALSE)


# 3. Tabla A3: coeficientes del modelo de Cox multivariante---------------------

# Se extraen los coeficientes directamente de summary(modelo_cox), en lugar
# de transcribirlos manualmente, para evitar errores de formato o de
# transcripción al pasar los resultados a la memoria.

resumen_cox <- summary(modelo_cox)
tabla_cox <- data.frame(
  gen          = rownames(resumen_cox$coefficients),
  coef         = resumen_cox$coefficients[, "coef"],
  hazard_ratio = resumen_cox$coefficients[, "exp(coef)"],
  se_coef      = resumen_cox$coefficients[, "se(coef)"],
  z            = resumen_cox$coefficients[, "z"],
  p_valor      = resumen_cox$coefficients[, "Pr(>|z|)"],
  ic95_inf     = resumen_cox$conf.int[, "lower .95"],
  ic95_sup     = resumen_cox$conf.int[, "upper .95"]
)
write.csv(tabla_cox,
          file.path("resultados", "tablas", "anexo_tabla3_modelo_cox_multivariante.csv"),
          row.names = FALSE)


# 4. Tabla A4: enriquecimiento funcional GO (procesos biológicos)---------------

tabla_go <- as.data.frame(enriquecimiento_GO)
write.csv(tabla_go,
          file.path("resultados", "tablas", "anexo_tabla4_enriquecimiento_GO.csv"),
          row.names = FALSE)

# 5. Tabla A5: enriquecimiento funcional KEGG (exploratorio, sin filtro)--------

tabla_kegg <- as.data.frame(enriquecimiento_KEGG_exploratorio)
write.csv(tabla_kegg,
          file.path("resultados", "tablas", "anexo_tabla5_enriquecimiento_KEGG.csv"),
          row.names = FALSE)

# Resumen de archivos generados-------------------------------------------------

cat("Archivos CSV generados en resultados/tablas/:\n")
cat("- anexo_tabla1_top50_genes_expresion_diferencial.csv\n")
cat("- anexo_tabla2_resultados_supervivencia.csv\n")
cat("- anexo_tabla3_modelo_cox_multivariante.csv\n")
cat("- anexo_tabla4_enriquecimiento_GO.csv\n")
cat("- anexo_tabla5_enriquecimiento_KEGG.csv\n")
