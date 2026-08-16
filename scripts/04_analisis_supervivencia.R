# ------------------------------------------------------------------------------
# 04_analisis_supervivencia.R
# TFM - Análisis de expresión génica y su asociación con la supervivencia 
# en cáncer de mama mediante datos públicos RNA-seq
# Fase 4: Kaplan-Meier y selección de genes pronósticos
# Corresponde al apartado 4.4 de la memoria (objetivo específico 3)
# Autora: Silvia Cuesta Cordón
# ------------------------------------------------------------------------------

library(survival)
library(survminer)


# 0. Recuperamos los objetos de las fases anteriores----------------------------

datos_completos <- readRDS(file.path("resultados", "rds", "datos_completos.rds"))
dds <- readRDS(file.path("resultados", "rds", "dds.rds"))
top_genes <- readRDS(file.path("resultados", "rds", "top_genes.rds"))
genes_candidatos <- rownames(top_genes)


# 1. Preparación de la variable de supervivencia--------------------------------

# Solo nos interesan las muestras de TUMOR (el tejido normal no tiene
# supervivencia que analizar)

datos_supervivencia <- datos_completos[datos_completos$sample_type == "Primary Tumor", ]

# Tiempo: days_to_death si falleció, days_to_last_follow_up si sigue viva

datos_supervivencia$tiempo <- ifelse(
  datos_supervivencia$vital_status == "Dead",
  datos_supervivencia$days_to_death,
  datos_supervivencia$days_to_last_follow_up
)

# Evento: 1 = fallecida, 0 = viva (censurada)

datos_supervivencia$evento <- ifelse(datos_supervivencia$vital_status == "Dead", 1, 0)

sum(is.na(datos_supervivencia$tiempo))
datos_supervivencia <- datos_supervivencia[!is.na(datos_supervivencia$tiempo), ]
dim(datos_supervivencia)
table(datos_supervivencia$evento)


# 2. Expresión normalizada de los genes candidatos------------------------------

norm_counts <- counts(dds, normalized = TRUE)
norm_counts_superv <- norm_counts[, datos_supervivencia$barcode]
norm_counts_top <- norm_counts_superv[genes_candidatos, ]
dim(norm_counts_top)  # debe ser 50 x nº de muestras tumorales con supervivencia

# Objeto de supervivencia (tiempo + evento), usado en todos los análisis siguientes
surv_obj <- Surv(time = datos_supervivencia$tiempo, event = datos_supervivencia$evento)


# 3. Kaplan-Meier automático para los 50 genes candidatos-----------------------

resultados_supervivencia <- data.frame(
  ensembl_id = character(),
  symbol = character(),
  p_valor = numeric(),
  stringsAsFactors = FALSE
)

for (gen in genes_candidatos) {
  
  # Se dicotomiza la expresión de cada gen según su mediana (criterio estándar
  # en estudios de supervivencia basados en expresión génica), generando dos
  # grupos de tamaño equilibrado para la comparación mediante Kaplan-Meier.
  
  expresion_gen <- norm_counts_top[gen, ]
  mediana_expr <- median(expresion_gen)
  grupo <- ifelse(expresion_gen > mediana_expr, "Alta", "Baja")
  
  # survdiff() calcula directamente el estadístico del test log-rank sin
  # necesidad de generar el gráfico, lo que agiliza el cribado automático
  # de los 50 genes candidatos.
  
  diff <- survdiff(surv_obj ~ grupo)
  p_valor <- 1 - pchisq(diff$chisq, df = 1)
  
  simbolo <- top_genes$symbol[top_genes$ensembl_id == gen]
  
  resultados_supervivencia <- rbind(
    resultados_supervivencia,
    data.frame(ensembl_id = gen, symbol = simbolo, p_valor = p_valor)
  )
}

resultados_supervivencia <- resultados_supervivencia[order(resultados_supervivencia$p_valor), ]
resultados_supervivencia

sum(resultados_supervivencia$p_valor < 0.05)

saveRDS(resultados_supervivencia, file.path("resultados", "rds", "resultados_supervivencia.rds"))

# Genes pronósticos significativos (p < 0.05)
genes_pronosticos <- resultados_supervivencia$ensembl_id[resultados_supervivencia$p_valor < 0.05]
genes_pronosticos


# 4. Curvas de Kaplan-Meier (figuras) para los genes pronósticos----------------

graficar_km <- function(gen_id, nombre_gen) {
  
  expresion_gen <- norm_counts_top[gen_id, ]
  mediana_expr <- median(expresion_gen)
  grupo <- ifelse(expresion_gen > mediana_expr, "Alta expresión", "Baja expresión")
  
  datos_temp <- datos_supervivencia
  datos_temp$grupo <- grupo
  
  fit <- survfit(surv_obj ~ grupo, data = datos_temp)
  
  plot <- ggsurvplot(
    fit,
    data = datos_temp,
    pval = TRUE,
    risk.table = TRUE,
    title = paste("Supervivencia según expresión de", nombre_gen),
    xlab = "Tiempo (días)",
    legend.title = "Grupo"
  )

  # Se emplea el dispositivo gráfico png() en lugar de ggsave(), ya que este
  # último no renderiza correctamente los objetos combinados (curva + tabla
  # de riesgo) que devuelve ggsurvplot().
  
  nombre_archivo <- file.path("resultados", "figuras", paste0("figura_KM_", nombre_gen, ".png"))
  png(nombre_archivo, width = 7, height = 6, units = "in", res = 300)
  print(plot)
  dev.off()
  
  return(plot)
}

# Genera una curva por cada gen en genes_pronosticos (ajustar nombres si cambian
# los resultados al re-ejecutar el pipeline completo)

km_SORBS1    <- graficar_km("ENSG00000095637.22", "SORBS1")
km_KIF4A     <- graficar_km("ENSG00000090889.12", "KIF4A")
km_UBE2C     <- graficar_km("ENSG00000175063.17", "UBE2C")
km_LINC01614 <- graficar_km("ENSG00000230838.1", "LINC01614")


# 5. Guardado de objetos para la siguiente fase (modelo de riesgo combinado)----

saveRDS(datos_supervivencia, file.path("resultados", "rds", "datos_supervivencia.rds"))
saveRDS(norm_counts_top, file.path("resultados", "rds", "norm_counts_top.rds"))

# Para recuperar en la siguiente sesión:
# datos_supervivencia <- readRDS(file.path("resultados", "rds", "datos_supervivencia.rds"))
# norm_counts_top <- readRDS(file.path("resultados", "rds", "norm_counts_top.rds"))
# genes_pronosticos <- readRDS(file.path("resultados", "rds", "resultados_supervivencia.rds"))
# y filtrar p < 0.05