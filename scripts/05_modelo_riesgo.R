# ------------------------------------------------------------------------------
# 05_modelo_riesgo.R
# TFM - Análisis de expresión génica y su asociación con la supervivencia 
# en cáncer de mama mediante datos públicos RNA-seq
# Fase 5: Modelo de riesgo pronóstico combinado (Cox multivariante)
# Corresponde al apartado 4.5 de la memoria (objetivo específico 4)
# Autora: Silvia Cuesta Cordón
# ------------------------------------------------------------------------------

library(survival)
library(survminer)
library(dplyr)


# 0. Recuperamos los objetos de la fase anterior--------------------------------

datos_supervivencia <- readRDS(file.path("resultados", "rds", "datos_supervivencia.rds"))
norm_counts_top <- readRDS(file.path("resultados", "rds", "norm_counts_top.rds"))
resultados_supervivencia <- readRDS(file.path("resultados", "rds", "resultados_supervivencia.rds"))
top_genes <- readRDS(file.path("resultados", "rds", "top_genes.rds"))

# Los 4 genes pronósticos significativos identificados en la fase anterior

genes_pronosticos <- resultados_supervivencia$ensembl_id[resultados_supervivencia$p_valor < 0.05]
genes_pronosticos

surv_obj <- Surv(time = datos_supervivencia$tiempo, event = datos_supervivencia$evento)


# 1. Preparamos una tabla con la expresión de los 4 genes, ya con
#    nombres de columna legibles (símbolo en vez de código Ensembl)-------------


expr_4genes <- as.data.frame(t(norm_counts_top[genes_pronosticos, ]))

# expr_4genes hereda el orden de norm_counts_top, que coincide con el orden
# de filas de datos_supervivencia (ambos proceden del mismo script 04).
# Se comprueba explícitamente antes de combinar ambas tablas:

all(rownames(expr_4genes) == datos_supervivencia$barcode)

# Sustituimos los nombres de columna (Ensembl) por el símbolo de cada gen

simbolos_4genes <- top_genes$symbol[match(genes_pronosticos, top_genes$ensembl_id)]
colnames(expr_4genes) <- simbolos_4genes

# Unimos esta expresión a la tabla de supervivencia (mismo orden de filas que 
# datos_supervivencia)

datos_modelo <- cbind(datos_supervivencia, expr_4genes)


# 2. Modelo de Cox multivariante con los 4 genes a la vez-----------------------

# La fórmula usa directamente los nombres de columna (símbolos de gen)

formula_cox <- as.formula(
  paste("surv_obj ~", paste(simbolos_4genes, collapse = " + "))
)

modelo_cox <- coxph(formula_cox, data = datos_modelo)
summary(modelo_cox)


# 3. Puntuación de riesgo combinada (firma génica)------------------------------

# predict(..., type = "lp") da el "linear predictor": combina la expresión
# de los 4 genes ponderada por su coeficiente en el modelo de Cox.
# Es, en esencia, la "firma génica" de riesgo para cada paciente.

datos_modelo$puntuacion_riesgo <- predict(modelo_cox, type = "lp")

# Se dicotomiza la puntuación de riesgo por su mediana, siguiendo el mismo
# criterio empleado en el análisis univariante (script 04), para mantener
# la comparabilidad entre ambos enfoques.

# Dividimos a las pacientes en dos grupos según la mediana de la puntuación

mediana_riesgo <- median(datos_modelo$puntuacion_riesgo)
datos_modelo$grupo_riesgo <- ifelse(
  datos_modelo$puntuacion_riesgo > mediana_riesgo,
  "Alto riesgo",
  "Bajo riesgo"
)

table(datos_modelo$grupo_riesgo)


# 4. Validación visual: Kaplan-Meier según el grupo de riesgo combinado---------

fit_riesgo <- survfit(surv_obj ~ grupo_riesgo, data = datos_modelo)

km_riesgo <- ggsurvplot(
  fit_riesgo,
  data = datos_modelo,
  pval = TRUE,
  risk.table = TRUE,
  title = "Supervivencia según firma génica de riesgo\n(4 genes combinados)",
  xlab = "Tiempo (días)",
  legend.title = "Grupo de riesgo"
)

km_riesgo
png(file.path("resultados", "figuras", "figura_KM_firma_riesgo.png"),
    width = 7, height = 6, units = "in", res = 300)
print(km_riesgo)
dev.off()


# 5. Guardado de resultados-----------------------------------------------------

saveRDS(modelo_cox, file.path("resultados", "rds", "modelo_cox.rds"))
saveRDS(datos_modelo, file.path("resultados", "rds", "datos_modelo.rds"))

# Para recuperar en la siguiente sesión:
# modelo_cox <- readRDS(file.path("resultados", "rds", "modelo_cox.rds"))
# datos_modelo <- readRDS(file.path("resultados", "rds", "datos_modelo.rds"))