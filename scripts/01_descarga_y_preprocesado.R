#-------------------------------------------------------------------------------
# 01_descarga_y_preprocesado.R
# TFM - Análisis de expresión génica y su asociación con la supervivencia 
# en cáncer de mama mediante datos públicos RNA-seq
# Fase 1: Descarga de datos (TCGA-BRCA) y preprocesado
# Corresponde al Apartado 3 (Datos) de la memoria
# Autora: Silvia Cuesta Cordón
#-------------------------------------------------------------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)


# 0. Creación de la estructura de carpetas de salida-----------------------------

# Se crean (si no existen ya) las carpetas donde se guardarán los objetos
# intermedios (.rds), las figuras y las tablas generadas por todo el pipeline,
# manteniendo la entrega organizada por tipo de resultado.

dir.create(file.path("resultados", "rds"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("resultados", "figuras"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("resultados", "tablas"), recursive = TRUE, showWarnings = FALSE)


# 1. Consulta y descarga de expresión génica (RNA-seq)--------------------------


# Se descargan conteos crudos (STAR - Counts), ya que DESeq2 requiere raw counts
# para su propio proceso de normalización interno. Se incluyen tanto muestras
# tumorales como de tejido normal adyacente, necesarias para el análisis de
# expresión diferencial del apartado 4.

query_expr <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# Vemos el detalle de las muestras encontradas antes de descargar

resultados <- getResults(query_expr)
nrow(resultados)
table(resultados$sample_type)

# Descarga real de los archivos (tarda varios minutos)

GDCdownload(query_expr)

# Cargamos los datos descargados en un objeto manejable

data_expr <- GDCprepare(query_expr)
dim(data_expr)

# Extraemos la matriz de conteos (raw counts)

counts <- assay(data_expr, "unstranded")
dim(counts)


# 2. Descarga de datos clínicos-------------------------------------------------

clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
dim(clinical)
colnames(clinical)

# De las 99 variables clínicas disponibles, se seleccionan únicamente las
# necesarias para el análisis de supervivencia (estado vital, tiempos) y
# como covariables clínicas de referencia (edad, estadio tumoral).

clinical_reducido <- clinical[, c(
  "submitter_id",
  "vital_status",
  "days_to_death",
  "days_to_last_follow_up",
  "age_at_index",
  "ajcc_pathologic_stage"
)]

summary(clinical_reducido)
table(clinical_reducido$vital_status, useNA = "ifany")
table(clinical_reducido$ajcc_pathologic_stage, useNA = "ifany")


# 3. Simplificación del estadio tumoral-----------------------------------------

clinical_reducido <- clinical_reducido %>%
  mutate(stage_simple = case_when(
    grepl("Stage IV", ajcc_pathologic_stage)  ~ "IV",
    grepl("Stage III", ajcc_pathologic_stage) ~ "III",
    grepl("Stage II", ajcc_pathologic_stage)  ~ "II",
    grepl("Stage I", ajcc_pathologic_stage)   ~ "I",
    grepl("Stage X", ajcc_pathologic_stage)   ~ "X",
    TRUE ~ NA_character_
  ))

table(clinical_reducido$stage_simple, useNA = "ifany")

clinical_reducido$stage_simple[is.na(clinical_reducido$stage_simple)] <- "Sin dato"

clinical_reducido$stage_simple <- factor(
  clinical_reducido$stage_simple,
  levels = c("I", "II", "III", "IV", "X", "Sin dato")
)


# 4. Figura 1: distribución de pacientes por estadio tumoral--------------------

ggplot(clinical_reducido, aes(x = stage_simple)) +
  geom_bar(fill = "#F4C0D1") +
  labs(
    title = "Distribución de pacientes por estadio tumoral",
    x = "Estadio (AJCC simplificado)",
    y = "Número de pacientes"
  ) +
  theme_minimal()

ggsave(file.path("resultados", "figuras", "figura_distribucion_estadio.png"), 
       width = 6, height = 4, dpi = 300)


# 5. Cruce de expresión génica y datos clínicos---------------------------------

sample_info <- as.data.frame(colData(data_expr))
sample_info$patient_id <- substr(sample_info$barcode, 1, 12)

datos_completos <- merge(
  sample_info,
  clinical_reducido,
  by.x = "patient_id",
  by.y = "submitter_id"
)

dim(datos_completos)
table(datos_completos$sample_type)

# La función merge() a veces duplica columnas con sufijos ".x"/".y" si ambas
# tablas comparten nombre de columna. Comprobamos si ha ocurrido y, en caso
# afirmativo, lo resolvemos renombrando las columnas procedentes de la tabla
# clínica y eliminando sus duplicados.


nombres <- colnames(datos_completos)
tiene_sufijos <- any(grepl("vital_status\\.", nombres))
print(tiene_sufijos)

if (tiene_sufijos) {
  datos_completos <- datos_completos %>%
    rename(
      vital_status = vital_status.y,
      days_to_death = days_to_death.y,
      days_to_last_follow_up = days_to_last_follow_up.y,
      age_at_index = age_at_index.y,
      ajcc_pathologic_stage = ajcc_pathologic_stage.y
    ) %>%
    select(-vital_status.x, -days_to_death.x, -days_to_last_follow_up.x,
           -age_at_index.x, -ajcc_pathologic_stage.x)
}


# 6. Filtrado de calidad--------------------------------------------------------

sum(is.na(datos_completos$vital_status))
datos_completos <- datos_completos[!is.na(datos_completos$vital_status), ]
dim(datos_completos)

# Nos quedamos solo con las columnas de counts que sobrevivieron al filtro

counts_filtrado <- counts[, colnames(counts) %in% datos_completos$barcode]


# Se emplea el criterio estándar en RNA-seq (mínimo 10 lecturas en al menos
# el 10% de las muestras) para eliminar genes con expresión residual o nula,
# reduciendo ruido técnico antes del análisis de expresión diferencial.

keep <- rowSums(counts_filtrado >= 10) >= (ncol(counts_filtrado) * 0.1)
counts_filtrado <- counts_filtrado[keep, ]
dim(counts_filtrado)


# 7. Guardado de objetos para las siguientes fases------------------------------

saveRDS(datos_completos, file.path("resultados", "rds", "datos_completos.rds"))
saveRDS(counts_filtrado, file.path("resultados", "rds", "counts_filtrado.rds"))

# Para recuperar en la siguiente sesión, sin repetir la descarga:
# datos_completos <- readRDS(file.path("resultados", "rds", "datos_completos.rds"))
# counts_filtrado <- readRDS(file.path("resultados", "rds", "counts_filtrado.rds"))
