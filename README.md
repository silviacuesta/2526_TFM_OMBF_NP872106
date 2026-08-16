# README — Anexo de código: Análisis de expresión génica y supervivencia en cáncer de mama (TFM)

Este documento describe el pipeline de análisis empleado en el trabajo, compuesto por 8 scripts de R que deben ejecutarse en el orden indicado desde la carpeta `scripts/`. Cada script crea automáticamente (si no existen ya) las carpetas de salida necesarias dentro de `resultados/` y parte de los resultados guardados por el script anterior, por lo que pueden ejecutarse de forma independiente sin necesidad de repetir todo el proceso desde el principio, siempre que los archivos `.rds` correspondientes ya existan.

## Estructura de carpetas

Los propios scripts crean y organizan automáticamente esta estructura de salida (no es necesario crearla a mano):

```
Codigo_TFM/
├── .gitignore
├── 00_README.md
├── session_info.txt
├── scripts/
│   ├── 01_descarga_y_preprocesado.R
│   ├── 02_expresion_diferencial.R
│   ├── 03_seleccion_genes_candidatos.R
│   ├── 04_analisis_supervivencia.R
│   ├── 05_modelo_riesgo.R
│   ├── 06_interpretacion_biologica.R
│   ├── 07_exportar_tablas.R
│   └── 08_exportar_tablas_anexoB.R
└── resultados/
    ├── rds/        (objetos intermedios .rds)
    ├── figuras/    (todas las figuras .png)
    └── tablas/     (tablas .xlsx y .csv)
```

Los scripts deben ejecutarse con el directorio de trabajo situado en la raíz de `Codigo_TFM/` (no dentro de `scripts/`), para que las rutas relativas `resultados/rds`, `resultados/figuras` y `resultados/tablas` se resuelvan correctamente.

## Requisitos previos

- R versión 4.4.3 (ver detalle completo de paquetes y versiones en `session_info.txt`).
- Conexión a internet (necesaria únicamente para el script 01, que descarga los datos desde el Genomic Data Commons, y para el script 06, que consulta las bases de datos KEGG).
- Espacio en disco: aproximadamente 3-5 GB libres para los datos descargados en el script 01.
- Los paquetes necesarios se cargan al principio de cada script mediante `library()`. Las instrucciones de instalación (`BiocManager::install(...)` / `install.packages(...)`) aparecen comentadas al inicio de los scripts que introducen un paquete nuevo por primera vez.

## Origen de los datos

Los datos de partida son públicos y no se incluyen como archivo en este anexo, dado su tamaño (varios GB). Proceden del proyecto **TCGA-BRCA** (*Breast Invasive Carcinoma*), disponible en el **Genomic Data Commons (GDC)** del National Cancer Institute: https://portal.gdc.cancer.gov/projects/TCGA-BRCA. El script `01_descarga_y_preprocesado.R` descarga automáticamente estos datos mediante el paquete `TCGAbiolinks`, en una carpeta `GDCdata/` que no se incluye en el repositorio (ver `.gitignore`).

## Descripción de los scripts

| Script | Descripción | Entradas (`resultados/rds/`) | Salidas | Apartado de la memoria |
|---|---|---|---|---|
| `01_descarga_y_preprocesado.R` | Descarga de expresión génica (RNA-seq) y datos clínicos de TCGA-BRCA; cruce y filtrado de calidad. Crea la estructura de carpetas `resultados/`. | — (descarga desde GDC) | `datos_completos.rds`, `counts_filtrado.rds` | Apartado 3 (Datos) |
| `02_expresion_diferencial.R` | Construcción del objeto DESeq2, PCA de control de calidad, análisis de expresión diferencial (tumor vs. normal) y volcano plot. | `datos_completos.rds`, `counts_filtrado.rds` | `dds.rds`, `res_deseq2.rds`, `sig_genes.rds` | Apartados 3.4 y 4.1-4.2 |
| `03_seleccion_genes_candidatos.R` | Selección de los 50 genes con menor p-valor ajustado y anotación de su símbolo génico. | `sig_genes.rds` | `top_genes.rds` | Apartado 4.3 |
| `04_analisis_supervivencia.R` | Preparación de la variable de supervivencia, curvas de Kaplan-Meier univariantes para los 50 genes candidatos y selección de los genes pronósticos significativos. | `datos_completos.rds`, `dds.rds`, `top_genes.rds` | `resultados_supervivencia.rds`, `datos_supervivencia.rds`, `norm_counts_top.rds` | Apartado 4.4 |
| `05_modelo_riesgo.R` | Modelo de Cox multivariante combinando los genes pronósticos y construcción de la puntuación de riesgo (firma génica). | `datos_supervivencia.rds`, `norm_counts_top.rds`, `resultados_supervivencia.rds`, `top_genes.rds` | `modelo_cox.rds`, `datos_modelo.rds` | Apartado 4.5 |
| `06_interpretacion_biologica.R` | Análisis de enriquecimiento funcional (Gene Ontology y KEGG) sobre los genes candidatos. | `top_genes.rds` | `enriquecimiento_GO.rds`, `enriquecimiento_KEGG.rds`, `enriquecimiento_KEGG_exploratorio.rds` | Apartado 5.4 |
| `07_exportar_tablas.R` | Exportación a Excel de las 4 tablas resumidas empleadas en el cuerpo del apartado 5 (genes pronósticos, modelo de Cox, top 15 genes, rutas KEGG exploratorias). | `top_genes.rds`, `resultados_supervivencia.rds`, `modelo_cox.rds`, `enriquecimiento_KEGG_exploratorio.rds` | `resultados/tablas/tablas_apartado5_TFM.xlsx` | Apartado 5 |
| `08_exportar_tablas_anexoB.R` | Exportación a CSV de las 5 tablas completas de resultados (Anexo B). | `top_genes.rds`, `resultados_supervivencia.rds`, `enriquecimiento_GO.rds`, `enriquecimiento_KEGG_exploratorio.rds`, `modelo_cox.rds` | 5 archivos `.csv` en `resultados/tablas/` (ver más abajo) | Anexo B |

## Figuras generadas (en `resultados/figuras/`)

| Archivo | Script que lo genera | Figura en la memoria |
|---|---|---|
| `figura_distribucion_estadio.png` | 01 | Figura 1 (Apartado 3.3) |
| `figura_pca_muestras.png` | 02 | Figura 2 (Apartado 3.3) |
| `figura_volcano_plot.png` | 02 | Figura 4 (Apartado 5.1) |
| `figura_KM_SORBS1.png` | 04 | Figura 5 (Apartado 5.2) |
| `figura_KM_KIF4A.png` | 04 | Figura 6 (Apartado 5.2) |
| `figura_KM_UBE2C.png` | 04 | Figura 7 (Apartado 5.2) |
| `figura_KM_LINC01614.png` | 04 | Figura 8 (Apartado 5.2) |
| `figura_KM_firma_riesgo.png` | 05 | Figura 9 (Apartado 5.3) |
| `figura_enriquecimiento_GO.png` | 06 | Figura 10 (Apartado 5.4) |
| `figura_enriquecimiento_KEGG.png` | 06 | (no significativo, no incluida en el cuerpo de la memoria) |

## Tablas del cuerpo del apartado 5 (script 07, en `resultados/tablas/tablas_apartado5_TFM.xlsx`)

| Pestaña del Excel | Contenido | Tabla en la memoria |
|---|---|---|
| `Tabla3_top15_genes` | Los 15 genes con menor padj en expresión diferencial, junto a su p-valor de supervivencia. | Tabla 2 (Apartado 5.1) |
| `Tabla1_genes_pronosticos` | Los 4 genes con asociación significativa a supervivencia, con sus estadísticos de expresión diferencial. | Tabla 3 (Apartado 5.2) |
| `Tabla2_modelo_cox` | Coeficientes del modelo de Cox multivariante (HR, IC 95%, p-valor). | Tabla 4 (Apartado 5.3) |
| `Tabla4_KEGG_exploratorio` | Las 5 rutas KEGG mejor posicionadas en el análisis exploratorio. | Tabla 6 (Apartado 5.4) |

## Tablas exportadas (Anexo B, script 08, en `resultados/tablas/`)

| Archivo CSV | Contenido |
|---|---|
| `anexo_tabla1_top50_genes_expresion_diferencial.csv` | Los 50 genes candidatos con sus estadísticos de expresión diferencial (baseMean, log2FoldChange, pvalue, padj). |
| `anexo_tabla2_resultados_supervivencia.csv` | P-valor del test de rangos logarítmicos para cada uno de los 50 genes candidatos. |
| `anexo_tabla3_modelo_cox_multivariante.csv` | Coeficientes, hazard ratio e intervalo de confianza del modelo de Cox multivariante. |
| `anexo_tabla4_enriquecimiento_GO.csv` | Resultado completo del enriquecimiento funcional de Gene Ontology (proceso biológico). |
| `anexo_tabla5_enriquecimiento_KEGG.csv` | Resultado completo del enriquecimiento funcional KEGG (exploratorio, sin filtro de significancia). |

## Archivos no incluidos en el repositorio (ver `.gitignore`)


- `GDCdata/`: datos brutos descargados por el script 01 (varios GB, de acceso público, no es necesario versionarlos).
- `resultados/rds/dds.rds`: objeto `DESeqDataSet` con la matriz de expresión completa, generado por el script 02. Se excluye por superar el límite de 100 MB por archivo de GitHub (~1 GB). Se regenera automáticamente y en pocos minutos al ejecutar `02_expresion_diferencial.R` a partir de `datos_completos.rds` y `counts_filtrado.rds`.
- `resultados/rds/counts_filtrado.rds`: matriz de conteos filtrada (60.660 → 26.133 genes), generada por el script 01. Se excluye por su tamaño (~59 MB), al superar el límite de 25 MB por archivo de la interfaz web de GitHub empleada para la subida del repositorio. Se regenera automáticamente al ejecutar `01_descarga_y_preprocesado.R`.
- `.Rproj.user/`, `.Rhistory`, `.RData`, `.RDataTmp`: archivos de sesión de RStudio, sin valor para la reproducibilidad del análisis.

## Reproducibilidad

La versión exacta de R y de todos los paquetes utilizados se encuentra documentada en `session_info.txt`, generado mediante `sessionInfo()` al finalizar el análisis. Para reproducir el pipeline completo desde cero, basta con ejecutar los 8 scripts en orden numérico, con el directorio de trabajo situado en la raíz del proyecto (no dentro de `scripts/`), partiendo de una carpeta `resultados/` vacía o inexistente (los scripts la crean automáticamente).
