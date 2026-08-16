# ------------------------------------------------------------------------------
# 03_seleccion_genes_candidatos.R
# TFM - Análisis de expresión génica y su asociación con la supervivencia 
# en cáncer de mama mediante datos públicos RNA-seq
# Fase 3: Selección del top de genes candidatos y anotación de símbolos
# Corresponde al apartado 4.3 de la memoria
# Autora: Silvia Cuesta Cordón
# ------------------------------------------------------------------------------

library(org.Hs.eg.db)
library(AnnotationDbi)


# 0. Recuperamos los resultados de 02_expresion_diferencial.R-------------------

sig_genes <- readRDS(file.path("resultados", "rds", "sig_genes.rds"))


# 1. Selección del top 50 genes por significancia (padj)------------------------

# Se seleccionan los 50 genes con menor padj como criterio pragmático de
# reducción, dado el elevado número de genes que superaron el filtro de
# significancia en la fase anterior (apartado 4.3 de la memoria).

top_genes <- head(sig_genes[order(sig_genes$padj), ], 50)

summary(top_genes$log2FoldChange)
range(top_genes$padj)

genes_candidatos <- rownames(top_genes)
head(genes_candidatos)


# 2. Conversión de identificadores Ensembl a nombres de gen (símbolo)-----------

# Quitamos el sufijo de versión (".10", ".13"...) antes de anotar

genes_candidatos_sin_version <- sub("\\..*", "", genes_candidatos)

# multiVals = "first" indica que, si un identificador Ensembl tuviera más
# de una correspondencia posible en la base de anotación, se conserva
# únicamente la primera, evitando duplicados en la tabla de resultados.

simbolos <- mapIds(
  org.Hs.eg.db,
  keys = genes_candidatos_sin_version,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

top_genes$ensembl_id <- rownames(top_genes)
top_genes$symbol <- simbolos[genes_candidatos_sin_version]

head(top_genes[, c("ensembl_id", "symbol", "log2FoldChange", "padj")], 10)


# 3. Guardado de objetos para la siguiente fase---------------------------------

saveRDS(top_genes, file.path("resultados", "rds", "top_genes.rds"))

# Para recuperar en la siguiente sesión:
# top_genes <- readRDS(file.path("resultados", "rds", "top_genes.rds"))