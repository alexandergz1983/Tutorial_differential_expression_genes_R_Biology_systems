######################################################################################################################################
# Biología de sistemas "EAFIT" 2024
# Análisis de expresión diferencial de genes
#
# Este script procesa datos de expresión génica, identifica genes diferencialmente expresados y 
# realiza análisis exploratorios (escalado, normalización y ordenación: PCA y PCoA).
#
# NOTA: DESeq2 normaliza automáticamente mediante size factors y se utiliza la transformación VST para 
# estabilizar la varianza (escalado) previo a la exploración mediante PCA y PCoA.
######################################################################################################################################

#########################################
# 1. Instalación y carga de librerías
#########################################
# (Descomente las líneas de instalación si es la primera vez que usa estos paquetes)
# if (!require('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('DESeq2')
# BiocManager::install('EnhancedVolcano')
# install.packages('pheatmap')
# install.packages('ggplot2')

library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)

#########################################
# 2. Carga de datos: Diseño experimental y recuentos sin procesar
#########################################
# Cargar archivo con el diseño experimental. Se asume que las primeras 18 muestras son tumores primarios (cáncer)
# y las siguientes 19 muestras son normales. Para este análisis se utilizan solo las primeras 37.
design <- read.delim('/Users/usuario/Documents/programas/oncology_diff_express/E-GEOD-50760-experiment-design.tsv')
View(design)

# Cargar el conjunto de datos de recuentos sin procesar
raw.counts <- read.delim('/Users/usuario/Documents/programas/oncology_diff_express/E-GEOD-50760-raw-counts.tsv')
View(raw.counts)

#########################################
# 3. Preparación de datos: Selección de muestras, etiquetas y metadata
#########################################
# Seleccionar únicamente las primeras 37 muestras (tumores primarios y normales)
design <- design[1:37,]

# Crear etiquetas para las muestras: 18 de "cancer" y 19 de "normal"
label1 <- rep("cancer", 18)
label2 <- rep("normal", 19)
labels <- c(label1, label2)

# Extraer la columna de nombres de genes y los recuentos correspondientes a las muestras de interés
dataset <- raw.counts[, c("Gene.Name", design$Run)]
genetable <- raw.counts[, design$Run]

# Crear el dataframe de metadatos para DESeq2
Metadata <- data.frame(id = design$Run, type = labels)

#########################################
# 4. Creación del objeto DESeq2 y ejecución del análisis diferencial
#########################################
# Crear el objeto DESeqDataSet; DESeq2 calculará los size factors (normalización) y estimará la dispersión
dds <- DESeqDataSetFromMatrix(genetable, Metadata, ~ type)
dds <- DESeq(dds)

#########################################
# 5. Análisis Exploratorio: Escalado y Ordenación (PCA y PCoA)
#########################################
# Aplicar transformación de varianza estabilizada (VST)
# Esta transformación "escala" los datos para estabilizar la varianza a lo largo de la gama de expresión
vsd <- vst(dds, blind = FALSE)

## Análisis de PCA
# Extraer datos del PCA (incluye porcentaje de varianza explicada) y graficar usando ggplot2
pcaData <- plotPCA(vsd, intgroup = "type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% varianza")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza")) +
  ggtitle("Análisis PCA de muestras RNA-seq") +
  theme_minimal()
print(pca_plot)

## Análisis de PCoA
# Calcular la matriz de distancias (usando los datos transformados por VST)
vsd_mat <- assay(vsd)       # Extraer la matriz de expresión transformada
d <- dist(t(vsd_mat))       # Calcular distancias entre muestras (transponer: cada fila es una muestra)

# Realizar el análisis de coordenadas principales (PCoA) utilizando cmdscale
pcoa <- cmdscale(d, k = 2, eig = TRUE)

# Preparar el dataframe para graficar, incluyendo las coordenadas PCoA y la información de tipo de muestra
pcoa_data <- data.frame(
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  type = colData(vsd)$type
)

# Calcular porcentaje de varianza explicado por cada eje de la PCoA
eig_sum <- sum(pcoa$eig)
var1 <- round(100 * pcoa$eig[1] / eig_sum, 1)
var2 <- round(100 * pcoa$eig[2] / eig_sum, 1)

pcoa_plot <- ggplot(pcoa_data, aes(PCoA1, PCoA2, color = type)) +
  geom_point(size = 3) +
  xlab(paste0("PCoA1: ", var1, "% varianza")) +
  ylab(paste0("PCoA2: ", var2, "% varianza")) +
  ggtitle("Análisis PCoA de muestras RNA-seq") +
  theme_minimal()
print(pcoa_plot)

#########################################
# 6. Obtención de resultados de expresión diferencial
#########################################
# Se extraen los resultados de DESeq2 comparando "cancer" vs "normal".
# Primero se usa la corrección de Benjamini-Hochberg (BH), luego se utiliza la corrección de Holm (más conservadora).

# Resultados con corrección BH
res <- results(object = dds, contrast = c("type", "cancer", "normal"),
               pAdjustMethod = "BH", alpha = 0.000001)
row.names(res) <- dataset$Gene.Name
summary(res)

# Resultados con corrección Holm (más conservadora)
res <- results(object = dds, contrast = c("type", "cancer", "normal"),
               pAdjustMethod = "holm", alpha = 0.000001)
row.names(res) <- dataset$Gene.Name
summary(res)

#########################################
# 7. Visualización de resultados: Gráficos de Volcán
#########################################
# Crear gráfico de volcán utilizando EnhancedVolcano. Se muestran dos versiones: 
# una usando el valor p original y otra usando el valor p ajustado (padj).

# Volcán usando pvalue
EnhancedVolcano(res,
                lab = dataset$Gene.Name,
                x = "log2FoldChange",
                y = "pvalue",
                pCutoff = 10e-5,
                FCcutoff = 1.333,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = "The results",
                subtitle = "Differential expression analysis",
                caption = "log2fc cutoff=1.333; p value cutoff=10e-5",
                legendPosition = "right",
                legendLabSize = 14,
                col = c("lightblue", "orange", "blue", "red2"),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)

# Volcán usando padj
EnhancedVolcano(res,
                lab = dataset$Gene.Name,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 10e-7,
                FCcutoff = 2.5,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = "The results",
                subtitle = "Differential expression analysis",
                caption = "log2fc cutoff=1.333; p value cutoff=10e-6",
                legendPosition = "right",
                legendLabSize = 14,
                col = c("lightblue", "orange", "blue", "red2"),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)

#########################################
# 8. Generación de tabla final y filtrado de genes significativos
#########################################
# Crear un dataframe final ordenado por log2FoldChange
resord <- as.data.frame(res)
finaltable <- cbind(dataset$Gene.Name, resord)
finaltable <- finaltable[order(finaltable$log2FoldChange), ]
write.table(finaltable, file = "/Users/usuario/Documents/programas/oncology_diff_express/finaltable_ori.csv", sep = ",", col.names = NA)

# Filtrar genes significativos: padj < 10e-7 y log2FoldChange >= 2.5 o <= -2.5
significativos <- resord[which(resord$padj < 10e-7 & (resord$log2FoldChange >= 2.5 | resord$log2FoldChange <= -2.5)), ]
significativos_table <- significativos
significativos_table$Gene.Name <- rownames(significativos)
colnames(significativos_table)[ncol(significativos_table)] <- "Gene.Name"
write.table(significativos_table, file = "/Users/usuario/Documents/programas/oncology_diff_express/significativo.csv", sep = ",", row.names = FALSE, col.names = TRUE)

#########################################
# 9. Heatmap de genes significativos
#########################################
# Cargar la tabla filtrada de genes significativos y generar un heatmap (utilizando log2FoldChange)
significativos_table <- read.csv("/Users/usuario/Documents/programas/oncology_diff_express/significativo.csv")
gene_log2FC <- data.frame(Gene.Name = significativos_table$Gene.Name, log2FoldChange = significativos_table$log2FoldChange)
matriz_log2FC <- as.matrix(gene_log2FC$log2FoldChange)
rownames(matriz_log2FC) <- gene_log2FC$Gene.Name

pheatmap(matriz_log2FC,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         border_color = NA,
         cluster_cols = FALSE,  # Sólo hay una columna
         cluster_rows = TRUE,   # Agrupar genes
         main = "Heatmap de log2FoldChange de Genes Significativos",
         fontsize_row = 8,
         show_colnames = FALSE)

#########################################
# 10. Sección final: Enriquecimiento de genes y herramientas complementarias
#########################################
# Para un análisis más profundo, se recomienda utilizar:
# 1. Enriquecimiento de genes: https://maayanlab.cloud/Enrichr/
# 2. Análisis de redes e interacciones: https://string-db.org/
# 3. Bases de datos adicionales: https://ngdc.cncb.ac.cn/databasecommons/database/id/3061
# 4. Visualización y análisis de redes: Cytoscape

######################################################################################################################################
# Conclusión:
# Este script utiliza DESeq2 para detectar genes diferencialmente expresados,
# aplica normalización (size factors) y escalado (transformación VST) para estabilizar la varianza,
# y realiza análisis de ordenación (PCA y PCoA) para explorar la estructura global de los datos.
# Finalmente, se visualizan los resultados con gráficos de volcán y heatmaps, facilitando la identificación de biomarcadores.
######################################################################################################################################
