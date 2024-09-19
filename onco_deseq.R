# Biologia de sistemas "EAFIT" 2024
# Analisis de expresión diferencial de genes

#install DESeq2 using Bioconductor
#if (!require('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('DESeq2')

#Install EnhancedVolcano using Bioconductor
#if (!require('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')

######################################################################################################################################

#Load the libraries
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)

######################################################################################################################################

# Cargar el archivo que contiene el diseño del experimento
#'path'= ubicación del archivo
#'# Como se puede observar, las primeras 18 muestras son tumores primarios (cáncer colorrectal) y las siguientes 19 son normales
# y el resto son muestras de metástasis hepática. Para este tutorial, utilizaremos solo las primeras 37, es decir, 
# muestras de tumores primarios y normales. Recuerde esta información, será necesaria para la manipulación de los datos.

design <- read.delim('/Users/usuario/Documents/programas/oncology_diff_express/E-GEOD-50760-experiment-design.tsv')
View(design)

######################################################################################################################################

# Ahora veamos el conjunto de datos de recuentos sin procesar y asegurémonos de que los datos de recuento estén allí.

# cargar "raw counts data"
#'path'= ubicación del archivo

raw.counts <- read.delim('/Users/usuario/Documents/programas/oncology_diff_express/E-GEOD-50760-raw-counts.tsv')
View(raw.counts)

# Ahora que todos los datos están confirmados, validados y examinados, se puede comenzar la parte de procesamiento de datos.

######################################################################################################################################

# Este archivo es importante porque contiene columnas que relacionan las muestras con los fenotipos
design=design[1:37,]

# Crear etiquetas para muestras normales y de cáncer primario
label1=rep(c('cancer'),18)
label2=rep(c('normal'),19)
labels=c(label1,label2)


# Cargar la tabla genética
# Seleccione los datos de las variables de interés y añada los nombres de los genes
dataset=raw.counts[,c('Gene.Name',design$Run)]
genetable=raw.counts[,design$Run]


# Crear el marco de datos de metadatos (Create the metadata dataframe)
Metadata = data.frame(id=design$Run,type=labels)

# Como se puede ver, creé manualmente las etiquetas utilizando la información del objeto de diseño y agregué las etiquetas al objeto 
# de metadatos.

# El objeto raw.counts, contiene los recuentos sin procesar de la expresión de ARN y los identificadores de genes, se utiliza 
# para crear una nueva tabla de genes que contiene solo las 37 muestras de interés (esto se logra utilizando design$Run, 
# que contiene los nombres de cadena para las muestras de interés). Los objetos de metadatos de la tabla de genes proporcionarán 
# los recuentos sin procesar y las etiquetas y la información de los identificadores de SRR al algoritmo DESeq durante el análisis.

######################################################################################################################################

# El código DESeq es bastante simple, solo requiere dos pasos: crear el DESeqDataSet y usar la función DESeq() en él; 
# sin embargo, DESeq es un algoritmo altamente complejo: realiza muchos cálculos en grandes conjuntos de datos (en este caso, más de 
# 50 000 genes).

# Dicho esto, puede llevar un tiempo hasta que R finalice todos los análisis y debería recibir mensajes en la consola de que todos 
# estos análisis se están ejecutando.

# Aquí hay algunos cálculos que está ejecutando DESeq…

# 1. Cálculo del factor de tamaño
# 2. Cálculo de dispersión por gen
# 3. Ajuste de dispersión
# 4. Encogiéndose hacia la estimación de la curva
# y otros detalles.

# Realizar análisis DESeq
dds = DESeqDataSetFromMatrix(genetable,Metadata,~type)
dds <- DESeq(dds)

######################################################################################################################################

# Después de eso, se utiliza la función results() para derivar y resumir los resultados de DESeq. Observe cómo establecí el umbral 
# alfa en p<0,000001. Analicemos los resultados...

# Crear objeto de resultados DESeq utilizando la corrección Benjamini-Hochberg
res=results(object = dds, contrast = c('type','cancer','normal'),
            pAdjustMethod = 'BH',alpha = 0.000001)
row.names(res)=dataset$Gene.Name
summary(res)

# De un total de 40.166 genes distintos de cero, hay 121 (el 0,3 %) de genes regulados positivamente y 30 (el 0,075 %) de genes regulados negativamente. 
# Por lo tanto, serían unos 151 genes de interés de más de 50.000. La primera conclusión es que hay muchos más genes regulados positivamente que podrían 
# ser factores causales potenciales de este tipo de cáncer, pero los 30 genes regulados negativamente tampoco son algo que se pueda ignorar.

# Una cosa es muy importante y es el procedimiento de corrección de pruebas múltiples. Más de 40 000 genes distintos de cero significa que el 
# análisis incluirá más de 40 000 pruebas de hipótesis. FWER o tasas de error por familia ocurren con tasas altas en estas situaciones, 
# encontrar resultados positivos al azar es casi seguro si se implementa el procedimiento de corrección múltiple. 
# Es por eso que DESeq siempre tiene algunos procedimientos de corrección múltiple.

######################################################################################################################################

# Agregué pAdjustMethod = 'BH', lo que significa usar la corrección Benjamini-Hochberg. Lo hice a propósito para que puedas ver cómo implementar 
# los procedimientos de prueba múltiples, aunque se implementen de manera predeterminada. Pero el procedimiento de prueba múltiple es muy 
# importante en este caso y quiero implementar un método más sólido. Por lo tanto, usaré el ajuste de prueba múltiple de Holm, que será más sólido 
# en FWER.

# Crear objeto de resultados DESeq utilizando la corrección Holmr

# results: Extrae los resultados del análisis de expresión diferencial a partir del objeto dds (posiblemente generado por DESeq2), 
# utilizando un contraste entre dos condiciones: cáncer vs normal.

# contrast = c('type', 'cancer', 'normal'): Especifica que se están comparando las muestras de cáncer contra las normales basándose en la columna type.
# pAdjustMethod = 'holm': Ajusta los valores p utilizando el método de Holm para corregir por pruebas múltiples, lo cual es una opción conservadora.
# alpha = 0.000001: Define un umbral de significancia muy estricto (1e-6) para considerar los genes como significativamente expresados diferencialmente.
# row.names(res) = dataset$Gene.Name: Asigna los nombres de los genes como identificadores de fila en el resultado (res).
# summary(res): Muestra un resumen de los resultados obtenidos, incluyendo el número de genes que cumplen los criterios de significancia.

res=results(object = dds, contrast = c('type','cancer','normal'),
            pAdjustMethod = 'holm', alpha = 0.000001)
row.names(res)=dataset$Gene.Name
summary(res)

# Bien, como podemos ver, el ajuste de las pruebas múltiples de Holm es mucho más sólido y ha eliminado casi 100 genes previamente marcados 
# como significativos. Ahora solo quedan 57 con diferencias significativas entre las muestras de cáncer y las normales 
# (46 regulados positivamente o 0,11 %, 11 regulados negativamente o 0,027 %). Ahora, como puedes ver, una de las experiencias del análisis 
# de expresión diferencial de ARN es encontrar tan solo un 0,027 % de genes regulados negativamente interesantes, lo que sería un poco más 
# de 1 en 3000 genes.

# Conclusión clave: cuando se estudia una gran cantidad de genes, como decenas de miles, las pruebas deben ser muy conservadoras 
# en términos de múltiples ajustes de pruebas.

######################################################################################################################################

# Crear un gráfico de volcán de grado de publicación con genes de interés marcados

# EnhancedVolcano: Función que genera un gráfico de volcán donde se visualizan los genes en función de su cambio de expresión 
# (log2FoldChange) y su significancia estadística (pvalue).

# lab = dataset$Gene.Name: Etiqueta los puntos en el gráfico con los nombres de los genes.
# x = 'log2FoldChange': El eje X representa el cambio logarítmico de la expresión entre las dos condiciones.
# y = 'pvalue': El eje Y representa el valor p de cada gen, lo que refleja su significancia.
# pCutoff = 10e-5: Se utiliza un umbral de valor p de 1e-5 para resaltar los genes más significativos.
# FCcutoff = 1.333: Define un umbral de cambio de expresión (fold change) para destacar genes con al menos 1.333 veces de cambio en su expresión.
# drawConnectors = TRUE: Dibuja conectores desde las etiquetas de los genes hacia los puntos en el gráfico.

EnhancedVolcano(res,
                lab = dataset$Gene.Name,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-5,
                FCcutoff = 1.333,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'The results',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=1.333; p value cutof=10e-5',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)

# El principio de interpretación de los gráficos de volcanes es bastante intuitivo. El gráfico está segmentado en áreas que tienen 
# observaciones significativas y no significativas según los valores p y los cambios log2fold.

######################################################################################################################################

# Mejor trama del volcán

# Esta sección repite el gráfico de volcán pero ajusta las condiciones:

# y = 'padj': Utiliza los valores p ajustados (padj) en lugar de los valores p sin corregir.
# pCutoff = 10e-7: Umbral de significancia más estricto (1e-7).
# FCcutoff = 2.5: Umbral más alto de cambio de expresión (fold change ≥ 2.5).


# Crear volcanoplot de grado de publicación con genes de interés marcados
EnhancedVolcano(res,
                lab = dataset$Gene.Name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-7,
                FCcutoff = 2.5,
                xlim = c(-5.7, 5.7),
                ylim = c(0, -log10(10.2e-12)),
                pointSize = 1.3,
                labSize = 2.6,
                title = 'The results',
                subtitle = 'Differential expression analysis',
                caption = 'log2fc cutoff=1.333; p value cutof=10e-6',
                legendPosition = "right",
                legendLabSize = 14,
                col = c('lightblue', 'orange', 'blue', 'red2'),
                colAlpha = 0.6,
                drawConnectors = TRUE,
                hline = c(10e-8),
                widthConnectors = 0.5)

# Conclusión clave: generalmente es mejor utilizar padj para el eje y cuando se trabaja con una gran cantidad de transcripciones.

######################################################################################################################################

# Crear el marco de datos final que consiste en resultados deseq ordenados basados en log2fc

# Create the final dataframe consisting of ordered DESeq results based on log2FoldChange
resord = as.data.frame(res)
finaltable = cbind(dataset$Gene.Name, resord)
finaltable = finaltable[order(finaltable$log2FoldChange), ]
write.table(finaltable, file = '/Users/usuario/Documents/programas/oncology_diff_express/finaltable_ori.csv', sep = ',', col.names = NA)

# Filtrar genes con cambios de expresión significativos (positivos o negativos)
# Usamos el criterio p-value ajustado < 10e-7 y log2FoldChange ≥ 2.5 o ≤ -2.5
significativos <- resord[which(resord$padj < 10e-7 & (resord$log2FoldChange >= 2.5 | resord$log2FoldChange <= -2.5)), ]

# Asegurarse de que los nombres de los genes estén alineados con las filas de resultados significativos
# Filtramos solo aquellos genes que están en `significativos`
significativos_table <- significativos
significativos_table$Gene.Name <- rownames(significativos)

# Renombrar la columna de los nombres de los genes para mayor claridad
colnames(significativos_table)[ncol(significativos_table)] <- "Gene.Name"

# Guardar el archivo CSV con los genes significativos
write.table(significativos_table, file = '/Users/usuario/Documents/programas/oncology_diff_express/significativo.csv', sep = ',', row.names = FALSE, col.names = TRUE)

######################################################################################################################################

# heatmap
# Cargar la tabla de genes significativos (assume that this is the filtered data)
significativos_table <- read.csv('/Users/usuario/Documents/programas/oncology_diff_express/significativo.csv')

# Seleccionar las columnas que nos interesan (Gene.Name, log2FoldChange)
# Asumiendo que log2FoldChange está en la columna log2FoldChange
gene_log2FC <- data.frame(Gene.Name = significativos_table$Gene.Name, log2FoldChange = significativos_table$log2FoldChange)

# Usamos los valores de log2FoldChange como una matriz para el heatmap
# Como no tenemos múltiples muestras, el heatmap será unidimensional basado en los valores de log2FoldChange

# Crear una matriz con Gene.Name como filas y log2FoldChange como valores
matriz_log2FC <- as.matrix(gene_log2FC$log2FoldChange)
rownames(matriz_log2FC) <- gene_log2FC$Gene.Name

# Crear el heatmap usando pheatmap
pheatmap(matriz_log2FC,
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Colores desde azul (bajo) a rojo (alto)
         border_color = NA,  # Sin bordes entre celdas
         cluster_cols = FALSE,  # No agrupar columnas, ya que solo hay una
         cluster_rows = TRUE,  # Agrupar filas (genes)
         main = "Heatmap de log2FoldChange de Genes Significativos",
         fontsize_row = 8,  # Tamaño de la fuente para los nombres de los genes
         show_colnames = FALSE  # Ocultar los nombres de columnas
)

######################################################################################################################################

# La última parte del tutorial es para practicar la interpretación.

# Los genes más interesantes se encuentran en la parte superior. La interpretación final sería que los genes interesantes se 
# identificaron comenzando con COL11A1, CEMIP, ADAM12, MMP1, OTOP3 y se necesitaría más investigación experimental para concluir 
# si alguno de ellos es un objetivo potencial para el tratamiento del cáncer.

######################################################################################################################################

# Conclusión final: DESeq es muy práctico y eficaz para tomar grandes conjuntos de datos complejos e identificar posibles 
# biomarcadores para la enfermedad de interés. 

# El análisis de expresión diferencial es específico para tener grandes cantidades de pruebas y los valores p ajustados desempeñan 
# un papel clave en la búsqueda de transcripciones expresadas diferencialmente.

# Y que mas podemos hacer

# analisis de enriquecimiento de genes:
# 1. https://maayanlab.cloud/Enrichr/
# 2. https://string-db.org/
# 3. https://ngdc.cncb.ac.cn/databasecommons/database/id/3061
# 4. Cytoscape



