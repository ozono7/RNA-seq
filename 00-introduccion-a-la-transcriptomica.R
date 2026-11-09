# ==========================================================================
# Lección 00. Introducción a la transcriptómica: qué mide y un primer
#   análisis de exploración con datos reales
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
# Serie RNA-seq — Introducción a la transcriptómica
#
# Descripción: Primera lección de transcriptómica de la serie. Antes de
#   buscar en bases de datos (Módulo 01) y antes de cualquier análisis
#   diferencial (Módulo 07, DESeq2 completo), aquí se responde la pregunta
#   básica: ¿qué es la transcriptómica, qué mide y cómo se ve un conjunto
#   de datos de RNA-seq? Se trabaja con datos reales del conjunto
#   `airway` de Bioconductor ExperimentData: se carga la matriz de
#   conteos, se explora, se normaliza de forma básica y se hacen tres
#   gráficos de exploración (tamaño de librería, distribución de conteos
#   y análisis de componentes principales, PCA).
#
#   Esta lección NO repite el Módulo 01 (no hay búsqueda en PubMed,
#   GenBank ni UniProt) y NO hace el análisis diferencial completo: no se
#   ejecuta DESeq() ni se calculan genes diferencialmente expresados.
#   Eso llega en el Módulo 07.
#
# Objetivos de aprendizaje:
#   - Explicar qué es el transcriptoma y qué mide la transcriptómica.
#   - Cargar un conjunto de RNA-seq con Bioconductor (paquete `airway`)
#     y leer su matriz de conteos.
#   - Calcular el tamaño de librería de cada muestra a partir de los
#     conteos y entender por qué las muestras no se comparan en crudo.
#   - Normalizar los conteos con los factores de tamaño de DESeq2 y
#     transformar con vst para explorar los datos.
#   - Leer un PCA de muestras: qué es cada eje y qué patrón buscar.
#
# Requisitos: R con Bioconductor y los paquetes airway, DESeq2,
#   SummarizedExperiment y ggplot2. Se instalan con:
#   BiocManager::install(c("airway", "DESeq2", "SummarizedExperiment",
#                          "ggplot2"))
#
# Organismo de trabajo de la serie: Chlorella vulgaris (microalga verde,
#   TaxID 3077), con énfasis en la acumulación de lípidos y la
#   biorremediación. El conjunto `airway` es de células humanas y se usa
#   aquí porque es el conjunto clásico de ejemplo de Bioconductor para
#   DESeq2; la estructura de los datos y los pasos son los mismos que se
#   aplicarán después a Chlorella.
# ==========================================================================

# --------------------------------------------------------------------------
# 1. Qué es la transcriptómica
# --------------------------------------------------------------------------
# El genoma es el conjunto completo de genes de un organismo: el catálogo
# de lo que una célula PUEDE hacer. El transcriptoma es el conjunto
# completo de moléculas de ARN presentes en una célula, tejido o muestra
# en un momento dado y en una condición dada: lo que la célula ESTÁ
# haciendo en ese momento.
#
# La transcriptómica mide los niveles de expresión de los genes: cuántas
# copias de ARN se producen a partir de cada gen. Un gen con muchas copias
# de ARN se está expresando mucho; un gen con pocas copias, poco. Cuando
# cambia la condición (por ejemplo, luz azul frente a luz roja en
# Chlorella vulgaris, o un tratamiento con un fármaco), algunos genes
# suben su expresión y otros la bajan. Comparar esos niveles entre
# condiciones es el corazón del análisis de RNA-seq.
#
# El RNA-seq secuencia el ARN de cada muestra y cuenta cuántas lecturas
# (reads) se asignan a cada gen. El resultado es una matriz de conteos:
# cada fila es un gen, cada columna es una muestra, y cada celda es el
# número de lecturas asignadas a ese gen en esa muestra. La unidad de
# cada celda es "número de lecturas" (conteos enteros, sin decimales).
# Toda la lección trabaja sobre esa matriz.

# --------------------------------------------------------------------------
# 2. Cargar paquetes
# --------------------------------------------------------------------------
library(airway)
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)

# --------------------------------------------------------------------------
# 3. Cargar los datos: el conjunto `airway`
# --------------------------------------------------------------------------
# El conjunto `airway` es un experimento de RNA-seq publicado: Himes
# y colaboradores (2014), PLoS One 9(6): e99625, PMID 24926665,
# serie GEO GSE52778. Se secuenció el ARN de cuatro líneas primarias de
# células humanas de músculo liso de las vías respiratorias. Cada línea
# se secuenció dos veces: una sin tratar y una tratada con dexametasona,
# un glucocorticoide sintético (1 micromolar durante 18 horas). En total
# son 8 muestras: 4 líneas celulares x 2 condiciones.
#
# En Bioconductor, estos datos vienen en un objeto
# RangedSummarizedExperiment: la matriz de conteos va en el "assay" y
# la información de las muestras (línea celular, tratamiento) va en
# colData. Es el mismo tipo de objeto que se usará con los datos de
# Chlorella vulgaris en los módulos siguientes.

data(airway)

# Dimensiones de la matriz: filas = genes, columnas = muestras.
# El número exacto de genes depende de la versión del paquete airway;
# no se escribe aquí a mano: se calcula con dim().
dim(airway)

# Información de las muestras: una fila por muestra.
# Columnas: cell (línea celular: N61311, N052611, N080611, N061011),
# dex (tratamiento: "untrt" sin tratar, "trt" tratado con dexametasona),
# albut, Run, avgLength, Experiment, Sample y BioSample.
colData(airway)

# --------------------------------------------------------------------------
# 4. La matriz de conteos
# --------------------------------------------------------------------------
# assay() extrae la matriz de conteos del objeto. Cada celda es el número
# de lecturas asignadas a un gen (fila) en una muestra (columna).
matriz_conteos <- assay(airway, "counts")

# Primeras filas y columnas: los identificadores de fila son códigos
# Ensembl de genes humanos (ENSG...). En Chlorella los identificadores
# serán otros, pero la forma de la tabla es la misma.
head(matriz_conteos[, 1:4])

# ¿Cuántos genes tienen cero lecturas en todas las muestras? Son genes
# que no se expresaron en ninguna muestra de este experimento. El número
# se calcula de la matriz, no se escribe a mano.
genes_sin_lecturas <- sum(rowSums(matriz_conteos) == 0)
genes_sin_lecturas

# --------------------------------------------------------------------------
# 5. Tamaño de librería: cuántas lecturas tiene cada muestra
# --------------------------------------------------------------------------
# El tamaño de librería de una muestra es la suma de todos sus conteos:
# el número total de lecturas secuenciadas en esa muestra. Se calcula
# sumando cada columna de la matriz. Su unidad es "número de lecturas";
# aquí se expresa en millones de lecturas para leerlo con comodidad.
tamano_libreria <- colSums(matriz_conteos)

resumen_muestras <- data.frame(
  muestra = colnames(matriz_conteos),
  linea_celular = colData(airway)$cell,
  tratamiento = colData(airway)$dex,
  tamano_libreria_lecturas = tamano_libreria,
  tamano_libreria_millones = tamano_libreria / 1e6,
  row.names = NULL,
  stringsAsFactors = FALSE
)

print(resumen_muestras)

# Las muestras no tienen el mismo tamaño de librería: unas se
# secuenciaron más a fondo que otras. Por eso los conteos crudos NO se
# comparan directamente entre muestras: una muestra con el doble de
# lecturas totales tendrá, en promedio, el doble de conteos en cada gen
# aunque la biología sea idéntica. La normalización corrige eso.

# --------------------------------------------------------------------------
# 6. Normalización básica con DESeq2 (solo cargar y normalizar)
# --------------------------------------------------------------------------
# Se construye un objeto DESeqDataSet a partir de la matriz de conteos
# y de la información de las muestras. El diseño ~ cell + dex indica
# que cada línea celular se midió en las dos condiciones (diseño
# pareado: la misma línea, sin tratar y tratada). Aquí NO se ejecuta
# DESeq(): no hay análisis diferencial en esta lección, solo
# normalización para poder explorar los datos.
dds <- DESeqDataSetFromMatrix(
  countData = matriz_conteos,
  colData = colData(airway),
  design = ~ cell + dex
)

# Factores de tamaño: un número por muestra, calculado a partir de los
# conteos por el método de la mediana de razones de DESeq2. Una muestra
# con factor de tamaño 2 se secuenció aproximadamente al doble de
# profundidad que una muestra con factor 1. Los factores de tamaño son
# adimensionales (un cociente, sin unidades).
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Conteos normalizados: cada conteo crudo se divide por el factor de
# tamaño de su muestra. Después de normalizar, las muestras quedan en
# la misma escala y sí se pueden comparar entre sí. Los conteos
# normalizados ya no son enteros: son números con decimales.
conteos_normalizados <- counts(dds, normalized = TRUE)
head(conteos_normalizados[, 1:4])

# --------------------------------------------------------------------------
# 7. Gráfico 1. Tamaño de librería por muestra
# --------------------------------------------------------------------------
# Barras: una por muestra. La altura es el tamaño de librería en
# millones de lecturas. El color indica el tratamiento (sin tratar o
# tratado con dexametasona) y las muestras se agrupan por línea celular
# (un panel por línea).
grafico_libreria <- ggplot(resumen_muestras,
                           aes(x = muestra,
                               y = tamano_libreria_millones,
                               fill = tratamiento)) +
  geom_col() +
  facet_wrap(~ linea_celular, scales = "free_x") +
  labs(
    title = "Tamaño de librería por muestra (conjunto airway)",
    x = "Muestra (código SRR)",
    y = "Tamaño de librería (millones de lecturas)",
    fill = "Tratamiento"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(grafico_libreria)

# Interpretación del gráfico 1:
# - Qué muestra: una barra por cada una de las 8 muestras del conjunto
#   airway. Eje x: el código de la muestra (SRR...), agrupada por línea
#   celular. Eje y: tamaño de librería en millones de lecturas (número
#   total de lecturas secuenciadas en esa muestra).
# - Qué patrón se ve: las barras no tienen la misma altura; unas muestras
#   tienen más lecturas totales que otras. Eso es profundidad de
#   secuenciación, no biología: una muestra no expresa más genes por
#   tener más lecturas, simplemente se secuenció más.
# - Qué significa: como los tamaños de librería difieren, los conteos
#   crudos no se comparan entre muestras sin normalizar. La normalización
#   de la sección 6 existe para corregir justo esta diferencia.

# --------------------------------------------------------------------------
# 8. Gráfico 2. Distribución de los conteos normalizados por muestra
# --------------------------------------------------------------------------
# Los conteos van de 0 a decenas de miles, así que se grafican en escala
# log10: log10(conteo normalizado + 1). Se suma 1 porque log10(0) no
# existe y muchos genes tienen conteo cero. La unidad del eje es
# log10(conteos normalizados + 1), una escala logarítmica sin unidades
# físicas: cada unidad del eje es un factor de 10 en conteos.
conteos_log <- log10(conteos_normalizados + 1)

datos_largos <- data.frame(
  muestra = rep(colnames(conteos_log), each = nrow(conteos_log)),
  tratamiento = rep(colData(airway)$dex, each = nrow(conteos_log)),
  conteo_log10 = as.vector(conteos_log),
  stringsAsFactors = FALSE
)

grafico_distribucion <- ggplot(datos_largos,
                               aes(x = muestra,
                                   y = conteo_log10,
                                   fill = tratamiento)) +
  geom_boxplot(outlier.shape = NA) +
  labs(
    title = "Distribución de conteos normalizados por muestra (conjunto airway)",
    x = "Muestra (código SRR)",
    y = "log10(conteo normalizado + 1)",
    fill = "Tratamiento"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(grafico_distribucion)

# Interpretación del gráfico 2:
# - Qué muestra: un diagrama de caja por muestra con la distribución de
#   los conteos normalizados de todos los genes, en escala log10.
#   Eje x: el código de la muestra. Eje y: log10(conteo normalizado + 1),
#   escala logarítmica adimensional; un valor de 3 en este eje
#   corresponde a unos 1000 conteos normalizados.
# - Qué patrón se ve: después de normalizar, las cajas de las 8 muestras
#   quedan a alturas parecidas, con la mediana en el mismo rango. La
#   mayoría de los genes tiene conteos bajos y unos pocos tienen conteos
#   muy altos; por eso la distribución es asimétrica y la escala log10
#   es necesaria para verla.
# - Qué significa: que las cajas queden alineadas confirma que la
#   normalización puso a las muestras en la misma escala. Si una muestra
#   quedara muy desplazada respecto a las demás, habría que revisarla
#   antes de seguir: podría tener un problema de calidad o de
#   secuenciación.

# --------------------------------------------------------------------------
# 9. Gráfico 3. PCA de las muestras
# --------------------------------------------------------------------------
# Para el PCA primero se estabiliza la varianza con vst (variance
# stabilizing transformation) de DESeq2: comprime los conteos altos y
# expande los bajos, para que los genes muy expresados no dominen el
# análisis. blind = FALSE usa el diseño del experimento (~ cell + dex).
# Después se hace el PCA con prcomp sobre las muestras (las muestras son
# las filas y los genes las columnas de la matriz transpuesta).
vsd <- vst(dds, blind = FALSE)

pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = FALSE)

# Porcentaje de varianza que explica cada componente: se calcula de los
# datos (desviaciones estándar del PCA). La suma de los porcentajes de
# todas las componentes es 100.
varianza_explicada <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

datos_pca <- data.frame(
  muestra = rownames(pca$x),
  linea_celular = colData(airway)$cell,
  tratamiento = colData(airway)$dex,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

print(datos_pca)
print(varianza_explicada[1:2])

grafico_pca <- ggplot(datos_pca,
                      aes(x = PC1, y = PC2,
                          color = tratamiento,
                          shape = linea_celular)) +
  geom_point(size = 4) +
  labs(
    title = "PCA de las muestras (conjunto airway, conteos con vst)",
    x = paste0("PC1 (", round(varianza_explicada[1], 1),
               "% de la varianza)"),
    y = paste0("PC2 (", round(varianza_explicada[2], 1),
               "% de la varianza)"),
    color = "Tratamiento",
    shape = "Línea celular"
  ) +
  theme_minimal()

print(grafico_pca)

# Interpretación del gráfico 3:
# - Qué muestra: un punto por muestra. Eje x: primera componente
#   principal (PC1); eje y: segunda componente principal (PC2). Los
#   ejes del PCA son adimensionales: cada componente es una combinación
#   de la expresión de miles de genes, no una magnitud medida en el
#   laboratorio, por eso no llevan unidades. El porcentaje entre
#   paréntesis es la parte de la variación total entre muestras que
#   resume esa componente, calculado de los datos.
# - Qué patrón buscar: puntos cercanos = muestras con perfiles de
#   expresión parecidos; puntos alejados = perfiles distintos. Si las
#   muestras tratadas con dexametasona quedan a un lado y las no
#   tratadas al otro, el tratamiento es la mayor fuente de variación
#   entre muestras. Si, en cambio, las muestras se agrupan por línea
#   celular (mismo símbolo), la línea celular pesa más que el
#   tratamiento. En un experimento bien diseñado, las réplicas de la
#   misma condición quedan juntas.
# - Qué significa: el PCA es un control de calidad antes de analizar.
#   Si una muestra quedara lejos de todas las demás, habría que
#   revisar si es un error de etiquetado o un problema técnico. En
#   Chlorella vulgaris el mismo gráfico se leerá igual: cada punto una
#   muestra (por ejemplo, un color de luz), y lo primero que se mira es
#   si las condiciones se separan.

# --------------------------------------------------------------------------
# 10. Del ejemplo humano a Chlorella vulgaris
# --------------------------------------------------------------------------
# El conjunto airway es de células humanas, pero la lección es la misma
# para Chlorella vulgaris: el RNA-seq de Chlorella también produce una
# matriz de conteos (filas = genes de Chlorella, columnas = muestras,
# por ejemplo luz azul, luz roja y luz blanca), también hay que mirar el
# tamaño de librería de cada muestra, también hay que normalizar antes
# de comparar, y el PCA también se lee buscando si las muestras se
# agrupan por condición de luz. Lo único que cambia es el organismo y
# la pregunta biológica: en Chlorella la pregunta es qué genes cambian
# su expresión con el color de la luz y cómo eso se relaciona con la
# acumulación de lípidos y la biorremediación. El análisis diferencial
# que responde esa pregunta es el Módulo 07.

# --------------------------------------------------------------------------
# 11. Guardar resultados
# --------------------------------------------------------------------------
# Resumen de muestras con el tamaño de librería calculado de la matriz.
write.csv(resumen_muestras, "00_resumen_muestras.csv",
          row.names = FALSE)

# Coordenadas del PCA por muestra, con línea celular y tratamiento.
write.csv(datos_pca, "00_pca_muestras.csv", row.names = FALSE)

# Archivos generados:
#   - 00_resumen_muestras.csv : muestra, línea celular, tratamiento y
#     tamaño de librería (en lecturas y en millones de lecturas).
#   - 00_pca_muestras.csv     : coordenadas PC1 y PC2 de cada muestra.

# --------------------------------------------------------------------------
# 12. Exportar los gráficos
# --------------------------------------------------------------------------
# Los tres gráficos se exportan al final del script, en PNG a 300 dpi.
ggsave("00_tamano_libreria.png", plot = grafico_libreria,
       width = 8, height = 5, dpi = 300)
ggsave("00_distribucion_conteos.png", plot = grafico_distribucion,
       width = 8, height = 5, dpi = 300)
ggsave("00_pca_muestras.png", plot = grafico_pca,
       width = 7, height = 5, dpi = 300)

# --------------------------------------------------------------------------
# Qué aprendimos
# --------------------------------------------------------------------------
# - La transcriptómica mide los niveles de expresión de los genes: cuántas
#   copias de ARN produce cada gen en cada muestra y en cada condición.
# - El resultado de un experimento de RNA-seq es una matriz de conteos:
#   filas = genes, columnas = muestras, celdas = número de lecturas
#   asignadas a cada gen en cada muestra.
# - El conjunto `airway` de Bioconductor es un experimento real de
#   RNA-seq (Himes et al. 2014): 8 muestras, 4 líneas celulares humanas
#   de músculo liso, cada una sin tratar y tratada con dexametasona.
# - El tamaño de librería (suma de conteos por muestra) difiere entre
#   muestras; por eso los conteos crudos no se comparan sin normalizar.
# - La normalización con los factores de tamaño de DESeq2 pone a las
#   muestras en la misma escala; la transformación vst prepara los
#   conteos para la exploración.
# - El PCA resume la variación entre muestras en dos ejes adimensionales
#   (PC1 y PC2): muestras cercanas tienen perfiles de expresión
#   parecidos, y lo primero que se mira es si las condiciones se
#   separan.
# - Los mismos pasos se aplicarán a Chlorella vulgaris: matriz de
#   conteos, tamaño de librería, normalización y PCA por condición de
#   luz. El análisis diferencial completo es el Módulo 07.
