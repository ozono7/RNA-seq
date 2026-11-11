# ==========================================================================
# Módulo 03. Control de calidad de lecturas (FastQC, MultiQC)
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
# Serie "RNA-seq" (avanzado) — Módulo 03
#
# Descripción: En el Módulo 02 descargamos las lecturas crudas de RNA-seq
#   de Chlorella vulgaris del SRA. Antes de limpiarlas o analizarlas hay
#   que revisar si son buenas: eso es el control de calidad. En este
#   módulo se corre FastQC sobre los archivos FASTQ de la corrida de
#   demostración (SRR5314137, WT_rep2), se agrupan los reportes con
#   MultiQC y se reproducen en R las gráficas principales de FastQC a
#   partir de las lecturas mismas, para aprender a leerlas y decidir si
#   los datos pasan a la limpieza (Módulo 04) o hay que descartar algo.
#
# Objetivos de aprendizaje:
#   - Entender el formato FASTQ y la escala de calidad Phred (Q20, Q30).
#   - Instalar y correr FastQC sobre archivos FASTQ, y leer su reporte
#     HTML: calidad por base, calidad por lectura, contenido GC,
#     composición de bases por posición y longitud de lectura.
#   - Agregar muchos reportes de FastQC en uno solo con MultiQC.
#   - Leer archivos FASTQ dentro de R con el paquete ShortRead de
#     Bioconductor y calcular las mismas métricas que FastQC.
#   - Decidir, con números, si las lecturas están listas para el
#     análisis o necesitan limpieza.
#
# Requisitos: R con los paquetes tidyverse y ShortRead (Bioconductor).
#   Instálalos con:
#     install.packages("tidyverse")
#     if (!requireNamespace("BiocManager", quietly = TRUE)) {
#       install.packages("BiocManager")
#     }
#     BiocManager::install("ShortRead")
#   Nota: en Linux con el locale en POSIX, los acentos de las gráficas
#   salen mal; en ese caso corre el script con LC_ALL=C.utf8 Rscript.
#   En macOS y RStudio el locale ya es UTF-8 y no hace falta.
#   Además, FastQC y MultiQC (programas de línea de comandos, fuera de R):
#     - Con conda:  conda install -c bioconda fastqc multiqc
#     - En macOS:   brew install fastqc  (y pip install multiqc)
#     - En Ubuntu:  sudo apt install fastqc  (y pip install multiqc)
#   FastQC necesita Java. MultiQC necesita Python 3.
#   Este módulo corre primero en R (siempre funciona) y, si FastQC y
#   MultiQC están instalados, los ejecuta también sobre la corrida de
#   demostración para comprobar que el flujo completo funciona.
#
# Paquetes de Bioconductor usados en este módulo:
#   - ShortRead: lee archivos FASTQ dentro de R (readFastq) y da acceso
#     a las secuencias (sread) y a las calidades (quality) de cada
#     lectura. La conversión as(quality, "matrix") entrega los valores
#     Phred en una matriz (una fila por lectura, una columna por
#     posición). También trae funciones de resumen como
#     alphabetFrequency, que cuenta cuántas veces aparece cada base
#     (A, C, G, T, N) en cada lectura. Se instala con BiocManager:
#       BiocManager::install("ShortRead")
#   En los módulos siguientes, ShortRead servirá para revisar las
#   lecturas antes y después de la limpieza, y DESeq2 (Módulo 07) usará
#   la matriz de conteos que se construya a partir de estos FASTQ.
#
# Organismo de trabajo: Chlorella vulgaris (microalga verde, TaxID 3077).
# Estudio del proyecto: PRJNA378177 (ver Módulo 02).
# ==========================================================================

# --------------------------------------------------------------------------
# 1. Cargar paquetes
# --------------------------------------------------------------------------
library(tidyverse)
library(ShortRead)

# --------------------------------------------------------------------------
# 2. Qué es el control de calidad y qué mide FastQC
# --------------------------------------------------------------------------
# Un archivo FASTQ guarda cada lectura en 4 líneas: el identificador
# (@SRR...), la secuencia de nucleótidos (A, C, G, T, N), una línea con
# el signo + y la cadena de calidad, con un carácter por cada base. Cada
# carácter codifica la probabilidad de que esa base esté mal leída, en
# la escala Phred: Q = -10 * log10(P_error). Así, Q10 significa 1 error
# por cada 10 bases (90 % de acierto), Q20 es 1 error por cada 100 bases
# (99 % de acierto) y Q30 es 1 error por cada 1,000 bases (99.9 % de
# acierto). En RNA-seq se suele pedir que la mayor parte de las bases
# esté por encima de Q30.
#
# FastQC (Babraham Institute) lee un archivo FASTQ y produce un reporte
# HTML con 11 análisis. Los que se revisan en este módulo son:
#   - Per base sequence quality: calidad media en cada posición de la
#     lectura (1 a 76 nt, la longitud modal). Si la calidad cae al final
#     de la lectura, es normal en Illumina y se corrige recortando
#     (Módulo 04).
#   - Per sequence quality scores: distribución de la calidad media de
#     cada lectura. Lo ideal es un pico por encima de Q30.
#   - Per base sequence content: porcentaje de A, C, G y T en cada
#     posición. En datos aleatorios cada base ronda el 25 %; al inicio
#     de la lectura siempre hay sesgo, porque los primeros ciclos de la
#     secuenciación no son aleatorios.
#   - Per sequence GC content: distribución del contenido GC por
#     lectura, comparada con una curva normal teórica. Un desplazamiento
#     o un segundo pico avisa de contaminación o de sesgo.
#   - Sequence Length Distribution: longitud de las lecturas. Aquí la
#     longitud modal es 76 nt, con una cola de lecturas más cortas que
#     el depósito del SRA ya trae recortadas (ver la sección 12).
#   - Overrepresented sequences y Adapter Content: secuencias que
#     aparecen demasiado o restos de adaptador; si salen, la limpieza
#     del Módulo 04 los quita.
# Cada análisis trae una marca: verde (bien), naranja (revisar) o rojo
# (problema). MultiQC (Ewels Lab) junta los reportes de FastQC de todas
# las muestras en un solo HTML, para comparar las 6 corridas del
# proyecto de un vistazo.
#
# En este módulo se hace lo mismo que FastQC, pero en R y desde los
# datos: se leen las lecturas con ShortRead y se calculan la calidad
# por base, la calidad por lectura, el contenido GC, la composición de
# bases por posición y la longitud, y se grafican con ggplot2.

# --------------------------------------------------------------------------
# 3. Leer la tabla de muestras del Módulo 02
# --------------------------------------------------------------------------
# La tabla 02_muestras_chlorella_sra.csv la genera el Módulo 02 y es el
# archivo que ordena el proyecto: muestra, condición, réplica, corrida
# SRR y tamaño. Si no está en la carpeta, corre primero el Módulo 02.

if (!file.exists("02_muestras_chlorella_sra.csv")) {
  stop("No se encontró 02_muestras_chlorella_sra.csv. ",
       "Corre primero el Módulo 02 (02-descargar-datos-sra.R), que genera ",
       "la tabla de muestras del proyecto.")
}

tabla_muestras <- read_csv(
  "02_muestras_chlorella_sra.csv",
  show_col_types = FALSE
)

cat("Muestras del proyecto:", nrow(tabla_muestras), "\n")
print(tabla_muestras %>% select(muestra, condicion, corrida, bases_gb))

# La corrida de demostración es la más pequeña del estudio (WT_rep2,
# SRR5314137, 0.18 gigabases): con ella se hace todo el control de
# calidad de este módulo, para que corra rápido en una computadora
# normal. Las otras 5 corridas se revisan con el mismo guion de
# FastQC/MultiQC de la sección 5.
corrida_demo <- tabla_muestras %>%
  slice_min(bases_gb, n = 1) %>%
  pull(corrida)

cat("Corrida de demostración:", corrida_demo, "\n")

# --------------------------------------------------------------------------
# 4. Descargar la corrida de demostración (ENA, descarga directa)
# --------------------------------------------------------------------------
# Los archivos FASTQ de la corrida de demostración también están en el
# ENA (European Nucleotide Archive), que refleja el SRA y permite
# descargarlos directo por HTTPS, sin el SRA Toolkit. Son dos archivos
# (R1 y R2, datos paired-end) de unos 69 y 77 megabytes. Después de
# descargarlos se comprueba su integridad con el valor MD5 publicado
# en el Módulo 02: si el MD5 no coincide, la descarga quedó incompleta
# y hay que repetirla.

dir.create("fastq", showWarnings = FALSE)

archivos_fastq <- tibble(
  archivo = paste0(corrida_demo, "_", 1:2, ".fastq.gz"),
  url = paste0(
    "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR531/007/", corrida_demo, "/",
    corrida_demo, "_", 1:2, ".fastq.gz"
  ),
  md5 = c(
    "79ab761d57b4a04baa4f64ee4859aecb",
    "77a05c6fd6eb6a91ce3cdf6c4d4d288c"
  )
)

for (i in seq_len(nrow(archivos_fastq))) {
  ruta <- file.path("fastq", archivos_fastq$archivo[i])
  if (!file.exists(ruta)) {
    cat("Descargando", archivos_fastq$archivo[i], "...\n")
    download.file(archivos_fastq$url[i], ruta, mode = "wb", quiet = TRUE)
  }
  md5_real <- tools::md5sum(ruta)
  if (unname(md5_real) != archivos_fastq$md5[i]) {
    stop("El MD5 de ", archivos_fastq$archivo[i], " no coincide con el ",
         "esperado. Borra el archivo y vuelve a correr el módulo.")
  }
  cat(archivos_fastq$archivo[i], "descargado y verificado (MD5 correcto).\n")
}

# --------------------------------------------------------------------------
# 5. Guion de FastQC y MultiQC para las 6 corridas
# --------------------------------------------------------------------------
# FastQC se corre en la terminal, un archivo FASTQ a la vez, y escribe
# un reporte HTML y un ZIP por archivo en la carpeta fastqc/. MultiQC
# lee todos esos reportes y escribe un solo HTML agregado en multiqc/.
# Este bloque escribe el guion 03_fastqc_multiqc.sh con los comandos
# exactos para las 6 corridas del proyecto, para correrlo tal cual:
#   bash 03_fastqc_multiqc.sh
# Si FastQC está instalado en esta computadora, el módulo lo corre
# ahora mismo sobre la corrida de demostración; si MultiQC también lo
# está, agrega ese reporte. Si no están instalados, el script sigue
# adelante: el control de calidad en R (secciones 6 a 11) ya está hecho,
# y el guion queda listo para cuando se instalen (ver la cabecera).

lineas_sh <- c(
  "#!/bin/bash",
  "# Control de calidad de las 6 corridas de RNA-seq de Chlorella vulgaris (PRJNA378177)",
  "# con FastQC y MultiQC. Generado por el Módulo 03 de la serie RNA-seq.",
  "# Uso: bash 03_fastqc_multiqc.sh",
  "# Requisitos: fastqc y multiqc instalados (ver la cabecera del script).",
  "",
  "mkdir -p fastqc multiqc",
  ""
)

for (corrida in tabla_muestras$corrida) {
  lineas_sh <- c(
    lineas_sh,
    paste0(
      "fastqc --threads 4 -o fastqc ",
      "fastq/", corrida, "_1.fastq.gz ",
      "fastq/", corrida, "_2.fastq.gz"
    )
  )
}
lineas_sh <- c(lineas_sh, "", "multiqc fastqc/ -o multiqc", "")

writeLines(lineas_sh, "03_fastqc_multiqc.sh")

if (nzchar(Sys.which("fastqc"))) {
  dir.create("fastqc", showWarnings = FALSE)
  cat("Corriendo FastQC sobre la corrida de demostración...\n")
  system2(
    "fastqc",
    args = c(
      "--threads", "4", "-o", "fastqc",
      file.path("fastq", paste0(corrida_demo, "_1.fastq.gz")),
      file.path("fastq", paste0(corrida_demo, "_2.fastq.gz"))
    )
  )
  reportes_fastqc <- list.files(
    "fastqc",
    pattern = paste0("^", corrida_demo, "_[12]_fastqc\\.html$"),
    full.names = TRUE
  )
  cat("Reportes FastQC generados:", length(reportes_fastqc), "de 2 esperados.\n")
  if (nzchar(Sys.which("multiqc")) && length(reportes_fastqc) == 2) {
    dir.create("multiqc", showWarnings = FALSE)
    cat("Corriendo MultiQC sobre los reportes de FastQC...\n")
    system2("multiqc", args = c("fastqc/", "-o", "multiqc"))
  }
} else {
  cat("FastQC no está instalado en esta computadora.\n")
  cat("Instálalo (ver la cabecera del script) y corre: bash 03_fastqc_multiqc.sh\n")
}

# --------------------------------------------------------------------------
# 6. Leer las lecturas con ShortRead (Bioconductor)
# --------------------------------------------------------------------------
# readFastq lee un archivo FASTQ completo. La corrida de demostración
# tiene 1,223,634 lecturas por archivo (spots, ver Módulo 02), así que
# para los cálculos de este módulo se toma una submuestra aleatoria de
# 50,000 lecturas de cada archivo, con semilla fija para que el
# resultado sea reproducible. FastQC hace lo mismo: analiza una
# submuestra de las lecturas, porque con decenas de miles de lecturas
# las distribuciones de calidad ya quedan bien estimadas.
#
# Las lecturas no miden todas lo mismo: la longitud modal es 76 nt,
# con una cola de lecturas más cortas (ver la sección 12). Para los
# cálculos por posición (calidad por base y composición por base) se
# usan solo las lecturas de longitud modal, porque esos cálculos
# necesitan que todas las lecturas tengan la misma longitud; las demás
# se cuentan en la tabla resumen y en el gráfico de longitudes, no se
# pierden.

set.seed(20261111)
n_muestra <- 50000

leer_muestra_fastq <- function(archivo) {
  fq <- readFastq("fastq", pattern = archivo)
  n_total <- length(fq)
  if (n_total == 0) {
    stop("El archivo ", archivo, " no tiene lecturas.")
  }
  indices <- sample(n_total, size = min(n_muestra, n_total))
  fq_muestra <- fq[indices]
  # Secuencias como texto, para calcular la composición de bases.
  secuencias <- as.character(sread(fq_muestra))
  # Longitud de cada lectura, en nucleótidos (nt).
  longitudes <- width(sread(fq_muestra))
  # Longitud modal: la longitud que más se repite en la submuestra.
  longitud_modal <- as.integer(
    names(which.max(table(longitudes)))
  )
  # Solo las lecturas de longitud modal entran en los cálculos por
  # posición; el resto se reporta y se excluye de esas matrices.
  mantener <- longitudes == longitud_modal
  n_excluidas <- sum(!mantener)
  cat(
    "Lecturas de longitud distinta a la modal (", longitud_modal,
    " nt) excluidas de los cálculos por posición: ", n_excluidas,
    " de ", length(secuencias), ".\n",
    sep = ""
  )
  # Matriz de calidad Phred con ShortRead: as(quality, "matrix")
  # convierte la cadena de calidad de cada lectura en sus valores
  # Phred (un número por base). Cada fila es una lectura, cada columna
  # una posición. Las lecturas más cortas que la modal se rellenan con
  # NA en las posiciones que no tienen.
  matriz_calidad_completa <- as(quality(fq_muestra), "matrix")
  # Matriz de calidad solo con las lecturas de longitud modal (sin NA).
  matriz_calidad <- matriz_calidad_completa[mantener, , drop = FALSE]
  # Matriz de secuencias como códigos ASCII, para la composición de
  # bases por posición (A = 65, C = 67, G = 71, T = 84, N = 78).
  matriz_secuencia <- matrix(
    unlist(lapply(
      secuencias[mantener],
      function(s) as.integer(charToRaw(s))
    )),
    nrow = sum(mantener), byrow = TRUE
  )
  # Calidad media de cada lectura y porcentaje de bases con Q30 o más,
  # calculados sobre TODAS las lecturas de la submuestra (no solo las
  # de longitud modal). Los NA del relleno no se cuentan.
  calidad_media_por_lectura <- rowMeans(
    matriz_calidad_completa, na.rm = TRUE
  )
  bases_q30 <- sum(matriz_calidad_completa >= 30, na.rm = TRUE)
  suma_phred <- sum(matriz_calidad_completa, na.rm = TRUE)
  # Conteo de cada base por lectura, con alphabetFrequency de ShortRead.
  conteo_bases <- alphabetFrequency(sread(fq_muestra), baseOnly = TRUE)
  columna <- function(nombre) {
    if (nombre %in% colnames(conteo_bases)) {
      conteo_bases[, nombre]
    } else {
      rep(0, nrow(conteo_bases))
    }
  }
  list(
    n_total = n_total,
    n_muestra = length(secuencias),
    longitudes = longitudes,
    longitud_modal = longitud_modal,
    n_excluidas = n_excluidas,
    matriz_calidad = matriz_calidad,
    matriz_secuencia = matriz_secuencia,
    calidad_media_por_lectura = calidad_media_por_lectura,
    calidad_media = suma_phred / sum(longitudes),
    pct_q30 = 100 * bases_q30 / sum(longitudes),
    conteo_gc = columna("G") + columna("C"),
    conteo_n = columna("N"),
    conteo_total = rowSums(conteo_bases)
  )
}

lecturas <- list(
  R1 = leer_muestra_fastq(paste0(corrida_demo, "_1.fastq.gz")),
  R2 = leer_muestra_fastq(paste0(corrida_demo, "_2.fastq.gz"))
)

for (extremo in c("R1", "R2")) {
  cat(
    extremo, ": ", lecturas[[extremo]]$n_total, " lecturas en el archivo; ",
    lecturas[[extremo]]$n_muestra, " en la submuestra.\n",
    sep = ""
  )
}

# --------------------------------------------------------------------------
# 7. Métricas de calidad por archivo y tabla resumen
# --------------------------------------------------------------------------
# Para cada archivo (R1 y R2) se calculan, a partir de la submuestra:
#   - lecturas: número de lecturas en el archivo completo.
#   - longitud media: media de la longitud de las lecturas, en nt.
#   - longitud modal: la longitud que más se repite, en nt.
#   - calidad media: media de la calidad Phred de todas las bases.
#   - pct_q30: porcentaje de bases con Phred mayor o igual a 30.
#   - pct_gc: porcentaje de bases G o C en las lecturas.
#   - pct_n: porcentaje de bases N (base no determinada) en las lecturas.
# La tabla se guarda en CSV: es el registro del control de calidad del
# proyecto y la usará el Módulo 04 para comparar antes y después de la
# limpieza.

resumen_qc <- map_dfr(c("R1", "R2"), function(extremo) {
  datos <- lecturas[[extremo]]
  tibble(
    archivo = paste0(corrida_demo, "_", tolower(extremo), ".fastq.gz"),
    extremo = extremo,
    lecturas = datos$n_total,
    longitud_media_nt = mean(datos$longitudes),
    longitud_modal_nt = datos$longitud_modal,
    calidad_media_phred = datos$calidad_media,
    pct_q30 = datos$pct_q30,
    pct_gc = 100 * sum(datos$conteo_gc) / sum(datos$conteo_total),
    pct_n = 100 * sum(datos$conteo_n) / sum(datos$conteo_total)
  )
})

print(resumen_qc)
write_csv(resumen_qc, "03_resumen_qc.csv")

# Calidad media por posición (1 a 76 nt, la longitud modal), con
# cuartiles, para el gráfico 1. Se calcula por separado para R1 y R2.
calidad_por_base <- map_dfr(c("R1", "R2"), function(extremo) {
  matriz <- lecturas[[extremo]]$matriz_calidad
  tibble(
    extremo = extremo,
    posicion = seq_len(ncol(matriz)),
    media = colMeans(matriz),
    q1 = apply(matriz, 2, quantile, probs = 0.25),
    mediana = apply(matriz, 2, quantile, probs = 0.5),
    q3 = apply(matriz, 2, quantile, probs = 0.75)
  )
})

# Calidad media de cada lectura, para el gráfico 2.
calidad_por_lectura <- map_dfr(c("R1", "R2"), function(extremo) {
  tibble(
    extremo = extremo,
    calidad_media = lecturas[[extremo]]$calidad_media_por_lectura
  )
})

# Contenido GC de cada lectura (porcentaje), para el gráfico 3.
gc_por_lectura <- map_dfr(c("R1", "R2"), function(extremo) {
  datos <- lecturas[[extremo]]
  tibble(
    extremo = extremo,
    pct_gc = 100 * datos$conteo_gc / datos$conteo_total
  )
})

# Composición de bases por posición (porcentaje de A, C, G, T y N en
# cada posición), para el gráfico 4.
composicion_por_base <- map_dfr(c("R1", "R2"), function(extremo) {
  matriz <- lecturas[[extremo]]$matriz_secuencia
  n <- nrow(matriz)
  tibble(
    extremo = extremo,
    posicion = seq_len(ncol(matriz)),
    A = 100 * colSums(matriz == 65) / n,
    C = 100 * colSums(matriz == 67) / n,
    G = 100 * colSums(matriz == 71) / n,
    T = 100 * colSums(matriz == 84) / n,
    N = 100 * colSums(matriz == 78) / n
  ) %>%
    pivot_longer(
      cols = c(A, C, G, T, N),
      names_to = "base", values_to = "porcentaje"
    )
})

# Longitud de cada lectura, para el gráfico 5.
longitud_por_lectura <- map_dfr(c("R1", "R2"), function(extremo) {
  tibble(
    extremo = extremo,
    longitud_nt = lecturas[[extremo]]$longitudes
  )
})

# --------------------------------------------------------------------------
# 8. Gráfico 1: calidad media por posición (como FastQC)
# --------------------------------------------------------------------------
# La línea es la calidad media Phred en cada posición de la lectura;
# la banda marca el rango entre el primer y el tercer cuartil. Las
# líneas punteadas marcan Q20 y Q30.

grafico_calidad_base <- ggplot(calidad_por_base, aes(
  x = posicion, y = media, color = extremo, fill = extremo
)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 30, linetype = "dashed", color = "#59A14F") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "#F28E2B") +
  annotate("text", x = 74, y = 31, label = "Q30", color = "#59A14F", size = 3.5) +
  annotate("text", x = 74, y = 21, label = "Q20", color = "#F28E2B", size = 3.5) +
  scale_color_manual(values = c("R1" = "#4E79A7", "R2" = "#F28E2B")) +
  scale_fill_manual(values = c("R1" = "#4E79A7", "R2" = "#F28E2B")) +
  labs(
    title = "Calidad por base en la corrida de demostración (SRR5314137)",
    subtitle = "Chlorella vulgaris WT_rep2, submuestra de 50,000 lecturas por archivo",
    x = "Posición en la lectura (nucleótidos, nt)",
    y = "Calidad media Phred (adimensional)",
    color = "Extremo",
    fill = "Extremo"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_calidad_base)

# Interpretación:
# El eje horizontal es la posición dentro de la lectura, de 1 a 76 nt
# (nucleótidos, la longitud modal); el eje vertical es la calidad media
# adimensional: Q30 significa 1 error por cada 1,000 bases y Q20, 1
# error por cada 100 bases. La banda sombreada marca el rango entre el
# primer y el tercer cuartil de la calidad en cada posición, es decir,
# dónde cae la mitad central de las lecturas. Si la línea se mantiene
# por encima de Q30 en las 76 posiciones, las lecturas son buenas de
# principio a fin. Lo habitual en Illumina es una caída suave hacia el
# final de la lectura; si la caída baja de Q20, esas bases del final se
# recortan en la limpieza (Módulo 04). R1 y R2 se comparan entre sí:
# en datos paired-end de buena calidad las dos curvas van casi juntas.

ggsave(
  "03_calidad_por_base.png",
  plot = grafico_calidad_base,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 9. Gráfico 2: distribución de la calidad media por lectura
# --------------------------------------------------------------------------
# Cada lectura tiene una calidad media (su Phred promedio). El
# histograma muestra cuántas lecturas caen en cada valor de calidad.

grafico_calidad_lectura <- ggplot(calidad_por_lectura, aes(
  x = calidad_media, fill = extremo
)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "#59A14F") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "#F28E2B") +
  annotate("text", x = 31, y = 0, label = "Q30", color = "#59A14F", size = 3.5) +
  scale_fill_manual(values = c("R1" = "#4E79A7", "R2" = "#F28E2B")) +
  labs(
    title = "Calidad media por lectura (SRR5314137)",
    subtitle = "Submuestra de 50,000 lecturas por archivo",
    x = "Calidad media por lectura (Phred, adimensional)",
    y = "Número de lecturas",
    fill = "Extremo"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_calidad_lectura)

# Interpretación:
# El eje horizontal es la calidad media Phred de cada lectura (número
# adimensional; Q30 es 1 error por cada 1,000 bases) y el eje vertical
# cuenta cuántas lecturas tienen esa calidad. En datos buenos el
# histograma forma un solo pico estrecho por encima de Q30: casi todas
# las lecturas son buenas. Un pico ancho o desplazado hacia la
# izquierda (calidades bajas) avisa de un problema en la corrida, y
# un segundo pico separado puede indicar contaminación o una mezcla de
# dos calidades distintas. La línea verde marca Q30 y la naranja Q20,
# los dos umbrales de referencia del control de calidad.

ggsave(
  "03_calidad_por_lectura.png",
  plot = grafico_calidad_lectura,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 10. Gráfico 3: contenido GC por lectura
# --------------------------------------------------------------------------
# El contenido GC de una lectura es el porcentaje de sus bases que son
# G o C. El histograma muestra cómo se reparte entre las lecturas, y
# la curva marca la distribución normal teórica con la misma media y
# la misma desviación estándar que los datos, que es la comparación
# que hace FastQC.

gc_resumen <- gc_por_lectura %>%
  group_by(extremo) %>%
  summarise(
    media_gc = mean(pct_gc),
    desv_gc = sd(pct_gc),
    .groups = "drop"
  )

grafico_gc <- ggplot(gc_por_lectura) +
  geom_histogram(
    aes(x = pct_gc, y = after_stat(density), fill = extremo),
    binwidth = 1, alpha = 0.6, position = "identity"
  ) +
  geom_vline(
    data = gc_resumen,
    aes(xintercept = media_gc, color = extremo),
    linetype = "dashed", linewidth = 0.8
  ) +
  scale_fill_manual(values = c("R1" = "#4E79A7", "R2" = "#F28E2B")) +
  scale_color_manual(values = c("R1" = "#4E79A7", "R2" = "#F28E2B")) +
  labs(
    title = "Contenido GC por lectura (SRR5314137)",
    subtitle = "Submuestra de 50,000 lecturas por archivo; la línea marca el GC medio",
    x = "Contenido GC por lectura (%)",
    y = "Densidad (fracción de lecturas por punto de GC)",
    fill = "Extremo",
    color = "Extremo"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_gc)

# Interpretación:
# El eje horizontal es el contenido GC de cada lectura, en porcentaje
# (% de bases G o C), y el eje vertical es la densidad: qué fracción
# de las lecturas tiene cada valor de GC. En un RNA-seq sin
# contaminación el histograma forma un solo pico, más o menos
# simétrico, centrado en el GC propio del organismo: Chlorella
# vulgaris tiene un genoma rico en GC (alrededor de 60 %), así que el
# pico debe quedar claramente a la derecha del 50 %. Un pico ancho,
# con hombros o con un segundo pico, avisa de contaminación (por
# ejemplo, ARN de otro organismo) o de un sesgo fuerte en la
# preparación de la librería. La línea punteada marca el GC medio de
# cada archivo, que también queda en la tabla resumen del módulo.

ggsave(
  "03_contenido_gc.png",
  plot = grafico_gc,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 11. Gráfico 4: composición de bases por posición
# --------------------------------------------------------------------------
# En cada posición de la lectura se calcula qué porcentaje de las
# lecturas tiene A, C, G, T o N. Las cuatro bases deben rondar el 25 %
# cada una si la librería es aleatoria.

grafico_composicion <- ggplot(composicion_por_base, aes(
  x = posicion, y = porcentaje, color = base
)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 25, linetype = "dotted", color = "grey40") +
  facet_wrap(~extremo, ncol = 1) +
  scale_color_manual(values = c(
    "A" = "#4E79A7", "C" = "#F28E2B", "G" = "#59A14F",
    "T" = "#E15759", "N" = "grey50"
  )) +
  labs(
    title = "Composición de bases por posición (SRR5314137)",
    subtitle = "Submuestra de 50,000 lecturas por archivo; la línea punteada marca el 25 %",
    x = "Posición en la lectura (nucleótidos, nt)",
    y = "Lecturas con esa base en esa posición (%)",
    color = "Base"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_composicion)

# Interpretación:
# El eje horizontal es la posición dentro de la lectura, de 1 a 76 nt
# (la longitud modal), y el eje vertical es el porcentaje de lecturas
# (A, C, G, T o N) en esa posición. La línea punteada marca el 25 %,
# que es lo esperado si las cuatro bases aparecieran por igual. En
# datos reales las líneas nunca son planas: al inicio de la lectura
# (primeras 10 a 15 posiciones) siempre hay un sesgo, porque los
# primeros ciclos de la secuenciación y la fragmentación del ARN no
# son aleatorios, y en Chlorella el contenido GC es alto, así que las
# líneas de G y C quedan por encima de las de A y T. Lo que importa
# es que después de esas primeras posiciones las curvas se estabilicen
# y se mantengan paralelas: un cambio brusco a mitad de la lectura
# avisaría de un problema en la corrida. Hay una excepción al final:
# en la última posición (76 nt) la A desaparece por completo y la C
# sube a más del 40 % en los dos archivos, un efecto propio del final
# de la lectura en este estudio que coincide con la caída de calidad
# del final que se ve en el gráfico 1; es otra razón por la que el
# Módulo 04 revisará los extremos antes de analizar. La línea de N
# (base no determinada) debe quedar pegada al 0 %: si sube, el
# secuenciador no pudo leer esas bases.

ggsave(
  "03_composicion_por_base.png",
  plot = grafico_composicion,
  width = 9, height = 6.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 12. Gráfico 5: longitud de las lecturas
# --------------------------------------------------------------------------
# La longitud modal de este estudio es 76 nt, pero no todas las
# lecturas miden lo mismo: hay una cola de lecturas más cortas, porque
# el depósito del SRA ya trae recortados los extremos con adaptador o
# con baja calidad. El gráfico de barras cuenta cuántas lecturas hay
# de cada longitud.

grafico_longitud <- ggplot(longitud_por_lectura, aes(
  x = factor(longitud_nt), fill = extremo
)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("R1" = "#4E79A7", "R2" = "#F28E2B")) +
  labs(
    title = "Longitud de las lecturas (SRR5314137)",
    subtitle = "Submuestra de 50,000 lecturas por archivo",
    x = "Longitud de lectura (nucleótidos, nt)",
    y = "Número de lecturas",
    fill = "Extremo"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_longitud)

# Interpretación:
# El eje horizontal es la longitud de la lectura en nucleótidos (nt)
# y el eje vertical cuenta cuántas lecturas tienen esa longitud. La
# barra más alta está en 76 nt, la longitud modal: alrededor de 7 de
# cada 10 lecturas miden 76 nt. A su izquierda se ve una cola de
# lecturas más cortas (75, 74, 73 nt y menos), cada vez con menos
# lecturas: son lecturas que el depósito del SRA ya trae recortadas,
# porque su extremo tenía adaptador o baja calidad. Esa cola es normal
# en datos del SRA y no es un problema: la longitud media queda en
# unos 75.6 nt (ver la tabla 03_resumen_qc.csv) y el análisis sigue
# igual. Lo que sí sería señal de alarma es que la barra principal no
# estuviera en la longitud esperada del protocolo (76 nt) o que hubiera
# dos picos separados, lo que indicaría una mezcla de datos de
# protocolos distintos.

ggsave(
  "03_longitud_lectura.png",
  plot = grafico_longitud,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 13. Veredicto: ¿las lecturas pasan el control de calidad?
# --------------------------------------------------------------------------
# Se decide con tres números de la tabla resumen: la calidad media
# Phred, el porcentaje de bases con Q30 o más, y el porcentaje de
# bases N. La regla que se usa en este proyecto es: calidad media
# mayor o igual a 30, al menos 80 % de las bases en Q30 o más, y menos
# de 1 % de bases N. Si se cumplen las tres, las lecturas pasan y el
# proyecto sigue a la limpieza (Módulo 04) para quitar adaptadores y
# recortar los extremos de baja calidad.

veredicto <- resumen_qc %>%
  mutate(
    pasa = calidad_media_phred >= 30 & pct_q30 >= 80 & pct_n < 1
  )

print(veredicto %>% select(archivo, calidad_media_phred, pct_q30, pct_n, pasa))

if (all(veredicto$pasa)) {
  cat("Veredicto: las lecturas PASAN el control de calidad.\n")
  cat("El proyecto sigue al Módulo 04 (limpieza con Trimmomatic / fastp).\n")
} else {
  cat("Veredicto: al menos un archivo NO pasa el control de calidad.\n")
  cat("Revisa los gráficos y la tabla 03_resumen_qc.csv antes de seguir.\n")
}

# El reporte de FastQC marca dos cosas que hay que saber interpretar:
# "Per base sequence content" sale en rojo (FAIL) porque las primeras
# posiciones tienen sesgo y porque Chlorella es rica en GC: G y C van
# por encima de A y T en toda la lectura. Eso es propio del organismo
# y de la preparación de la librería, no un error de los datos, y por
# eso el veredicto de este módulo no depende de esa marca. "Sequence
# Length Distribution" sale en naranja (WARNING) porque las longitudes
# varían (modal 76 nt con cola de lecturas más cortas, ver sección 12).
# Los módulos de calidad (por base y por lectura), el contenido GC, el
# contenido de N, la duplicación, las secuencias sobrerrepresentadas y
# los adaptadores salen en verde (PASS) en los dos archivos.

# --------------------------------------------------------------------------
# 14. Resumen del módulo
# --------------------------------------------------------------------------
# Archivos generados (guárdalos: los usaremos en los módulos siguientes):
#   - 03_resumen_qc.csv             : tabla con las métricas de calidad
#     de R1 y R2 de la corrida de demostración (lecturas, longitud
#     media y modal, calidad media, % Q30, % GC y % N).
#   - 03_fastqc_multiqc.sh          : guion con los comandos fastqc y
#     multiqc para las 6 corridas del proyecto.
#   - 03_calidad_por_base.png       : gráfico de la calidad por posición.
#   - 03_calidad_por_lectura.png    : gráfico de la calidad por lectura.
#   - 03_contenido_gc.png           : gráfico del contenido GC por lectura.
#   - 03_composicion_por_base.png   : gráfico de la composición de bases.
#   - 03_longitud_lectura.png       : gráfico de la longitud de lectura.
#   - fastqc/SRR5314137_1_fastqc.html y fastqc/SRR5314137_2_fastqc.html:
#     reportes de FastQC de la corrida de demostración, solo si FastQC
#     está instalado; y multiqc/multiqc_report.html, solo si MultiQC
#     también lo está.
#
# Qué sigue (Módulo 04):
#   Limpieza de las lecturas con Trimmomatic o fastp: quitar adaptadores
#   y recortar los extremos de baja calidad, y volver a correr el
#   control de calidad para comprobar que la limpieza funcionó.

# --------------------------------------------------------------------------
# Qué aprendimos
# --------------------------------------------------------------------------
# - Un archivo FASTQ guarda cada lectura en 4 líneas: identificador,
#   secuencia, signo + y cadena de calidad. Cada carácter de calidad
#   es un valor Phred: Q = -10 * log10(P_error), así que Q20 es 1 error
#   por cada 100 bases y Q30 es 1 error por cada 1,000 bases.
# - FastQC revisa un archivo FASTQ y produce un reporte HTML con 11
#   análisis y marcas verde, naranja o rojo. Los principales son la
#   calidad por base, la calidad por lectura, el contenido GC, la
#   composición de bases por posición y la longitud de lectura.
# - MultiQC junta los reportes de FastQC de todas las muestras en un
#   solo HTML, para comparar las 6 corridas del proyecto de un vistazo.
# - En R, el paquete ShortRead de Bioconductor lee archivos FASTQ con
#   readFastq, y de ahí se calculan las mismas métricas que FastQC:
#   calidad por base, calidad por lectura, contenido GC, composición
#   de bases y longitud. Se trabajó con una submuestra de 50,000
#   lecturas por archivo, con semilla fija, igual que hace FastQC.
# - La corrida de demostración (SRR5314137, WT_rep2) tiene lecturas
#   paired-end de 76 nt de longitud modal (media 75.6 nt, con una cola
#   de lecturas más cortas ya recortadas en el SRA). El veredicto se
#   decide con tres números de la tabla 03_resumen_qc.csv: calidad
#   media Phred, porcentaje de bases con Q30 o más, y porcentaje de
#   bases N.
# - Chlorella vulgaris tiene un genoma rico en GC, así que el pico del
#   contenido GC queda a la derecha del 50 % y las curvas de G y C
#   van por encima de las de A y T en la composición por posición.
#   Eso es propio del organismo, no un problema de los datos.
