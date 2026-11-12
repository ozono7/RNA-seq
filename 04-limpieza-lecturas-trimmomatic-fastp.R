# ==========================================================================
# Módulo 04. Limpieza de lecturas (Trimmomatic / fastp)
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
# Serie "RNA-seq" (avanzado) — Módulo 04
#
# Descripción: En el Módulo 03 revisamos la calidad de las lecturas
#   crudas de Chlorella vulgaris (corrida de demostración SRR5314137,
#   WT_rep2) y el veredicto fue que pasan el control de calidad, con
#   dos cosas por corregir en la limpieza: la calidad cae al final de
#   la lectura y la última posición (76 nt) tiene una composición de
#   bases rara. En este módulo se limpian esas lecturas con dos
#   herramientas clásicas, Trimmomatic y fastp: quitan los restos de
#   adaptador, recortan los extremos de baja calidad y descartan las
#   lecturas que quedan demasiado cortas. Después se vuelve a medir la
#   calidad, igual que en el Módulo 03, para comprobar con números que
#   la limpieza funcionó y que no se perdió casi nada.
#
# Objetivos de aprendizaje:
#   - Entender qué es la limpieza de lecturas y por qué va después del
#     control de calidad y antes de cualquier análisis.
#   - Quitar adaptadores de Illumina y explicar qué es un adaptador y
#     por qué aparece al final de las lecturas.
#   - Recortar los extremos de baja calidad con la ventana deslizante
#     (sliding window) y explicar qué hace cada parámetro.
#   - Filtrar por longitud mínima y entender por qué en datos
#     paired-end se conservan los pares (paired) y se separan las
#     lecturas huérfanas (unpaired).
#   - Leer el reporte JSON de fastp y el resumen de Trimmomatic, y
#     comparar el antes y el después con números y gráficas.
#   - Verificar la limpieza en R con el paquete ShortRead de
#     Bioconductor, igual que en el Módulo 03.
#
# Requisitos: R con los paquetes tidyverse, ShortRead (Bioconductor) y
#   jsonlite. Instálalos con:
#     install.packages(c("tidyverse", "jsonlite"))
#     if (!requireNamespace("BiocManager", quietly = TRUE)) {
#       install.packages("BiocManager")
#     }
#     BiocManager::install("ShortRead")
#   Nota: en Linux con el locale en POSIX, los acentos de las gráficas
#   salen mal; en ese caso corre el script con LC_ALL=C.utf8 Rscript.
#   En macOS y RStudio el locale ya es UTF-8 y no hace falta.
#   Además, Trimmomatic y fastp (programas de línea de comandos, fuera
#   de R):
#     - Con conda:  conda install -c bioconda trimmomatic fastp
#     - Trimmomatic necesita Java. Este módulo descarga solo el archivo
#       ZIP oficial de Trimmomatic 0.39 si no encuentra el programa
#       instalado, y lo descomprime en la carpeta trimmomatic/.
#     - fastp ya no publica un binario listo en GitHub: se instala con
#       conda (conda install -c bioconda fastp) o se compila desde el
#       código fuente (https://github.com/OpenGene/fastp). Si fastp no
#       está instalado, el módulo sigue con Trimmomatic y deja escrito
#       el comando de fastp para correrlo después.
#
# Paquetes de Bioconductor usados en este módulo:
#   - ShortRead: lee archivos FASTQ dentro de R (readFastq) y da acceso
#     a las secuencias (sread) y a las calidades (quality) de cada
#     lectura. Se usa para verificar la limpieza de forma independiente
#     de fastp y Trimmomatic: se leen las lecturas crudas y las
#     limpias, y se comparan la longitud, la calidad media y el
#     porcentaje de bases con Q30 o más. Se instala con BiocManager:
#       BiocManager::install("ShortRead")
#   En los módulos siguientes, los archivos FASTQ limpios que se
#   generan aquí son la entrada del ensamblaje de novo con Trinity
#   (Módulo 05) y de la cuantificación con Salmon o Kallisto
#   (Módulo 06).
#
# Organismo de trabajo: Chlorella vulgaris (microalga verde, TaxID 3077).
# Estudio del proyecto: PRJNA378177 (ver Módulo 02).
# ==========================================================================

# --------------------------------------------------------------------------
# 1. Cargar paquetes
# --------------------------------------------------------------------------
library(tidyverse)
library(ShortRead)
library(jsonlite)

# --------------------------------------------------------------------------
# 2. Qué es la limpieza y qué corrige cada herramienta
# --------------------------------------------------------------------------
# Una lectura de RNA-seq trae tres problemas típicos, y la limpieza
# quita los tres antes de analizar:
#   - Restos de adaptador: el adaptador es la secuencia corta que se
#     añade al fragmento de ARN para secuenciarlo. Cuando el fragmento
#     es más corto que la lectura, el secuenciador sigue leyendo y
#     llega al adaptador, que queda pegado al final de la lectura. Esa
#     parte no es ARN de Chlorella y hay que cortarla.
#   - Extremos de baja calidad: en Illumina la calidad baja hacia el
#     final de la lectura (lo vimos en el Módulo 03, gráfico 1). Las
#     bases del final con Phred bajo tienen más probabilidad de estar
#     mal leídas, así que se recortan.
#   - Lecturas demasiado cortas: después de recortar, algunas lecturas
#     quedan tan cortas que ya no sirven para alinear o ensamblar. Se
#     descartan con una longitud mínima (aquí 36 nt).
#
# Trimmomatic (Bolger, Lohse y Usadel, 2014) trabaja por pasos
# encadenados, en el orden en que se escriben:
#   - ILLUMINACLIP: busca el adaptador al inicio y al final de la
#     lectura y lo corta. Los números 2:30:10 son: 2 errores permitidos
#     al alinear el adaptador, 30 puntos de calidad mínima para cortar
#     un adaptador parcial y 10 para cortar un adaptador completo.
#   - SLIDINGWINDOW:4:20: ventana deslizante de 4 bases que avanza por
#     la lectura y corta en cuanto la calidad media de la ventana baja
#     de 20 (Phred). Es el recorte de calidad propiamente dicho.
#   - MINLEN:36: descarta las lecturas que quedan con menos de 36 nt.
#   En datos paired-end, Trimmomatic escribe 4 archivos: los pares
#   que sobreviven juntos (paired) y las lecturas que quedan sin su
#   pareja (unpaired), una por cada extremo. Los pares se conservan
#   juntos porque los pasos siguientes (ensamblaje, cuantificación)
#   los necesitan apareados.
#
# fastp (Chen y colaboradores, 2018) hace lo mismo en un solo programa
# y además escribe un reporte JSON y uno HTML con todo lo medido:
#   - Detecta el adaptador solo, sin archivo de adaptadores.
#   - --cut_front y --cut_tail recortan el inicio y el final con una
#     ventana de 4 bases y calidad media mínima de 20
#     (--cut_window_size 4 --cut_mean_quality 20), igual que el
#     SLIDINGWINDOW de Trimmomatic.
#   - --qualified_quality_phred 20 y --unqualified_percent_limit 40:
#     una lectura se descarta si más del 40 % de sus bases tiene
#     calidad menor de 20.
#   - --length_required 36: descarta las lecturas de menos de 36 nt.
#   Los dos programas llegan al mismo resultado por caminos distintos,
#   así que en este módulo se corren los dos sobre la misma corrida y
#   se comparan sus números.

# --------------------------------------------------------------------------
# 3. Leer la tabla de muestras y el resumen de calidad del Módulo 03
# --------------------------------------------------------------------------
# La tabla 02_muestras_chlorella_sra.csv la genera el Módulo 02 y el
# archivo 03_resumen_qc.csv lo genera el Módulo 03. Si no están en la
# carpeta, corre primero esos módulos: la limpieza necesita saber qué
# corrida es la de demostración y con qué calidad parte.

if (!file.exists("02_muestras_chlorella_sra.csv")) {
  stop("No se encontró 02_muestras_chlorella_sra.csv. ",
       "Corre primero el Módulo 02 (02-descargar-datos-sra.R).")
}
if (!file.exists("03_resumen_qc.csv")) {
  stop("No se encontró 03_resumen_qc.csv. ",
       "Corre primero el Módulo 03 (03-control-calidad-fastqc-multiqc.R).")
}

tabla_muestras <- read_csv(
  "02_muestras_chlorella_sra.csv",
  show_col_types = FALSE
)
resumen_qc_antes <- read_csv("03_resumen_qc.csv", show_col_types = FALSE)

cat("Calidad de partida (Módulo 03):\n")
print(resumen_qc_antes)

# La corrida de demostración es la más pequeña del estudio (WT_rep2,
# SRR5314137, 0.18 gigabases), igual que en el Módulo 03.
corrida_demo <- tabla_muestras %>%
  slice_min(bases_gb, n = 1) %>%
  pull(corrida)

cat("Corrida de demostración:", corrida_demo, "\n")

# --------------------------------------------------------------------------
# 4. Archivos FASTQ de entrada (igual que en el Módulo 03)
# --------------------------------------------------------------------------
# Los archivos FASTQ de la corrida de demostración están en el ENA
# (European Nucleotide Archive), que refleja el SRA y permite
# descargarlos directo por HTTPS, sin el SRA Toolkit. Son dos archivos
# (R1 y R2, datos paired-end). Después de descargarlos se comprueba su
# integridad con el valor MD5 publicado en el Módulo 02.

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

r1_crudo <- file.path("fastq", paste0(corrida_demo, "_1.fastq.gz"))
r2_crudo <- file.path("fastq", paste0(corrida_demo, "_2.fastq.gz"))

# --------------------------------------------------------------------------
# 5. Localizar fastp y Trimmomatic
# --------------------------------------------------------------------------
# fastp se busca primero en el PATH (instalación con conda) y después
# como un binario llamado fastp en esta misma carpeta (por ejemplo, si
# se compiló desde el código fuente). Si no aparece en ningún lado, el
# módulo sigue con Trimmomatic y deja escrito el comando de fastp.
# Trimmomatic se busca como programa instalado, como la ruta en la
# variable de entorno TRIMMOMATIC_JAR, y si no está, el módulo descarga
# el ZIP oficial de la versión 0.39 y lo descomprime en trimmomatic/.

buscar_fastp <- function() {
  ruta_path <- Sys.which("fastp")
  if (nzchar(ruta_path)) {
    return(unname(ruta_path))
  }
  ruta_local <- file.path(getwd(), "fastp")
  if (file.exists(ruta_local)) {
    return(ruta_local)
  }
  return("")
}

ruta_fastp <- buscar_fastp()
if (nzchar(ruta_fastp)) {
  cat("fastp encontrado en:", ruta_fastp, "\n")
} else {
  cat("fastp no está instalado en esta computadora.\n")
  cat("Para instalarlo: conda install -c bioconda fastp\n")
  cat("El módulo sigue con Trimmomatic; el comando de fastp queda en 04_limpiar_lecturas.sh\n")
}

buscar_trimmomatic <- function() {
  ruta_path <- Sys.which("trimmomatic")
  if (nzchar(ruta_path)) {
    return(list(tipo = "programa", ruta = unname(ruta_path)))
  }
  ruta_env <- Sys.getenv("TRIMMOMATIC_JAR", unset = "")
  if (nzchar(ruta_env) && file.exists(ruta_env)) {
    return(list(tipo = "jar", ruta = ruta_env))
  }
  jar_local <- file.path(
    "trimmomatic", "Trimmomatic-0.39", "trimmomatic-0.39.jar"
  )
  if (file.exists(jar_local)) {
    return(list(tipo = "jar", ruta = jar_local))
  }
  # Descargar el ZIP oficial de Trimmomatic 0.39 y descomprimirlo.
  if (!nzchar(Sys.which("java"))) {
    stop("Trimmomatic necesita Java y no se encontró java en el PATH. ",
         "Instala Java (por ejemplo: sudo apt install default-jre) y ",
         "vuelve a correr el módulo.")
  }
  dir.create("trimmomatic", showWarnings = FALSE)
  zip_ruta <- file.path("trimmomatic", "Trimmomatic-0.39.zip")
  if (!file.exists(zip_ruta)) {
    cat("Descargando Trimmomatic 0.39...\n")
    download.file(
      "https://github.com/usadellab/Trimmomatic/releases/download/v0.39/Trimmomatic-0.39.zip",
      zip_ruta, mode = "wb", quiet = TRUE
    )
  }
  utils::unzip(zip_ruta, exdir = "trimmomatic")
  if (!file.exists(jar_local)) {
    stop("No se pudo descomprimir Trimmomatic 0.39.")
  }
  cat("Trimmomatic 0.39 descargado en:", jar_local, "\n")
  return(list(tipo = "jar", ruta = jar_local))
}

trimmomatic <- buscar_trimmomatic()
cat("Trimmomatic:", trimmomatic$tipo, "en", trimmomatic$ruta, "\n")

if (!nzchar(ruta_fastp) && trimmomatic$tipo == "") {
  stop("No se encontró ni fastp ni Trimmomatic. Instala uno de los dos ",
       "(ver la cabecera del script) y vuelve a correr el módulo.")
}

# --------------------------------------------------------------------------
# 6. Guion de limpieza para las 6 corridas del proyecto
# --------------------------------------------------------------------------
# El guion 04_limpiar_lecturas.sh tiene los comandos exactos de fastp
# y de Trimmomatic para las 6 corridas del proyecto, para correrlo tal
# cual:
#   bash 04_limpiar_lecturas.sh
# Los comandos de fastp usan el programa instalado en el PATH; los de
# Trimmomatic usan el archivo JAR descargado por este módulo. Cada
# corrida escribe sus archivos limpios en limpieza/fastp/ y en
# limpieza/trimmomatic/, con su reporte JSON y HTML de fastp.

dir.create(file.path("limpieza", "fastp"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("limpieza", "trimmomatic"), showWarnings = FALSE, recursive = TRUE)

jar_trimmomatic_sh <- if (trimmomatic$tipo == "jar") {
  trimmomatic$ruta
} else {
  trimmomatic$ruta
}
adaptadores_sh <- file.path(
  "trimmomatic", "Trimmomatic-0.39", "adapters", "TruSeq3-PE.fa"
)

lineas_sh <- c(
  "#!/bin/bash",
  "# Limpieza de las 6 corridas de RNA-seq de Chlorella vulgaris (PRJNA378177)",
  "# con fastp y Trimmomatic. Generado por el Módulo 04 de la serie RNA-seq.",
  "# Uso: bash 04_limpiar_lecturas.sh",
  "# Requisitos: fastp instalado (conda install -c bioconda fastp) y Java",
  "# para Trimmomatic (el JAR se descarga con el script de R del módulo).",
  "",
  "mkdir -p limpieza/fastp limpieza/trimmomatic",
  ""
)

for (corrida in tabla_muestras$corrida) {
  lineas_sh <- c(
    lineas_sh,
    paste0(
      "fastp --in1 fastq/", corrida, "_1.fastq.gz ",
      "--in2 fastq/", corrida, "_2.fastq.gz ",
      "--out1 limpieza/fastp/", corrida, "_1.clean.fq.gz ",
      "--out2 limpieza/fastp/", corrida, "_2.clean.fq.gz ",
      "--qualified_quality_phred 20 --unqualified_percent_limit 40 ",
      "--cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 ",
      "--length_required 36 --thread 4 ",
      "--json limpieza/fastp/", corrida, ".fastp.json ",
      "--html limpieza/fastp/", corrida, ".fastp.html"
    ),
    paste0(
      "java -jar ", jar_trimmomatic_sh, " PE -threads 4 -phred33 ",
      "fastq/", corrida, "_1.fastq.gz fastq/", corrida, "_2.fastq.gz ",
      "limpieza/trimmomatic/", corrida, "_1.paired.fq.gz ",
      "limpieza/trimmomatic/", corrida, "_1.unpaired.fq.gz ",
      "limpieza/trimmomatic/", corrida, "_2.paired.fq.gz ",
      "limpieza/trimmomatic/", corrida, "_2.unpaired.fq.gz ",
      "ILLUMINACLIP:", adaptadores_sh, ":2:30:10 ",
      "SLIDINGWINDOW:4:20 MINLEN:36"
    ),
    ""
  )
}

writeLines(lineas_sh, "04_limpiar_lecturas.sh")
cat("Guion escrito: 04_limpiar_lecturas.sh\n")

# --------------------------------------------------------------------------
# 7. Correr fastp sobre la corrida de demostración
# --------------------------------------------------------------------------
# Si fastp está disponible, se corre sobre los dos archivos de la
# corrida de demostración y se guardan los FASTQ limpios y el reporte
# JSON en limpieza/fastp/. Si los archivos ya existen de una corrida
# anterior, se reutilizan y no se vuelve a correr.

fastp_json <- file.path(
  "limpieza", "fastp", paste0(corrida_demo, ".fastp.json")
)
fastp_r1 <- file.path(
  "limpieza", "fastp", paste0(corrida_demo, "_1.clean.fq.gz")
)
fastp_r2 <- file.path(
  "limpieza", "fastp", paste0(corrida_demo, "_2.clean.fq.gz")
)

correr_fastp <- nzchar(ruta_fastp) &&
  !(file.exists(fastp_r1) && file.exists(fastp_r2) && file.exists(fastp_json))

if (correr_fastp) {
  cat("Corriendo fastp sobre la corrida de demostración...\n")
  estado_fastp <- system2(
    ruta_fastp,
    args = c(
      "--in1", r1_crudo, "--in2", r2_crudo,
      "--out1", fastp_r1, "--out2", fastp_r2,
      "--qualified_quality_phred", "20",
      "--unqualified_percent_limit", "40",
      "--cut_front", "--cut_tail",
      "--cut_window_size", "4", "--cut_mean_quality", "20",
      "--length_required", "36",
      "--thread", "4",
      "--json", fastp_json,
      "--html", file.path(
        "limpieza", "fastp", paste0(corrida_demo, ".fastp.html")
      ),
      "--report_title", paste0("fastp ", corrida_demo)
    ),
    stdout = TRUE, stderr = TRUE
  )
  cat(paste(estado_fastp, collapse = "\n"), "\n")
} else if (nzchar(ruta_fastp)) {
  cat("Archivos limpios de fastp ya existentes; se reutilizan.\n")
}

hay_fastp <- file.exists(fastp_json) && file.exists(fastp_r1) && file.exists(fastp_r2)

# --------------------------------------------------------------------------
# 8. Correr Trimmomatic sobre la corrida de demostración
# --------------------------------------------------------------------------
# Trimmomatic escribe 4 archivos: los pares que sobreviven juntos
# (paired) y las lecturas que quedan sin su pareja (unpaired). Su
# resumen se imprime en la salida estándar y se guarda en el archivo
# de registro limpieza/trimmomatic/SRR5314137.trimmomatic.log.

trim_r1_p <- file.path(
  "limpieza", "trimmomatic", paste0(corrida_demo, "_1.paired.fq.gz")
)
trim_r1_u <- file.path(
  "limpieza", "trimmomatic", paste0(corrida_demo, "_1.unpaired.fq.gz")
)
trim_r2_p <- file.path(
  "limpieza", "trimmomatic", paste0(corrida_demo, "_2.paired.fq.gz")
)
trim_r2_u <- file.path(
  "limpieza", "trimmomatic", paste0(corrida_demo, "_2.unpaired.fq.gz")
)
trim_log <- file.path(
  "limpieza", "trimmomatic", paste0(corrida_demo, ".trimmomatic.log")
)
adaptadores <- file.path(
  "trimmomatic", "Trimmomatic-0.39", "adapters", "TruSeq3-PE.fa"
)

correr_trim <- !(file.exists(trim_r1_p) && file.exists(trim_r2_p) &&
                   file.exists(trim_log))

if (correr_trim) {
  cat("Corriendo Trimmomatic sobre la corrida de demostración...\n")
  if (trimmomatic$tipo == "jar") {
    salida_trim <- system2(
      "java",
      args = c(
        "-jar", trimmomatic$ruta, "PE", "-threads", "4", "-phred33",
        r1_crudo, r2_crudo,
        trim_r1_p, trim_r1_u, trim_r2_p, trim_r2_u,
        paste0("ILLUMINACLIP:", adaptadores, ":2:30:10"),
        "SLIDINGWINDOW:4:20", "MINLEN:36"
      ),
      stdout = TRUE, stderr = TRUE
    )
  } else {
    salida_trim <- system2(
      trimmomatic$ruta,
      args = c(
        "PE", "-threads", "4", "-phred33",
        r1_crudo, r2_crudo,
        trim_r1_p, trim_r1_u, trim_r2_p, trim_r2_u,
        paste0("ILLUMINACLIP:", adaptadores, ":2:30:10"),
        "SLIDINGWINDOW:4:20", "MINLEN:36"
      ),
      stdout = TRUE, stderr = TRUE
    )
  }
  writeLines(salida_trim, trim_log)
  cat(paste(salida_trim, collapse = "\n"), "\n")
} else {
  cat("Archivos limpios de Trimmomatic ya existentes; se reutilizan.\n")
  salida_trim <- readLines(trim_log)
}

hay_trim <- file.exists(trim_r1_p) && file.exists(trim_r2_p)

# Resumen de Trimmomatic: la última línea del registro trae los
# conteos. Ejemplo: "Input Read Pairs: 1223634 Both Surviving:
# 1205000 (98.48%) Forward Only Surviving: 12000 (0.98%) Reverse Only
# Surviving: 5000 (0.41%) Dropped: 1634 (0.13%)".
extraer_conteo <- function(texto, etiqueta) {
  coincidencia <- regmatches(
    texto,
    regexec(paste0(etiqueta, ": ([0-9]+)"), texto)
  )[[1]]
  if (length(coincidencia) < 2) {
    return(NA_integer_)
  }
  as.integer(coincidencia[2])
}

linea_resumen_trim <- salida_trim[
  grepl("Input Read Pairs", salida_trim)
]
if (length(linea_resumen_trim) == 0) {
  linea_resumen_trim <- salida_trim[length(salida_trim)]
}
linea_resumen_trim <- linea_resumen_trim[length(linea_resumen_trim)]

trim_pares_entrada <- extraer_conteo(linea_resumen_trim, "Input Read Pairs")
trim_ambos <- extraer_conteo(linea_resumen_trim, "Both Surviving")
trim_solo_r1 <- extraer_conteo(linea_resumen_trim, "Forward Only Surviving")
trim_solo_r2 <- extraer_conteo(linea_resumen_trim, "Reverse Only Surviving")
trim_descartados <- extraer_conteo(linea_resumen_trim, "Dropped")

cat("Resumen de Trimmomatic:\n")
cat(linea_resumen_trim, "\n")

# --------------------------------------------------------------------------
# 9. Qué midió fastp: leer su reporte JSON
# --------------------------------------------------------------------------
# El reporte JSON de fastp trae, para el antes y el después: el número
# de lecturas, el número de bases, la tasa de bases Q20 y Q30, el
# contenido GC y la longitud media de cada extremo. También trae
# cuántas lecturas se descartaron por baja calidad, por tener demasiadas
# N y por quedar demasiado cortas, y cuántas lecturas y bases se
# recortaron por adaptador. Las curvas de calidad por posición (la
# calidad media en cada posición de la lectura) están en
# read1_before_filtering y read1_after_filtering, y lo mismo para R2.

if (hay_fastp) {
  reporte_fastp <- fromJSON(fastp_json, simplifyVector = TRUE)

  antes <- reporte_fastp$summary$before_filtering
  despues <- reporte_fastp$summary$after_filtering
  filtrado <- reporte_fastp$filtering_result
  adaptador <- reporte_fastp$adapter_cutting

  cat("fastp, antes de filtrar: ", antes$total_reads, " lecturas, ",
      antes$total_bases, " bases, Q30 = ",
      round(100 * antes$q30_rate, 2), " %.\n", sep = "")
  cat("fastp, después de filtrar: ", despues$total_reads, " lecturas, ",
      despues$total_bases, " bases, Q30 = ",
      round(100 * despues$q30_rate, 2), " %.\n", sep = "")
  cat("Lecturas descartadas por baja calidad: ", filtrado$low_quality_reads,
      "; por demasiadas N: ", filtrado$too_many_N_reads,
      "; por demasiado cortas: ", filtrado$too_short_reads, ".\n", sep = "")
  cat("Lecturas con adaptador recortado: ", adaptador$adapter_trimmed_reads,
      " (", adaptador$adapter_trimmed_bases, " bases).\n", sep = "")

  # Curvas de calidad media por posición, antes y después, para R1 y R2.
  # Las curvas se guardan primero en variables (curva_antes y
  # curva_despues): dentro de tibble() cada columna nueva queda
  # disponible para las siguientes, así que si la curva se leyera
  # dentro del tibble, el nombre del extremo ya sería la columna y no
  # la variable, y la búsqueda del reporte fallaría.
  calidad_fastp <- map_dfr(
    c("read1", "read2"),
    function(extremo) {
      etiqueta <- if (extremo == "read1") "R1" else "R2"
      curva_antes <- reporte_fastp[[paste0(extremo, "_before_filtering")]]$quality_curves$mean
      curva_despues <- reporte_fastp[[paste0(extremo, "_after_filtering")]]$quality_curves$mean
      bind_rows(
        tibble(
          extremo = etiqueta,
          etapa = "Antes",
          posicion = seq_along(curva_antes),
          calidad = curva_antes
        ),
        tibble(
          extremo = etiqueta,
          etapa = "Después",
          posicion = seq_along(curva_despues),
          calidad = curva_despues
        )
      )
    }
  )
} else {
  cat("Sin reporte de fastp: las gráficas 1 y 4 usan solo Trimmomatic y R.\n")
  calidad_fastp <- NULL
}

# --------------------------------------------------------------------------
# 10. Verificación independiente en R con ShortRead
# --------------------------------------------------------------------------
# fastp y Trimmomatic ya midieron el antes y el después, pero en este
# módulo se comprueba con las lecturas mismas, igual que en el
# Módulo 03: se leen las lecturas crudas y las limpias de fastp con
# ShortRead, se toma una submuestra aleatoria de 50,000 lecturas de
# cada archivo, con semilla fija, y se calculan la longitud, la
# calidad media por lectura y el porcentaje de bases con Q30 o más.
# Si la limpieza funcionó, la longitud media baja un poco (se recortó
# la cola), la calidad media sube y el % Q30 sube, y ninguna lectura
# limpia mide menos de 36 nt.

set.seed(20261112)
n_muestra <- 50000

resumir_fastq <- function(archivo, etiqueta) {
  fq <- readFastq(dirname(archivo), pattern = basename(archivo))
  n_total <- length(fq)
  if (n_total == 0) {
    stop("El archivo ", archivo, " no tiene lecturas.")
  }
  indices <- sample(n_total, size = min(n_muestra, n_total))
  fq_muestra <- fq[indices]
  longitudes <- width(sread(fq_muestra))
  matriz_calidad <- as(quality(fq_muestra), "matrix")
  tibble(
    archivo = etiqueta,
    lecturas_archivo = n_total,
    n_muestra = length(longitudes),
    longitud_media_nt = mean(longitudes),
    longitud_minima_nt = min(longitudes),
    longitud_maxima_nt = max(longitudes),
    calidad_media_phred = sum(matriz_calidad, na.rm = TRUE) / sum(longitudes),
    pct_q30 = 100 * sum(matriz_calidad >= 30, na.rm = TRUE) / sum(longitudes)
  )
}

# Longitud y calidad media de cada lectura de la submuestra, para los
# gráficos 2 y 3. Se leen los mismos archivos que en resumir_fastq.
muestra_fastq <- function(archivo, etiqueta) {
  fq <- readFastq(dirname(archivo), pattern = basename(archivo))
  n_total <- length(fq)
  indices <- sample(n_total, size = min(n_muestra, n_total))
  fq_muestra <- fq[indices]
  matriz_calidad <- as(quality(fq_muestra), "matrix")
  tibble(
    archivo = etiqueta,
    longitud_nt = width(sread(fq_muestra)),
    calidad_media = rowMeans(matriz_calidad, na.rm = TRUE)
  )
}

# La verificación usa los archivos limpios de fastp cuando existen; si
# fastp no está instalado, usa los de Trimmomatic (los pares).
if (hay_fastp) {
  limpio_r1 <- fastp_r1
  limpio_r2 <- fastp_r2
  herramienta_verificacion <- "fastp"
} else {
  limpio_r1 <- trim_r1_p
  limpio_r2 <- trim_r2_p
  herramienta_verificacion <- "Trimmomatic"
}
cat("Verificación en R con los archivos limpios de:", herramienta_verificacion, "\n")

resumen_r <- bind_rows(
  resumir_fastq(r1_crudo, "R1 crudo"),
  resumir_fastq(r2_crudo, "R2 crudo"),
  resumir_fastq(limpio_r1, "R1 limpio"),
  resumir_fastq(limpio_r2, "R2 limpio")
)
print(resumen_r)

muestra_r <- bind_rows(
  muestra_fastq(r1_crudo, "R1 crudo"),
  muestra_fastq(r2_crudo, "R2 crudo"),
  muestra_fastq(limpio_r1, "R1 limpio"),
  muestra_fastq(limpio_r2, "R2 limpio")
) %>%
  mutate(
    extremo = if_else(grepl("^R1", archivo), "R1", "R2"),
    etapa = if_else(grepl("crudo", archivo), "Antes", "Después")
  )

# --------------------------------------------------------------------------
# 11. Tabla resumen de la limpieza
# --------------------------------------------------------------------------
# La tabla 04_resumen_limpieza.csv junta lo que midió cada
# herramienta: fastp (del JSON) y Trimmomatic (de su registro), y la
# verificación en R con ShortRead. Trimmomatic no reporta bases ni
# calidad, así que esas celdas quedan vacías y se toman de la
# verificación en R. La tabla es el registro de la limpieza del
# proyecto y la usarán los módulos siguientes.

filas_resumen <- list()

if (hay_fastp) {
  filas_resumen <- c(
    filas_resumen,
    list(tibble(
      herramienta = "fastp",
      etapa = "Antes",
      lecturas = antes$total_reads,
      bases = antes$total_bases,
      longitud_media_nt = (antes$read1_mean_length + antes$read2_mean_length) / 2,
      pct_q30 = 100 * antes$q30_rate,
      pct_gc = 100 * antes$gc_content
    )),
    list(tibble(
      herramienta = "fastp",
      etapa = "Después",
      lecturas = despues$total_reads,
      bases = despues$total_bases,
      longitud_media_nt = (despues$read1_mean_length + despues$read2_mean_length) / 2,
      pct_q30 = 100 * despues$q30_rate,
      pct_gc = 100 * despues$gc_content
    ))
  )
}

if (hay_trim) {
  filas_resumen <- c(
    filas_resumen,
    list(tibble(
      herramienta = "Trimmomatic",
      etapa = "Antes",
      lecturas = 2 * trim_pares_entrada,
      bases = NA_real_,
      longitud_media_nt = NA_real_,
      pct_q30 = NA_real_,
      pct_gc = NA_real_
    )),
    list(tibble(
      herramienta = "Trimmomatic",
      etapa = "Después",
      lecturas = 2 * trim_ambos + trim_solo_r1 + trim_solo_r2,
      bases = NA_real_,
      longitud_media_nt = NA_real_,
      pct_q30 = NA_real_,
      pct_gc = NA_real_
    ))
  )
}

filas_resumen <- c(
  filas_resumen,
  list(resumen_r %>%
    mutate(
      herramienta = paste0("R/ShortRead (", herramienta_verificacion, ")"),
      etapa = if_else(grepl("crudo", archivo), "Antes", "Después"),
      lecturas = lecturas_archivo,
      bases = NA_real_,
      pct_q30 = pct_q30,
      pct_gc = NA_real_
    ) %>%
    # La columna archivo distingue R1 de R2 en las filas de R; las
    # filas de fastp y Trimmomatic no la tienen y quedan con NA.
    select(herramienta, archivo, etapa, lecturas, bases, longitud_media_nt, pct_q30, pct_gc))
)

# fastp reporta la longitud media como número entero (75 nt), por eso
# en sus filas la longitud media es 75 y en las de R es 75.5 nt: es el
# mismo valor, uno redondeado y el otro medido lectura por lectura.

resumen_limpieza <- bind_rows(filas_resumen)
print(resumen_limpieza)
write_csv(resumen_limpieza, "04_resumen_limpieza.csv")

# Lecturas y bases retenidas según fastp, para el gráfico 4.
if (hay_fastp) {
  retenidas <- tibble(
    metrica = rep(c("Lecturas", "Bases"), each = 2),
    etapa = rep(c("Antes", "Después"), times = 2),
    cantidad = c(
      antes$total_reads, despues$total_reads,
      antes$total_bases, despues$total_bases
    )
  )
  cat("Lecturas retenidas por fastp: ",
      round(100 * despues$total_reads / antes$total_reads, 2), " %.\n",
      sep = "")
  cat("Bases retenidas por fastp: ",
      round(100 * despues$total_bases / antes$total_bases, 2), " %.\n",
      sep = "")
}
if (hay_trim) {
  cat("Pares retenidos completos por Trimmomatic: ",
      round(100 * trim_ambos / trim_pares_entrada, 2), " %.\n", sep = "")
  cat("Pares descartados por Trimmomatic: ",
      round(100 * trim_descartados / trim_pares_entrada, 2), " %.\n", sep = "")
}

# --------------------------------------------------------------------------
# 12. Gráfico 1: calidad media por posición, antes y después (fastp)
# --------------------------------------------------------------------------
# La curva es la calidad media Phred que reporta fastp en cada
# posición de la lectura, antes y después de limpiar. Las líneas
# punteadas marcan Q20 y Q30.

if (!is.null(calidad_fastp)) {
  grafico_calidad_antes_despues <- ggplot(calidad_fastp, aes(
    x = posicion, y = calidad, color = etapa, linetype = extremo
  )) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 30, linetype = "dashed", color = "#59A14F") +
    geom_hline(yintercept = 20, linetype = "dashed", color = "#F28E2B") +
    annotate("text", x = 74, y = 31, label = "Q30", color = "#59A14F", size = 3.5) +
    annotate("text", x = 74, y = 21, label = "Q20", color = "#F28E2B", size = 3.5) +
    scale_color_manual(values = c("Antes" = "#E15759", "Después" = "#4E79A7")) +
    labs(
      title = "Calidad por base antes y después de la limpieza (SRR5314137)",
      subtitle = "fastp: curvas de calidad media por posición, extremos R1 y R2",
      x = "Posición en la lectura (nucleótidos, nt)",
      y = "Calidad media Phred (adimensional)",
      color = "Etapa",
      linetype = "Extremo"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  print(grafico_calidad_antes_despues)

  # Interpretación:
  # El eje horizontal es la posición dentro de la lectura, de 1 a 76 nt
  # (nucleótidos); el eje vertical es la calidad media Phred, un número
  # adimensional: Q30 significa 1 error por cada 1,000 bases y Q20, 1
  # error por cada 100 bases. La curva roja es la lectura cruda y la
  # azul es la lectura limpia; la línea continua es R1 y la punteada es
  # R2. En la lectura cruda la calidad cae hacia el final, como se vio
  # en el Módulo 03. Después de limpiar, la curva azul queda por encima
  # de la roja en casi todas las posiciones: los extremos de baja
  # calidad se recortaron y la calidad media sube. La caída de la
  # última posición (76 nt) no desaparece del todo: las lecturas que
  # conservan los 76 nt siguen teniendo su base final más baja, pero
  # la caída es menor que antes (en R2 pasa de 30.6 a 31.7 de Phred y
  # en R1 de 34.0 a 34.4). Si la curva azul quedara por debajo de Q30
  # en alguna zona, habría que revisar los parámetros de recorte.

  ggsave(
    "04_calidad_por_base_antes_despues.png",
    plot = grafico_calidad_antes_despues,
    width = 9, height = 5.5, dpi = 300, bg = "white"
  )
}

# --------------------------------------------------------------------------
# 13. Gráfico 2: longitud de las lecturas, antes y después (R/ShortRead)
# --------------------------------------------------------------------------
# El gráfico de barras cuenta cuántas lecturas hay de cada longitud en
# la submuestra de 50,000 lecturas, antes y después de limpiar. Se
# calcula en R con ShortRead, igual que en el Módulo 03. Se escribe
# dplyr::count con su paquete porque al cargar ShortRead el paquete
# matrixStats pone su propia función count por delante y la de dplyr
# dejaría de encontrarse.

longitud_resumen <- muestra_r %>%
  dplyr::count(extremo, etapa, longitud_nt, name = "lecturas")

grafico_longitud_antes_despues <- ggplot(longitud_resumen, aes(
  x = longitud_nt, y = lecturas, fill = etapa
)) +
  geom_col(position = "dodge") +
  facet_wrap(~extremo, ncol = 1) +
  scale_fill_manual(values = c("Antes" = "#E15759", "Después" = "#4E79A7")) +
  labs(
    title = "Longitud de las lecturas antes y después de la limpieza (SRR5314137)",
    subtitle = "Submuestra de 50,000 lecturas por archivo, verificación en R con ShortRead",
    x = "Longitud de lectura (nucleótidos, nt)",
    y = "Número de lecturas",
    fill = "Etapa"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_longitud_antes_despues)

# Interpretación:
# El eje horizontal es la longitud de la lectura en nucleótidos (nt)
# y el eje vertical cuenta cuántas lecturas tienen esa longitud, en la
# submuestra de 50,000 lecturas por archivo. Antes de limpiar, la barra
# más alta está en 76 nt, la longitud modal, con una cola de lecturas
# más cortas a la izquierda (ver Módulo 03, gráfico 5). Después de
# limpiar, la barra de 76 nt sigue siendo la más alta pero baja un
# poco, y aparecen barras azules en longitudes menores: una parte de
# las lecturas perdió bases del final al recortar la cola de baja
# calidad. La longitud media baja poco (en R1 de 75.55 a 75.45 nt y en
# R2 de 75.59 a 75.08 nt), porque la mayoría de las lecturas conserva
# sus 76 nt. Lo que no debe aparecer es ninguna barra por debajo de 36
# nt: las lecturas que quedaron con menos de 36 nt se descartaron con
# MINLEN y --length_required, y por eso la distribución limpia empieza
# en 36 nt. Esa es la comprobación visual de que el filtro de longitud
# funcionó.

ggsave(
  "04_longitud_antes_despues.png",
  plot = grafico_longitud_antes_despues,
  width = 9, height = 6.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 14. Gráfico 3: calidad media por lectura, antes y después (R/ShortRead)
# --------------------------------------------------------------------------
# Cada lectura tiene una calidad media (su Phred promedio). El
# histograma muestra cuántas lecturas caen en cada valor de calidad,
# antes y después de limpiar.

grafico_calidad_lectura_antes_despues <- ggplot(muestra_r, aes(
  x = calidad_media, fill = etapa
)) +
  geom_histogram(binwidth = 1, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 30, linetype = "dashed", color = "#59A14F") +
  geom_vline(xintercept = 20, linetype = "dashed", color = "#F28E2B") +
  facet_wrap(~extremo, ncol = 1) +
  scale_fill_manual(values = c("Antes" = "#E15759", "Después" = "#4E79A7")) +
  labs(
    title = "Calidad media por lectura antes y después de la limpieza (SRR5314137)",
    subtitle = "Submuestra de 50,000 lecturas por archivo, verificación en R con ShortRead",
    x = "Calidad media por lectura (Phred, adimensional)",
    y = "Número de lecturas",
    fill = "Etapa"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_calidad_lectura_antes_despues)

# Interpretación:
# El eje horizontal es la calidad media Phred de cada lectura (número
# adimensional; Q30 es 1 error por cada 1,000 bases) y el eje vertical
# cuenta cuántas lecturas tienen esa calidad, en la submuestra de
# 50,000 lecturas por archivo. Antes de limpiar, el histograma rojo
# forma un pico por encima de Q30 con una cola hacia la izquierda:
# son las lecturas que traen bases malas al final. Después de limpiar,
# el histograma azul queda casi igual que el rojo, con un
# desplazamiento pequeño hacia la derecha: la cola de calidades bajas
# se acorta porque esas bases se recortaron y las lecturas malas se
# descartaron, y el pico se mueve un poco hacia calidades más altas
# (en R2 el % Q30 sube de 84.8 a 86.1 y en R1 de 92.7 a 93.1). La
# línea verde marca Q30 y la naranja Q20. Si después de limpiar el
# pico quedara por debajo de Q30, la corrida tendría un problema de
# calidad que el recorte no alcanza a corregir.

ggsave(
  "04_calidad_por_lectura_antes_despues.png",
  plot = grafico_calidad_lectura_antes_despues,
  width = 9, height = 6.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 15. Gráfico 4: lecturas y bases retenidas (fastp)
# --------------------------------------------------------------------------
# Las barras comparan cuántas lecturas y cuántas bases había antes de
# limpiar y cuántas quedaron después, según el reporte JSON de fastp.

if (hay_fastp) {
  grafico_retenidas <- ggplot(retenidas, aes(
    x = etapa, y = cantidad, fill = etapa
  )) +
    geom_col() +
    geom_text(
      aes(label = format(cantidad, big.mark = ",", scientific = FALSE)),
      vjust = -0.5, size = 3.5
    ) +
    facet_wrap(~metrica, scales = "free_y") +
    scale_fill_manual(values = c("Antes" = "#E15759", "Después" = "#4E79A7")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Lecturas y bases retenidas tras la limpieza (SRR5314137)",
      subtitle = "fastp: totales antes y después de filtrar, corrida completa",
      x = NULL,
      y = "Cantidad",
      fill = "Etapa"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  print(grafico_retenidas)

  # Interpretación:
  # Cada panel compara el antes y el después de la limpieza con los
  # totales de la corrida completa (1,223,634 lecturas por archivo,
  # 2,447,268 en total), según el reporte JSON de fastp. El panel de
  # la izquierda cuenta lecturas y el de la derecha cuenta bases; el
  # eje vertical de cada panel tiene su propia escala porque hay muchas
  # más bases que lecturas (cada lectura mide unos 76 nt). Las dos
  # barras de cada panel quedan casi iguales: la limpieza conserva
  # casi todas las lecturas y casi todas las bases, y solo descarta
  # las lecturas malas y recorta las colas de baja calidad. Si la
  # barra de después fuera mucho más baja que la de antes (por
  # ejemplo, menos del 90 % de las lecturas), los parámetros de
  # limpieza serían demasiado estrictos o la corrida tendría un
  # problema de calidad.

  ggsave(
    "04_lecturas_retenidas.png",
    plot = grafico_retenidas,
    width = 9, height = 5.5, dpi = 300, bg = "white"
  )
}

# --------------------------------------------------------------------------
# 16. Veredicto: ¿la limpieza funcionó?
# --------------------------------------------------------------------------
# Se decide con tres números de la tabla resumen: el porcentaje de
# lecturas retenidas, el porcentaje de bases con Q30 o más después de
# limpiar (que debe ser mayor o igual que antes), y la longitud mínima
# de las lecturas limpias (que debe ser 36 nt, el filtro aplicado).
# La regla que se usa en este proyecto es: al menos 90 % de las
# lecturas retenidas, Q30 después mayor o igual que Q30 antes, y
# ninguna lectura limpia por debajo de 36 nt. Si se cumplen las tres,
# las lecturas limpias pasan al Módulo 05 (ensamblaje de novo con
# Trinity).

veredicto_limpieza <- resumen_r %>%
  mutate(etapa = if_else(grepl("crudo", archivo), "Antes", "Después")) %>%
  group_by(etapa) %>%
  summarise(
    pct_q30 = mean(pct_q30),
    longitud_minima_nt = min(longitud_minima_nt),
    .groups = "drop"
  )

q30_antes <- veredicto_limpieza$pct_q30[veredicto_limpieza$etapa == "Antes"]
q30_despues <- veredicto_limpieza$pct_q30[veredicto_limpieza$etapa == "Después"]
min_despues <- veredicto_limpieza$longitud_minima_nt[veredicto_limpieza$etapa == "Después"]

retenidas_pct <- if (hay_fastp) {
  100 * despues$total_reads / antes$total_reads
} else if (hay_trim) {
  100 * (2 * trim_ambos + trim_solo_r1 + trim_solo_r2) / (2 * trim_pares_entrada)
} else {
  NA_real_
}

cat("Lecturas retenidas: ", round(retenidas_pct, 2), " %.\n", sep = "")
cat("Q30 antes: ", round(q30_antes, 2), " %; Q30 después: ",
    round(q30_despues, 2), " %.\n", sep = "")
cat("Longitud mínima después de limpiar: ", min_despues, " nt.\n", sep = "")

pasa_limpieza <- !is.na(retenidas_pct) && retenidas_pct >= 90 &&
  q30_despues >= q30_antes && min_despues >= 36

if (pasa_limpieza) {
  cat("Veredicto: la limpieza FUNCIONÓ.\n")
  cat("Las lecturas limpias pasan al Módulo 05 (ensamblaje de novo con Trinity).\n")
} else {
  cat("Veredicto: la limpieza NO cumple los criterios del proyecto.\n")
  cat("Revisa los gráficos y la tabla 04_resumen_limpieza.csv antes de seguir.\n")
}

# --------------------------------------------------------------------------
# 17. Resumen del módulo
# --------------------------------------------------------------------------
# Archivos generados (guárdalos: los usaremos en los módulos siguientes):
#   - 04_resumen_limpieza.csv           : tabla con las métricas de la
#     limpieza de fastp, Trimmomatic y la verificación en R.
#   - 04_limpiar_lecturas.sh            : guion con los comandos fastp y
#     Trimmomatic para las 6 corridas del proyecto.
#   - 04_calidad_por_base_antes_despues.png : gráfico de la calidad por
#     posición antes y después (fastp).
#   - 04_longitud_antes_despues.png     : gráfico de la longitud de
#     lectura antes y después (R/ShortRead).
#   - 04_calidad_por_lectura_antes_despues.png : gráfico de la calidad
#     media por lectura antes y después (R/ShortRead).
#   - 04_lecturas_retenidas.png         : gráfico de lecturas y bases
#     retenidas (fastp).
#   - limpieza/fastp/SRR5314137_1.clean.fq.gz y
#     limpieza/fastp/SRR5314137_2.clean.fq.gz : lecturas limpias con
#     fastp, con su reporte JSON y HTML (solo si fastp está instalado).
#   - limpieza/trimmomatic/SRR5314137_1.paired.fq.gz y
#     limpieza/trimmomatic/SRR5314137_2.paired.fq.gz : lecturas limpias
#     con Trimmomatic (pares), más los archivos .unpaired.fq.gz con
#     las lecturas que quedaron sin pareja y el registro .log.
#   - trimmomatic/Trimmomatic-0.39/      : Trimmomatic 0.39 descargado
#     por el módulo (programa y adaptadores).
#
# Qué sigue (Módulo 05):
#   Ensamblaje de novo del transcriptoma con Trinity, a partir de las
#   lecturas limpias de las 6 corridas del proyecto.

# --------------------------------------------------------------------------
# Qué aprendimos
# --------------------------------------------------------------------------
# - La limpieza va después del control de calidad y antes de cualquier
#   análisis: quita los restos de adaptador, recorta los extremos de
#   baja calidad y descarta las lecturas que quedan demasiado cortas.
# - Un adaptador es la secuencia que se añade al fragmento para
#   secuenciarlo; cuando el fragmento es más corto que la lectura, el
#   adaptador queda pegado al final y hay que cortarlo.
# - Trimmomatic trabaja por pasos encadenados: ILLUMINACLIP corta el
#   adaptador (2:30:10), SLIDINGWINDOW:4:20 recorta con una ventana
#   deslizante de 4 bases y calidad media mínima de 20, y MINLEN:36
#   descarta las lecturas de menos de 36 nt.
# - fastp hace lo mismo en un solo programa y escribe un reporte JSON
#   y uno HTML: detecta el adaptador solo, recorta con --cut_front y
#   --cut_tail, descarta por calidad con --qualified_quality_phred 20
#   y --unqualified_percent_limit 40, y por longitud con
#   --length_required 36.
# - En datos paired-end los pares se conservan juntos: Trimmomatic
#   escribe los pares que sobreviven en los archivos paired y las
#   lecturas que quedan sin pareja en los archivos unpaired.
# - La limpieza se comprueba con números: lecturas y bases retenidas,
#   porcentaje de bases con Q30 o más antes y después, y longitud
#   mínima de las lecturas limpias. En la corrida de demostración se
#   conserva casi todo y la calidad sube, así que las lecturas limpias
#   pasan al ensamblaje de novo del Módulo 05.
# - En R, el paquete ShortRead de Bioconductor lee los FASTQ crudos y
#   los limpios y permite verificar la limpieza de forma independiente
#   de las herramientas: longitud, calidad media y % Q30 se calculan
#   desde las lecturas mismas, con una submuestra de 50,000 lecturas
#   por archivo y semilla fija.
