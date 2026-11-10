# ==========================================================================
# Módulo 02. Descargar datos de RNA-seq del SRA (Sequence Read Archive)
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
# Serie "RNA-seq" (avanzado) — Módulo 02
#
# Descripción: En el Módulo 01 buscamos literatura y secuencias. Ahora
#   descargamos datos reales de RNA-seq: las lecturas crudas de un
#   experimento de Chlorella vulgaris guardado en el SRA (Sequence Read
#   Archive, el archivo público de datos de secuenciación del NCBI).
#   Primero buscamos el experimento desde R, luego lo descargamos con
#   el SRA Toolkit (prefetch y fasterq-dump) y dejamos lista la tabla
#   de muestras que usarán todos los módulos siguientes.
#
# Objetivos de aprendizaje:
#   - Entender qué es el SRA y cómo se organiza: BioProject, experimento
#     (SRX) y corrida (SRR).
#   - Buscar experimentos de RNA-seq en el SRA desde R con rentrez.
#   - Instalar y usar el SRA Toolkit: prefetch y fasterq-dump.
#   - Descargar archivos FASTQ y comprobar que la descarga quedó completa.
#   - Construir la tabla de muestras del proyecto (muestra, condición,
#     réplica, corrida) y guardarla en CSV.
#
# Requisitos: R con los paquetes tidyverse, rentrez, httr y jsonlite.
#   Instálalos con: install.packages(c("tidyverse", "rentrez", "httr", "jsonlite"))
#   Nota: en Linux con el locale en POSIX, los acentos de las gráficas
#   salen mal; en ese caso corre el script con LC_ALL=C.utf8 Rscript.
#   En macOS y RStudio el locale ya es UTF-8 y no hace falta.
#   Además, el SRA Toolkit (programa de línea de comandos, fuera de R):
#     - Con conda:  conda install -c bioconda sra-tools
#     - En macOS:   brew install sra-tools
#     - En Ubuntu:  sudo apt install sra-toolkit
#   La primera vez se configura la carpeta de descarga con:
#     vdb-config --interactive
#
# Paquetes de Bioconductor relacionados con este módulo:
#   - SRAdb: consulta los metadatos del SRA desde R, con una base de datos
#     local (SQLite) que se descarga una vez. Sirve para buscar estudios
#     y corridas sin depender de la conexión en cada búsqueda.
#   - ShortRead: lee archivos FASTQ dentro de R (se usará en módulos
#     siguientes para revisar las lecturas descargadas).
#   Se instalan con BiocManager:
#     if (!requireNamespace("BiocManager", quietly = TRUE)) {
#       install.packages("BiocManager")
#     }
#     BiocManager::install(c("SRAdb", "ShortRead"))
#   En este módulo la búsqueda se hace con rentrez (como en el Módulo 01)
#   y la descarga con el SRA Toolkit, porque descargar FASTQ es trabajo
#   de la terminal, no de R.
#
# Organismo de trabajo: Chlorella vulgaris (microalga verde, TaxID 3077).
# ==========================================================================

# --------------------------------------------------------------------------
# 1. Cargar paquetes
# --------------------------------------------------------------------------
library(tidyverse)
library(rentrez)
library(httr)
library(jsonlite)

# Las consultas al NCBI a veces fallan por un corte momentáneo de la red.
# Esta función reintenta la consulta hasta 3 veces antes de rendirse,
# para que un fallo pasajero no detenga el módulo.
con_reintento <- function(expr, intentos = 3) {
  for (i in seq_len(intentos)) {
    resultado <- tryCatch(expr, error = function(e) e)
    if (!inherits(resultado, "error")) {
      return(resultado)
    }
    if (i < intentos) {
      message("La consulta al NCBI falló (intento ", i, " de ", intentos,
              "). Reintentando en 5 segundos...")
      Sys.sleep(5)
    }
  }
  stop(resultado)
}

# --------------------------------------------------------------------------
# 2. ¿Qué es el SRA y qué es el SRA Toolkit?
# --------------------------------------------------------------------------
# El SRA (Sequence Read Archive) es el archivo público del NCBI donde se
# depositan las lecturas crudas de secuenciación de casi todos los
# artículos que publican RNA-seq. Cada depósito se organiza en tres
# niveles, cada uno con su propio código de acceso:
#   - BioProject (PRJNA...): el proyecto completo, por ejemplo un artículo.
#   - Experimento (SRX...): una muestra dentro del proyecto (una condición
#     y una réplica, por ejemplo "WT_rep1").
#   - Corrida o run (SRR...): los archivos de lecturas de ese experimento.
#     Una corrida PAIRED produce dos archivos FASTQ (_1 y _2); una SINGLE,
#     uno solo.
# Dos números describen el tamaño de una corrida:
#   - spots: número de fragmentos secuenciados (en datos paired-end, cada
#     spot produce dos lecturas, una por cada extremo del fragmento).
#   - bases: número total de nucleótidos (nt) secuenciados en la corrida.
# El SRA Toolkit es el conjunto de programas del NCBI para bajar esos
# datos a tu computadora. Los dos que usamos aquí son:
#   - prefetch: descarga la corrida del SRA a la carpeta local (.sra).
#   - fasterq-dump: convierte el archivo .sra en archivos FASTQ (.fastq.gz),
#     que son los que usan FastQC, Trimmomatic, Trinity y Salmon en los
#     módulos siguientes. Es la versión actual de la herramienta; el
#     antiguo fastq-dump ya no se recomienda.

# --------------------------------------------------------------------------
# 3. Buscar en el SRA experimentos de RNA-seq de Chlorella vulgaris
# --------------------------------------------------------------------------
# Buscamos en la base de datos "sra" todos los experimentos de RNA-seq
# del organismo. El número que devuelve la búsqueda cambia con el tiempo,
# porque el SRA crece cada semana: lo que se imprime aquí es el conteo
# real del día en que se corre el script.

busqueda_sra <- con_reintento(entrez_search(
  db = "sra",
  term = "Chlorella vulgaris[Organism] AND RNA-Seq[Strategy]",
  retmax = 20
))

cat("Experimentos de RNA-seq de Chlorella vulgaris en el SRA:",
    busqueda_sra$count, "\n")

# --------------------------------------------------------------------------
# 4. Elegir el estudio de trabajo: PRJNA378177
# --------------------------------------------------------------------------
# De todos los experimentos del SRA elegimos uno para todo el proyecto.
# Los criterios son prácticos, para que el análisis corra en una
# computadora normal:
#   - Pequeño: cada corrida tiene entre 0.18 y 0.79 gigabases, así que el
#     proyecto completo ocupa unas 3.2 gigabases (unas 2.1 GB en FASTQ
#     comprimido), no decenas de gigas.
#   - Diseño completo: 6 muestras, 3 condiciones con 2 réplicas cada una.
#   - Tema del proyecto: mutantes de Chlorella vulgaris que producen más
#     lípidos que la cepa silvestre (WT), obtenidos por mutagénesis con
#     luz UV (UV) y con EMS (metanosulfonato de etilo). Es el mismo tema
#     de los artículos del Módulo 01: acumulación de lípidos.
# El BioProject es PRJNA378177, con título "Understanding lipid metabolism
# in high-lipid-producing Chlorella vulgaris mutants at the genome-wide
# level" (estudio de 2017, secuenciación Illumina MiSeq, paired-end).
#
# Buscamos las 6 muestras del estudio por su BioProject y pedimos el
# registro completo de cada experimento en formato XML. De ese XML
# sacamos, para cada muestra: el nombre (WT_rep1, EMS_rep2, ...), el
# experimento SRX, la corrida SRR, el diseño (PAIRED), el instrumento
# y el tamaño (spots y bases).

busqueda_estudio <- con_reintento(entrez_search(
  db = "sra",
  term = "PRJNA378177",
  retmax = 10
))

if (length(busqueda_estudio$ids) == 0) {
  stop("La búsqueda del BioProject PRJNA378177 no devolvió experimentos.")
}

xml_estudio <- con_reintento(entrez_fetch(
  db = "sra",
  id = busqueda_estudio$ids,
  rettype = "full"
))

# El XML trae un bloque <EXPERIMENT_PACKAGE> por muestra. Separamos los
# bloques y extraemos los campos de cada uno con expresiones regulares.
# El primer trozo que resulta de la separación es la cabecera del XML
# y se descarta.
bloques <- str_split(xml_estudio, "<EXPERIMENT_PACKAGE>")[[1]][-1]

extraer_muestra <- function(bloque) {
  titulo <- str_match(bloque, "<TITLE>([^<]+)</TITLE>")[, 2]
  # El título tiene la forma "GSM2522405: WT_rep1; Chlorella vulgaris; RNA-Seq"
  partes <- str_match(titulo, ": ([A-Za-z]+)_rep([0-9]+)")
  muestra <- paste0(partes[, 2], "_rep", partes[, 3])
  condicion <- partes[, 2]
  replica <- partes[, 3]
  experimento <- str_match(bloque, 'accession="(SRX[0-9]+)"')[, 2]
  corrida <- str_match(
    bloque,
    '<RUN[^>]*accession="(SRR[0-9]+)"[^>]*total_spots="([0-9]+)"[^>]*total_bases="([0-9]+)"'
  )
  disenho <- str_match(bloque, "<LIBRARY_LAYOUT>\\s*<(\\w+)")[, 2]
  instrumento <- str_match(bloque, "<INSTRUMENT_MODEL>([^<]+)</INSTRUMENT_MODEL>")[, 2]
  # Extraer los valores a escalares antes de tibble(): dentro de tibble()
  # cada columna puede referirse a las anteriores, así que una columna
  # llamada "corrida" ocultaría la matriz "corrida" de aquí arriba.
  acceso_corrida <- corrida[, 2]
  n_spots <- as.numeric(corrida[, 3])
  n_bases <- as.numeric(corrida[, 4])
  tibble(
    muestra = muestra,
    condicion = condicion,
    replica = as.integer(replica),
    experimento = experimento,
    corrida = acceso_corrida,
    disenho = disenho,
    instrumento = instrumento,
    spots = n_spots,
    bases = n_bases
  )
}

tabla_muestras <- map_dfr(bloques, extraer_muestra) %>%
  arrange(condicion, replica) %>%
  mutate(
    # bases en gigabases (1 gigabase = 1,000 millones de nucleótidos)
    bases_gb = bases / 1e9,
    # lecturas: en datos paired-end cada spot produce 2 lecturas
    lecturas = if_else(disenho == "PAIRED", spots * 2, spots),
    # longitud media de cada lectura, en nucleótidos (nt):
    # bases totales divididas entre el número de lecturas
    longitud_lectura_nt = bases / lecturas
  )

print(tabla_muestras)

# Guardar la tabla de muestras: es el archivo que ordena todo el
# proyecto y la usarán los módulos 03 a 13.
write_csv(tabla_muestras, "02_muestras_chlorella_sra.csv")

cat("Muestras:", nrow(tabla_muestras), "\n")
cat("Spots totales:", format(sum(tabla_muestras$spots), big.mark = ","), "\n")
cat("Bases totales (gigabases):", round(sum(tabla_muestras$bases_gb), 2), "\n")

# --------------------------------------------------------------------------
# 5. Escribir el guion de descarga con el SRA Toolkit
# --------------------------------------------------------------------------
# La descarga se hace en la terminal, con prefetch y fasterq-dump. Este
# bloque escribe el guion 02_descargar_sra.sh con los comandos exactos
# para las 6 corridas, para correrlo tal cual:
#   bash 02_descargar_sra.sh
# Cada corrida se descarga con prefetch y se convierte a FASTQ con
# fasterq-dump --split-files (separa R1 y R2 en datos paired-end) y
# --gzip (comprime la salida). Los archivos quedan en la carpeta
# fastq/ con los nombres SRR5314136_1.fastq.gz y SRR5314136_2.fastq.gz.
# Se descargan en total unas 2.1 GB, así que conviene conexión estable
# y espacio libre en disco antes de correr el guion completo.

lineas_sh <- c(
  "#!/bin/bash",
  "# Descarga de las 6 corridas de RNA-seq de Chlorella vulgaris (PRJNA378177)",
  "# con el SRA Toolkit. Generado por el Módulo 02 de la serie RNA-seq.",
  "# Uso: bash 02_descargar_sra.sh",
  "",
  "mkdir -p fastq",
  ""
)

for (corrida in tabla_muestras$corrida) {
  lineas_sh <- c(
    lineas_sh,
    paste0("prefetch ", corrida),
    paste0("fasterq-dump --split-files --gzip --outdir fastq ", corrida),
    ""
  )
}

writeLines(lineas_sh, "02_descargar_sra.sh")

# Si el SRA Toolkit está instalado en esta computadora, descargamos
# ahora mismo la corrida más pequeña (SRR5314137, WT_rep2, 0.18
# gigabases, unos 146 MB en FASTQ comprimido) para comprobar que la
# descarga funciona. Si no está instalado, el script sigue adelante:
# la tabla de muestras y las gráficas ya están hechas, y el guion
# 02_descargar_sra.sh queda listo para cuando se instale.
corrida_demo <- tabla_muestras %>%
  slice_min(bases_gb, n = 1) %>%
  pull(corrida)

cat("Corrida de demostración (la más pequeña):", corrida_demo, "\n")

if (nzchar(Sys.which("fasterq-dump"))) {
  dir.create("fastq", showWarnings = FALSE)
  system2("prefetch", args = corrida_demo)
  system2(
    "fasterq-dump",
    args = c("--split-files", "--gzip", "--outdir", "fastq", corrida_demo)
  )
  archivos_fastq <- list.files(
    "fastq",
    pattern = paste0("^", corrida_demo, "_[12]\\.fastq\\.gz$"),
    full.names = TRUE
  )
  if (length(archivos_fastq) == 2) {
    cat("Descarga de demostración completa:\n")
    print(file.info(archivos_fastq)[, c("size")])
  } else {
    cat("La descarga de demostración no produjo los 2 archivos FASTQ esperados.\n")
  }
} else {
  cat("El SRA Toolkit (fasterq-dump) no está instalado en esta computadora.\n")
  cat("Instálalo (ver la cabecera del script) y corre: bash 02_descargar_sra.sh\n")
}

# Ruta alternativa: los mismos archivos FASTQ también están en el ENA
# (European Nucleotide Archive), que refleja el SRA y ofrece descarga
# directa por FTP. Para la corrida de demostración:
#   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR531/007/SRR5314137/SRR5314137_1.fastq.gz
#   wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR531/007/SRR5314137/SRR5314137_2.fastq.gz
# Los valores MD5 de esos dos archivos (para comprobar la descarga con
# md5sum) son:
#   SRR5314137_1.fastq.gz: 79ab761d57b4a04baa4f64ee4859aecb
#   SRR5314137_2.fastq.gz: 77a05c6fd6eb6a91ce3cdf6c4d4d288c

# --------------------------------------------------------------------------
# 6. Gráfico 1: tamaño de cada corrida, en gigabases
# --------------------------------------------------------------------------
# Cada barra es una corrida (SRR) y su altura es el total de bases
# secuenciadas, en gigabases. El color marca la condición: WT (cepa
# silvestre), EMS (mutante por EMS) y UV (mutante por luz UV).

grafico_bases <- ggplot(tabla_muestras, aes(
  x = reorder(corrida, bases_gb),
  y = bases_gb,
  fill = condicion
)) +
  geom_col() +
  geom_text(
    aes(label = paste0(round(bases_gb, 2), " Gb")),
    hjust = -0.1, size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c(
    "WT" = "#4E79A7", "EMS" = "#F28E2B", "UV" = "#59A14F"
  )) +
  labs(
    title = "Tamaño de cada corrida de RNA-seq (PRJNA378177)",
    subtitle = "Chlorella vulgaris: cepa silvestre (WT) y mutantes de alta producción de lípidos (EMS y UV)",
    x = "Corrida (código SRR del SRA)",
    y = "Bases secuenciadas (gigabases, Gb)",
    fill = "Condición"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_bases)

# Interpretación:
# El eje horizontal mide las bases secuenciadas en gigabases (Gb), donde
# 1 Gb son 1,000 millones de nucleótidos. El eje vertical lista las 6
# corridas del estudio, cada una con su código SRR, y el color indica la
# condición: azul para la cepa silvestre (WT), naranja para el mutante
# EMS y verde para el mutante UV. Todas las corridas están entre 0.18 y
# 0.79 Gb, es decir, en el mismo orden de magnitud: ninguna muestra quedó
# con una cantidad de datos mucho menor que las demás, así que ninguna
# se descarta por falta de profundidad antes del análisis. La corrida
# más grande es UV_rep1 (0.79 Gb) y la más pequeña es WT_rep2 (0.18 Gb),
# que es la que se usa como demostración de descarga en este módulo.
# En un análisis de expresión diferencial, esta diferencia de tamaño
# entre muestras se corrige después con la normalización (Módulo 07,
# DESeq2), así que no es un problema, pero sí es algo que hay que
# conocer antes de analizar.

ggsave(
  "02_bases_por_corrida.png",
  plot = grafico_bases,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 7. Gráfico 2: longitud media de lectura de cada corrida
# --------------------------------------------------------------------------
# La longitud media de lectura se calcula de los datos: bases totales
# divididas entre el número de lecturas. En datos paired-end cada spot
# produce 2 lecturas, así que lecturas = spots x 2.

grafico_longitud <- ggplot(tabla_muestras, aes(
  x = reorder(corrida, longitud_lectura_nt),
  y = longitud_lectura_nt,
  fill = condicion
)) +
  geom_col() +
  geom_text(
    aes(label = paste0(round(longitud_lectura_nt, 1), " nt")),
    hjust = -0.1, size = 3.5
  ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c(
    "WT" = "#4E79A7", "EMS" = "#F28E2B", "UV" = "#59A14F"
  )) +
  labs(
    title = "Longitud media de lectura por corrida (PRJNA378177)",
    subtitle = "Calculada como bases totales entre número de lecturas",
    x = "Corrida (código SRR del SRA)",
    y = "Longitud media de lectura (nucleótidos, nt)",
    fill = "Condición"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_longitud)

# Interpretación:
# El eje horizontal mide la longitud media de lectura en nucleótidos
# (nt), es decir, cuántas bases tiene en promedio cada lectura
# secuenciada. El eje vertical lista las 6 corridas y el color indica
# la condición. Las 6 corridas miden alrededor de 75 nt por lectura,
# prácticamente lo mismo en todas: es lo esperado, porque las 6 se
# secuenciaron en el mismo instrumento (Illumina MiSeq) con el mismo
# protocolo de lecturas paired-end de 75 nt por extremo. Que la
# longitud sea igual en todas las muestras confirma que los datos son
# comparables entre sí y que ninguna corrida se secuenció con un
# protocolo distinto. En el Módulo 03 (FastQC) se revisará la calidad
# base por base de estas lecturas de 75 nt.

ggsave(
  "02_longitud_lectura.png",
  plot = grafico_longitud,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 8. Resumen del módulo
# --------------------------------------------------------------------------
# Archivos generados (guárdalos: los usaremos en los módulos siguientes):
#   - 02_muestras_chlorella_sra.csv : tabla de las 6 muestras del proyecto
#     (muestra, condición, réplica, experimento SRX, corrida SRR, diseño,
#     instrumento, spots, bases, gigabases, lecturas y longitud de lectura).
#   - 02_descargar_sra.sh           : guion con los comandos prefetch y
#     fasterq-dump para descargar las 6 corridas.
#   - 02_bases_por_corrida.png      : gráfico del tamaño de cada corrida.
#   - 02_longitud_lectura.png       : gráfico de la longitud de lectura.
#   - fastq/SRR5314137_1.fastq.gz y fastq/SRR5314137_2.fastq.gz:
#     archivos FASTQ de la corrida de demostración, solo si el SRA
#     Toolkit está instalado.
#
# Qué sigue (Módulo 03):
#   Control de calidad de las lecturas descargadas con FastQC y MultiQC,
#   antes de limpiarlas y analizarlas.

# --------------------------------------------------------------------------
# Qué aprendimos
# --------------------------------------------------------------------------
# - El SRA guarda las lecturas crudas de secuenciación en tres niveles:
#   BioProject (PRJNA, el proyecto), experimento (SRX, una muestra) y
#   corrida (SRR, los archivos de lecturas de esa muestra).
# - Los datos de RNA-seq se buscan desde R con rentrez en la base de
#   datos "sra", igual que en el Módulo 01 se buscaron artículos y
#   secuencias. El registro completo de cada experimento se pide con
#   entrez_fetch en formato XML y de ahí se extraen la muestra, la
#   corrida y su tamaño (spots y bases).
# - Para todo el proyecto se eligió el estudio PRJNA378177: 6 muestras
#   de Chlorella vulgaris (WT, EMS y UV, 2 réplicas por condición),
#   pequeñas (0.18 a 0.79 gigabases por corrida) y sobre producción de
#   lípidos, el tema del proyecto.
# - La descarga se hace en la terminal con el SRA Toolkit: prefetch
#   descarga la corrida y fasterq-dump la convierte a archivos FASTQ
#   con --split-files (separa R1 y R2) y --gzip (comprime la salida).
# - Cada corrida se describe con dos números calculados de los datos:
#   spots (fragmentos secuenciados) y bases (nucleótidos totales). De
#   ellos se derivan las gigabases por corrida y la longitud media de
#   lectura (alrededor de 75 nt en las 6 corridas, porque se usó el
#   mismo instrumento y el mismo protocolo).
# - La tabla de muestras (02_muestras_chlorella_sra.csv) es el archivo
#   que ordena el proyecto: cada módulo siguiente la usará para saber
#   qué archivo corresponde a cada condición y cada réplica.
