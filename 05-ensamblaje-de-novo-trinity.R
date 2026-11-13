# ==========================================================================
# Módulo 05. Ensamblaje de novo del transcriptoma (Trinity)
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
# Serie "RNA-seq" (avanzado) — Módulo 05
#
# Descripción: En los módulos 02 a 04 descargamos las lecturas de
#   RNA-seq de Chlorella vulgaris (estudio PRJNA378177, 6 corridas:
#   WT, EMS y UV), revisamos su calidad y las limpiamos con fastp y
#   Trimmomatic. En este módulo se ensambla el transcriptoma de novo
#   con Trinity: a partir de las lecturas limpias se reconstruyen los
#   transcritos sin usar un genoma de referencia. Se trabaja con una
#   submuestra de 100,000 pares de lecturas de la corrida de
#   demostración (SRR5314137, WT_rep2, la más pequeña del estudio),
#   para que el ensamblaje corra en una computadora normal, y se
#   evalúa el ensamblaje con sus propias estadísticas (número de
#   transcritos, longitud total, N50, contenido GC) y con el
#   porcentaje de lecturas que vuelven a alinearse sobre el
#   ensamblaje. El guion 05_ensamblar_trinity.sh trae el comando para
#   ensamblar las 6 corridas completas del proyecto.
#
# Objetivos de aprendizaje:
#   - Entender qué es el ensamblaje de novo y cuándo se usa en lugar
#     de alinear contra un genoma de referencia.
#   - Explicar cómo funciona Trinity: Inchworm, Chrysalis y Butterfly,
#     y qué es un k-mero.
#   - Tomar una submuestra reproducible de lecturas paired-end con
#     semilla fija y explicar por qué se submuestrea en un equipo
#     normal.
#   - Correr Trinity sobre la submuestra y leer su archivo FASTA de
#     transcritos con el paquete Biostrings de Bioconductor.
#   - Calcular las estadísticas del ensamblaje desde los transcritos:
#     número de transcritos y de genes Trinity, longitud total, N50,
#     longitud media y mediana, y contenido GC.
#   - Medir la calidad del ensamblaje alineando las lecturas de vuelta
#     sobre los transcritos con bowtie2 (incluido en Trinity).
#   - Interpretar la distribución de longitudes y la curva acumulada
#     de bases con el N50 marcado.
#
# Requisitos: R con los paquetes tidyverse, Biostrings y ShortRead
#   (Bioconductor). Instálalos con:
#     install.packages("tidyverse")
#     if (!requireNamespace("BiocManager", quietly = TRUE)) {
#       install.packages("BiocManager")
#     }
#     BiocManager::install(c("Biostrings", "ShortRead"))
#   Además, Trinity 2.15.2 (programa de línea de comandos, fuera de
#   R) y Java para su paso Butterfly. Si Trinity no está instalado, el
#   módulo descarga el paquete oficial FULL y lo compila en trinity/;
#   sin autoconf se compila con make no_bamsifter (ver sección 6).
#   Si Java no está, hay que instalarlo antes, por ejemplo:
#     sudo apt install default-jre-headless   (Ubuntu/Debian)
#     brew install openjdk                    (macOS)
#   jellyfish se usa para la normalización in silico de Trinity; si
#   no está instalado, el módulo corre Trinity con
#   --no_normalize_reads (ver sección 6). Se instala con:
#     conda install -c bioconda jellyfish
#   bowtie2 se usa para alinear las lecturas de vuelta sobre el
#   ensamblaje (sección 10); si no está instalado, el módulo descarga
#   el binario precompilado para Linux, y en macOS se instala con:
#     conda install -c bioconda bowtie2
#   samtools lo exige Trinity al arrancar (para manejar archivos SAM
#   y BAM). Se instala con:
#     conda install -c bioconda samtools
#   salmon lo exige Trinity al arrancar (estima la abundancia de los
#   transcritos para el filtrado final). Se instala con:
#     conda install -c bioconda salmon
#   En Linux, salmon necesita además el locale en_US.UTF-8 generado
#   en el sistema (sin él, su indexado falla); se genera con:
#     sudo localedef -i en_US -f UTF-8 en_US.UTF-8
#     (o sudo locale-gen en_US.UTF-8, si está el paquete locales).
#   En macOS el locale ya existe y no hace falta.
#   Nota: en Linux con el locale en POSIX, los acentos de las gráficas
#   salen mal; en ese caso corre el script con LC_ALL=C.utf8 Rscript.
#   En macOS y RStudio el locale ya es UTF-8 y no hace falta.
#
# Paquetes de Bioconductor usados en este módulo:
#   - Biostrings: lee el archivo FASTA del ensamblaje
#     (readDNAStringSet) y da acceso a las secuencias de los
#     transcritos, sus longitudes (width) y su composición de bases
#     (letterFrequency, alphabetFrequency). Con eso se calculan todas
#     las estadísticas del ensamblaje desde los datos, sin copiar
#     números de ningún reporte. Se instala con BiocManager:
#       BiocManager::install("Biostrings")
#   - ShortRead: lee los archivos FASTQ limpios (readFastq) para
#     tomar la submuestra cuando seqtk no está instalado. Se instala
#     con BiocManager:
#       BiocManager::install("ShortRead")
#   En los módulos siguientes, el archivo FASTA del ensamblaje
#   (Trinity.fasta) es la referencia para cuantificar los transcritos
#   con Salmon o Kallisto (Módulo 06) y para la anotación funcional
#   (Módulo 09).
#
# Organismo de trabajo: Chlorella vulgaris (microalga verde, TaxID 3077).
# Estudio del proyecto: PRJNA378177 (ver Módulo 02).
# ==========================================================================

# --------------------------------------------------------------------------
# 1. Cargar paquetes
# --------------------------------------------------------------------------
library(tidyverse)
library(Biostrings)
library(ShortRead)

# --------------------------------------------------------------------------
# 2. Qué es el ensamblaje de novo y cómo trabaja Trinity
# --------------------------------------------------------------------------
# Un experimento de RNA-seq produce millones de lecturas cortas (aquí
# de 76 nt). El transcriptoma es el conjunto de los ARN mensajeros de
# la célula, y para estudiarlo hay que reconstruir esos ARN completos
# a partir de las lecturas. Hay dos caminos:
#   - Alineamiento contra un genoma de referencia: cada lectura se
#     coloca en el lugar del genoma del que proviene. Es el camino
#     habitual cuando el genoma está bien anotado.
#   - Ensamblaje de novo: no se usa ningún genoma. Los transcritos se
#     reconstruyen uniendo las lecturas que se traslapan entre sí.
#     Es el camino cuando no hay genoma de referencia, cuando el
#     genoma está incompleto o cuando se trabaja con una cepa distinta
#     de la de referencia. Chlorella vulgaris tiene genomas publicados,
#     pero la cepa del estudio (OW-01, ver Módulo 02) no es la cepa de
#     referencia, así que el ensamblaje de novo es el camino que
#     enseña esta serie: el transcriptoma se construye desde las
#     lecturas mismas.
#
# Trinity (Grabherr y colaboradores, 2011; Haas y colaboradores,
# 2013) es el ensamblador de novo de transcriptomas más usado. Trabaja
# en tres pasos, y cada uno tiene un nombre propio:
#   - Inchworm: rompe cada lectura en k-meros, que son todas las
#     subcadenas de longitud k (aquí k = 25 nucleótidos). Cuenta
#     cuántas veces aparece cada k-mero y une los k-meros que se
#     traslapan en contigs, que son las primeras piezas del
#     transcriptoma. Un k-mero es la unidad mínima con la que se
#     trabaja: dos lecturas que comparten un k-mero probablemente
#     vienen del mismo transcrito.
#   - Chrysalis: agrupa los contigs que comparten k-meros y construye
#     un grafo de de Bruijn por cada grupo. Un grafo de de Bruijn es
#     una red donde cada nodo es un k-mero y cada flecha une dos
#     k-meros que se traslapan en k - 1 bases; recorrer el grafo es
#     recorrer el transcrito. Cada grupo de contigs que comparte
#     k-meros es un componente, que corresponde a un gen con sus
#     isoformas.
#   - Butterfly: recorre cada grafo y escribe los caminos posibles
#     como transcritos. Cuando un gen tiene empalmes alternativos, el
#     grafo se ramifica y Butterfly escribe varias isoformas del mismo
#     gen. Por eso el archivo final tiene más transcritos que genes:
#     cada gen (TRINITY_DN*_c*_g*) puede tener varias isoformas
#     (TRINITY_DN*_c*_g*_i*).
# Los nombres de los transcritos de Trinity dicen de qué componente
# vienen: TRINITY_DN1_c0_g1_i1 es el transcrito 1 (i1) del gen 1 (g1)
# del componente 0 (c0) del primer grupo de Inchworm (DN1).
#
# --------------------------------------------------------------------------
# 3. Leer las tablas de los módulos anteriores
# --------------------------------------------------------------------------
# La tabla 02_muestras_chlorella_sra.csv la genera el Módulo 02 y el
# archivo 04_resumen_limpieza.csv lo genera el Módulo 04. Si no están
# en la carpeta, corre primero esos módulos: el ensamblaje necesita
# saber qué corrida es la de demostración y con qué lecturas limpias
# se trabaja.

if (!file.exists("02_muestras_chlorella_sra.csv")) {
  stop("No se encontró 02_muestras_chlorella_sra.csv. ",
       "Corre primero el Módulo 02 (02-descargar-datos-sra.R).")
}
if (!file.exists("04_resumen_limpieza.csv")) {
  stop("No se encontró 04_resumen_limpieza.csv. ",
       "Corre primero el Módulo 04 (04-limpieza-lecturas-trimmomatic-fastp.R).")
}

tabla_muestras <- read_csv(
  "02_muestras_chlorella_sra.csv",
  show_col_types = FALSE
)
resumen_limpieza <- read_csv("04_resumen_limpieza.csv", show_col_types = FALSE)

cat("Limpieza del Módulo 04 (filas de fastp):\n")
print(resumen_limpieza %>% filter(herramienta == "fastp"))

# La corrida de demostración es la más pequeña del estudio (WT_rep2,
# SRR5314137, 0.18 gigabases), igual que en los módulos 03 y 04.
corrida_demo <- tabla_muestras %>%
  slice_min(bases_gb, n = 1) %>%
  pull(corrida)

cat("Corrida de demostración:", corrida_demo, "\n")

# --------------------------------------------------------------------------
# 4. Lecturas limpias de entrada (Módulo 04)
# --------------------------------------------------------------------------
# Las lecturas limpias con fastp están en limpieza/fastp/. Si no
# están, corre primero el Módulo 04: este módulo no vuelve a limpiar,
# usa lo que ya se limpió.

limpio_r1 <- file.path(
  "limpieza", "fastp", paste0(corrida_demo, "_1.clean.fq.gz")
)
limpio_r2 <- file.path(
  "limpieza", "fastp", paste0(corrida_demo, "_2.clean.fq.gz")
)

if (!file.exists(limpio_r1) || !file.exists(limpio_r2)) {
  stop("No se encontraron las lecturas limpias de fastp en limpieza/fastp/. ",
       "Corre primero el Módulo 04 (04-limpieza-lecturas-trimmomatic-fastp.R).")
}
cat("Lecturas limpias de entrada:\n  ", limpio_r1, "\n  ", limpio_r2, "\n")

# --------------------------------------------------------------------------
# 5. Submuestra de 100,000 pares para ensamblar en un equipo normal
# --------------------------------------------------------------------------
# La corrida de demostración tiene 1,223,634 pares de lecturas
# limpias. Ensamblar los 1.2 millones de pares con Trinity necesita
# bastante memoria RAM (varios gigabytes) y más de una hora; en un
# equipo normal eso no siempre cabe. Para aprender el método se
# ensambla una submuestra de 100,000 pares, tomados al azar con
# semilla fija (20261113): cualquiera que repita el módulo con la
# misma semilla obtiene exactamente los mismos 100,000 pares. La
# submuestra conserva los pares juntos (el par i de R1 va con el par
# i de R2) porque Trinity necesita los dos extremos apareados. El
# guion 05_ensamblar_trinity.sh trae el comando para ensamblar las 6
# corridas completas, sin submuestrear, en un equipo con al menos
# 10 GB de RAM libres.
#
# La submuestra se toma con seqtk si está instalado (seqtk sample
# lee el FASTQ en flujo, sin cargar todo en memoria) y, si no, con
# ShortRead en R. seqtk se instala con:
#   conda install -c bioconda seqtk
#   sudo apt install seqtk              (Ubuntu/Debian)

dir.create("ensamblaje", showWarnings = FALSE)

n_pares_submuestra <- 100000
semilla_submuestra <- 20261113

sub_r1 <- file.path(
  "ensamblaje", paste0(corrida_demo, "_sub100k_1.fq.gz")
)
sub_r2 <- file.path(
  "ensamblaje", paste0(corrida_demo, "_sub100k_2.fq.gz")
)

ruta_seqtk <- Sys.which("seqtk")

if (!file.exists(sub_r1) || !file.exists(sub_r2)) {
  if (nzchar(ruta_seqtk)) {
    cat("Tomando la submuestra con seqtk (",
        formatC(n_pares_submuestra, format = "d", big.mark = ","),
        " pares, semilla ", semilla_submuestra, ")...\n", sep = "")
    # seqtk sample escribe sin comprimir; se comprime después con
    # gzip para que Trinity reciba archivos .fq.gz. La salida va a un
    # archivo (stdout = ruta) y el estado de salida es el valor que
    # devuelve system2; no se pone stderr = TRUE porque eso haría que
    # system2 devolviera el FASTQ como texto en lugar del estado.
    sub_r1_sin_gz <- sub(".gz$", "", sub_r1)
    sub_r2_sin_gz <- sub(".gz$", "", sub_r2)
    estado_seqtk_1 <- system2(
      ruta_seqtk,
      args = c(
        "sample", "-s", as.character(semilla_submuestra),
        limpio_r1, formatC(n_pares_submuestra, format = "d")
      ),
      stdout = sub_r1_sin_gz
    )
    estado_seqtk_2 <- system2(
      ruta_seqtk,
      args = c(
        "sample", "-s", as.character(semilla_submuestra),
        limpio_r2, formatC(n_pares_submuestra, format = "d")
      ),
      stdout = sub_r2_sin_gz
    )
    if (estado_seqtk_1 != 0 || estado_seqtk_2 != 0) {
      stop("seqtk falló al tomar la submuestra. Revisa los mensajes de arriba.")
    }
    system2("gzip", args = c("-f", sub_r1_sin_gz, sub_r2_sin_gz))
  } else {
    cat("seqtk no está instalado; la submuestra se toma con ShortRead en R.\n")
    cat("Para instalar seqtk: conda install -c bioconda seqtk\n")
    fq_r1 <- readFastq(dirname(limpio_r1), pattern = basename(limpio_r1))
    fq_r2 <- readFastq(dirname(limpio_r2), pattern = basename(limpio_r2))
    n_total <- length(fq_r1)
    if (length(fq_r2) != n_total) {
      stop("R1 y R2 no tienen el mismo número de lecturas.")
    }
    set.seed(semilla_submuestra)
    indices <- sort(sample(n_total, size = min(n_pares_submuestra, n_total)))
    sub_r1_sin_gz <- sub(".gz$", "", sub_r1)
    sub_r2_sin_gz <- sub(".gz$", "", sub_r2)
    writeFastq(fq_r1[indices], sub_r1_sin_gz, mode = "w")
    writeFastq(fq_r2[indices], sub_r2_sin_gz, mode = "w")
    system2("gzip", args = c("-f", sub_r1_sin_gz, sub_r2_sin_gz))
    rm(fq_r1, fq_r2)
  }
  cat("Submuestra lista:",
      formatC(n_pares_submuestra, format = "d", big.mark = ","),
      "pares en\n  ", sub_r1,
      "\n  ", sub_r2, "\n")
} else {
  cat("Submuestra ya existente; se reutiliza.\n")
}

# Comprobación: los dos archivos de la submuestra tienen el mismo
# número de lecturas y son 100,000 pares. Se lee solo la submuestra
# (100,000 lecturas por archivo), no la corrida completa.
n_sub_r1 <- length(readFastq(
  dirname(sub_r1), pattern = basename(sub_r1)
))
n_sub_r2 <- length(readFastq(
  dirname(sub_r2), pattern = basename(sub_r2)
))
cat("Lecturas en la submuestra: R1 =", n_sub_r1, ", R2 =", n_sub_r2, "\n")
if (n_sub_r1 != n_pares_submuestra || n_sub_r2 != n_pares_submuestra) {
  stop("La submuestra no tiene ",
       formatC(n_pares_submuestra, format = "d", big.mark = ","),
       " pares. ",
       "Borra ", sub_r1, " y ", sub_r2, " y vuelve a correr el módulo.")
}

# --------------------------------------------------------------------------
# 6. Localizar o instalar Trinity
# --------------------------------------------------------------------------
# Trinity se busca primero en el PATH (instalación con conda: conda
# install -c bioconda trinity), después como un programa llamado
# Trinity en la carpeta trinity/ de este proyecto (compilado por una
# corrida anterior de este módulo), y si no aparece en ningún lado, el
# módulo descarga el paquete oficial FULL de la versión 2.15.2 y lo
# compila. La compilación necesita Java (para Butterfly, que ya viene
# compilado como Butterfly.jar), gcc, make y cmake; tarda varios
# minutos y escribe su registro en trinity/compilar-trinity.log. La
# compilación completa (make) necesita además autoconf (autoheader)
# para el programa bamsifter; si autoconf no está instalado, se usa
# make no_bamsifter, que compila lo mismo sin bamsifter. bamsifter es
# una utilidad para archivos BAM que el ensamblaje básico no usa; para
# la compilación completa instala autoconf, por ejemplo:
#   sudo apt install autoconf automake libtool   (Ubuntu/Debian)
#   brew install autoconf automake libtool      (macOS)
# Los programas acompañantes se compilan con make plugins (slclust y
# collectl). jellyfish, bowtie2 y salmon NO vienen en el paquete FULL
# de Trinity 2.15.2: se instalan aparte, con conda
# (conda install -c bioconda jellyfish bowtie2 salmon) o, para
# bowtie2, con el binario precompilado que descarga este módulo.

buscar_trinity <- function() {
  ruta_path <- Sys.which("Trinity")
  if (nzchar(ruta_path)) {
    return(unname(ruta_path))
  }
  candidatos <- list.files(
    "trinity", pattern = "^Trinity$",
    recursive = TRUE, full.names = TRUE
  )
  if (length(candidatos) > 0) {
    return(candidatos[1])
  }
  return("")
}

ruta_trinity <- buscar_trinity()

if (!nzchar(ruta_trinity)) {
  if (!nzchar(Sys.which("java"))) {
    stop("Trinity necesita Java y no se encontró java en el PATH. ",
         "Instala Java (por ejemplo: sudo apt install default-jre-headless) ",
         "y vuelve a correr el módulo.")
  }
  if (!nzchar(Sys.which("make")) || !nzchar(Sys.which("cmake"))) {
    stop("Trinity necesita make y cmake para compilarse y no se ",
         "encontraron en el PATH. Instálalos (por ejemplo: ",
         "sudo apt install build-essential cmake) y vuelve a correr el módulo.")
  }
  dir.create("trinity", showWarnings = FALSE)
  tar_ruta <- file.path("trinity", "trinityrnaseq-v2.15.2.FULL.tar.gz")
  if (!file.exists(tar_ruta)) {
    cat("Descargando Trinity 2.15.2 (paquete FULL)...\n")
    download.file(
      paste0(
        "https://github.com/trinityrnaseq/trinityrnaseq/releases/download/",
        "Trinity-v2.15.2/trinityrnaseq-v2.15.2.FULL.tar.gz"
      ),
      tar_ruta, mode = "wb", quiet = TRUE
    )
  }
  cat("Descomprimiendo Trinity 2.15.2...\n")
  utils::untar(tar_ruta, exdir = "trinity")
  dir_trinity <- list.dirs("trinity", recursive = FALSE) %>%
    grep("trinityrnaseq-", ., value = TRUE)
  if (length(dir_trinity) == 0) {
    stop("No se pudo descomprimir Trinity 2.15.2 en trinity/.")
  }
  dir_trinity <- dir_trinity[1]
  cat("Compilando Trinity en", dir_trinity, "(varios minutos)...\n")
  registro_compilacion <- file.path("trinity", "compilar-trinity.log")
  # La compilación completa necesita autoconf (autoheader) para
  # bamsifter. Sin autoconf se usa make no_bamsifter, que compila el
  # ensamblador completo (Inchworm y Chrysalis) sin esa utilidad.
  objetivo_make <- if (nzchar(Sys.which("autoheader"))) "all" else "no_bamsifter"
  if (objetivo_make == "no_bamsifter") {
    cat("Sin autoconf (autoheader): se compila con make no_bamsifter; ",
        "bamsifter queda fuera, el ensamblaje no lo necesita.\n", sep = "")
  }
  estado_make <- system2(
    "make",
    args = c("-C", dir_trinity, objetivo_make, "-j4"),
    stdout = registro_compilacion
  )
  if (estado_make != 0) {
    stop("La compilación de Trinity falló. Revisa ", registro_compilacion, ".")
  }
  cat("Compilando los programas acompañantes (slclust, collectl)...\n")
  registro_plugins <- file.path("trinity", "compilar-trinity-plugins.log")
  estado_plugins <- system2(
    "make",
    args = c("-C", dir_trinity, "plugins", "-j4"),
    stdout = registro_plugins
  )
  if (estado_plugins != 0) {
    stop("La compilación de los programas acompañantes falló. Revisa ",
         registro_plugins, ".")
  }
  ruta_trinity <- buscar_trinity()
  if (!nzchar(ruta_trinity)) {
    stop("Trinity se compiló pero no se encontró el programa Trinity en ",
         dir_trinity, ".")
  }
  cat("Trinity 2.15.2 instalado en:", ruta_trinity, "\n")
} else {
  cat("Trinity encontrado en:", ruta_trinity, "\n")
}

# Versión de Trinity, para dejarla registrada.
version_trinity <- system2(
  ruta_trinity, args = "--version", stdout = TRUE, stderr = TRUE
)
cat("Versión de Trinity:", paste(version_trinity, collapse = " "), "\n")

# samtools lo exige Trinity al arrancar; sin él, Trinity se detiene
# antes de ensamblar. Se instala con conda install -c bioconda samtools.
if (!nzchar(Sys.which("samtools"))) {
  stop("Trinity necesita samtools en el PATH y no se encontró. ",
       "Instálalo (conda install -c bioconda samtools) y vuelve a ",
       "correr el módulo.")
}
cat("samtools encontrado en:", unname(Sys.which("samtools")), "\n")

# salmon también lo exige Trinity al arrancar (lo usa para estimar la
# abundancia de los transcritos). Se instala con conda install -c
# bioconda salmon.
if (!nzchar(Sys.which("salmon"))) {
  stop("Trinity necesita salmon en el PATH y no se encontró. ",
       "Instálalo (conda install -c bioconda salmon) y vuelve a ",
       "correr el módulo.")
}
cat("salmon encontrado en:", unname(Sys.which("salmon")), "\n")

# salmon necesita el locale en_US.UTF-8 generado en el sistema: sin
# él, su paso de indexado falla con el error
# "locale::facet::_S_create_c_locale name not valid" y Trinity se
# detiene después de ensamblar, en el filtrado por expresión. En
# macOS el locale ya existe; en Linux se genera con:
#   sudo locale-gen en_US.UTF-8     (si está instalado el paquete locales)
#   sudo localedef -i en_US -f UTF-8 en_US.UTF-8
if (Sys.info()[["sysname"]] == "Linux") {
  locales_sistema <- system2("locale", args = "-a", stdout = TRUE, stderr = TRUE)
  if (!any(grepl("^en_US\\.utf8$", locales_sistema))) {
    stop("salmon necesita el locale en_US.UTF-8 y no está generado. ",
         "Genéralo con: sudo localedef -i en_US -f UTF-8 en_US.UTF-8 ",
         "(o sudo locale-gen en_US.UTF-8) y vuelve a correr el módulo.")
  }
  cat("Locale en_US.UTF-8 disponible para salmon.\n")
}

# jellyfish (Marçais y Kingsford, 2011) cuenta los k-meros y lo usa
# Trinity para la normalización in silico de las lecturas: antes de
# ensamblar, Trinity reduce los transcritos muy expresados a una
# cobertura máxima (30x por defecto) para que no dominen el
# ensamblaje y para que quepa en memoria. Si jellyfish no está
# instalado, Trinity se corre con --no_normalize_reads y la
# normalización se omite; en la submuestra de demostración eso no
# cambia el resultado de forma apreciable, porque 100,000 pares no
# tienen transcritos tan dominantes. Para instalar jellyfish:
#   conda install -c bioconda jellyfish
# El guion 05_ensamblar_trinity.sh usa la normalización (el
# comportamiento por defecto de Trinity) y supone jellyfish instalado.
hay_jellyfish <- nzchar(Sys.which("jellyfish"))
if (hay_jellyfish) {
  cat("jellyfish encontrado: la normalización in silico queda activada.\n")
} else {
  cat("jellyfish no está instalado: el ensamblaje se corre con ",
      "--no_normalize_reads (sin normalización in silico).\n", sep = "")
  cat("Para instalarlo: conda install -c bioconda jellyfish\n")
}

# bowtie2 viene incluido en el paquete FULL de Trinity y se usa en la
# sección 10 para alinear las lecturas de vuelta sobre el ensamblaje.
# Se busca primero en el PATH y después dentro de la carpeta de
# Trinity (trinity-plugins/bowtie2-*/).
buscar_bowtie2 <- function(dir_trinity) {
  ruta_path <- Sys.which("bowtie2")
  if (nzchar(ruta_path)) {
    return(list(
      bowtie2 = unname(ruta_path),
      bowtie2_build = unname(Sys.which("bowtie2-build"))
    ))
  }
  if (!is.na(dir_trinity) && nzchar(dir_trinity)) {
    binarios <- list.files(
      file.path(dir_trinity, "trinity-plugins"),
      pattern = "^bowtie2$", recursive = TRUE, full.names = TRUE
    )
    constructores <- list.files(
      file.path(dir_trinity, "trinity-plugins"),
      pattern = "^bowtie2-build$", recursive = TRUE, full.names = TRUE
    )
    if (length(binarios) > 0 && length(constructores) > 0) {
      return(list(bowtie2 = binarios[1], bowtie2_build = constructores[1]))
    }
  }
  # bowtie2 no viene en el paquete FULL de Trinity 2.15.2, así que si
  # no está instalado se descarga el binario precompilado para Linux
  # x86_64 (versión 2.5.4) de su página oficial de releases. En macOS
  # este binario no sirve: instala bowtie2 con conda
  # (conda install -c bioconda bowtie2).
  if (Sys.info()[["sysname"]] == "Linux") {
    dir.create("bowtie2", showWarnings = FALSE)
    zip_ruta <- file.path("bowtie2", "bowtie2-2.5.4-linux-x86_64.zip")
    if (!file.exists(zip_ruta)) {
      cat("Descargando bowtie2 2.5.4 (binario precompilado para Linux)...\n")
      download.file(
        paste0(
          "https://github.com/BenLangmead/bowtie2/releases/download/",
          "v2.5.4/bowtie2-2.5.4-linux-x86_64.zip"
        ),
        zip_ruta, mode = "wb", quiet = TRUE
      )
    }
    utils::unzip(zip_ruta, exdir = "bowtie2")
    binario <- file.path(
      "bowtie2", "bowtie2-2.5.4-linux-x86_64", "bowtie2"
    )
    constructor <- file.path(
      "bowtie2", "bowtie2-2.5.4-linux-x86_64", "bowtie2-build"
    )
    if (file.exists(binario) && file.exists(constructor)) {
      return(list(bowtie2 = binario, bowtie2_build = constructor))
    }
  }
  return(list(bowtie2 = "", bowtie2_build = ""))
}

dir_trinity_actual <- if (nzchar(ruta_trinity)) {
  dirname(ruta_trinity)
} else {
  ""
}
bowtie2 <- buscar_bowtie2(dir_trinity_actual)
if (nzchar(bowtie2$bowtie2)) {
  cat("bowtie2 encontrado en:", bowtie2$bowtie2, "\n")
} else {
  cat("bowtie2 no disponible: la sección 10 se omite y el porcentaje de ",
      "lecturas mapeadas queda sin medir.\n", sep = "")
}

# --------------------------------------------------------------------------
# 7. Guion de ensamblaje para las 6 corridas del proyecto
# --------------------------------------------------------------------------
# El guion 05_ensamblar_trinity.sh tiene los comandos exactos de
# Trinity para las 6 corridas del proyecto, con las lecturas limpias
# completas (sin submuestrear) y con los parámetros del proyecto:
#   --seqType fq      : las lecturas son archivos FASTQ.
#   --left / --right  : los dos extremos paired-end (R1 y R2).
#   --CPU 4           : 4 hilos de procesador.
#   --max_memory 10G  : hasta 10 GB de RAM; ajústalo a tu equipo.
#   --full_cleanup    : borra los archivos intermedios al terminar.
# Uso:
#   bash 05_ensamblar_trinity.sh
# El ensamblaje completo de las 6 corridas necesita un equipo con al
# menos 10 GB de RAM libres y tarda varias horas; el de la submuestra
# de este módulo es el que se corre y se analiza aquí.

lineas_sh <- c(
  "#!/bin/bash",
  "# Ensamblaje de novo con Trinity de las 6 corridas de RNA-seq de",
  "# Chlorella vulgaris (PRJNA378177), con las lecturas limpias del",
  "# Módulo 04. Generado por el Módulo 05 de la serie RNA-seq.",
  "# Uso: bash 05_ensamblar_trinity.sh",
  "# Requisitos: Trinity 2.15.2 en el PATH (conda install -c bioconda",
  "# trinity) o en trinity/trinityrnaseq-2.15.2/Trinity, Java,",
  "# jellyfish (conda install -c bioconda jellyfish) para la",
  "# normalización in silico que Trinity hace por defecto, samtools",
  "# (conda install -c bioconda samtools), y al menos 10 GB de RAM",
  "# libres. Sin autoconf, Trinity se compila con",
  "# make no_bamsifter (ver Módulo 05).",
  "",
  "mkdir -p ensamblaje",
  ""
)

ruta_trinity_sh <- if (grepl("^trinity/", ruta_trinity)) {
  ruta_trinity
} else {
  "Trinity"
}

for (corrida in tabla_muestras$corrida) {
  lineas_sh <- c(
    lineas_sh,
    paste0(
      ruta_trinity_sh, " --seqType fq ",
      "--left limpieza/fastp/", corrida, "_1.clean.fq.gz ",
      "--right limpieza/fastp/", corrida, "_2.clean.fq.gz ",
      "--CPU 4 --max_memory 10G --full_cleanup ",
      "--output ensamblaje/trinity_", corrida
    ),
    ""
  )
}

writeLines(lineas_sh, "05_ensamblar_trinity.sh")
cat("Guion escrito: 05_ensamblar_trinity.sh\n")

# --------------------------------------------------------------------------
# 8. Correr Trinity sobre la submuestra de demostración
# --------------------------------------------------------------------------
# Trinity se corre sobre los 100,000 pares de la submuestra, con
# --max_memory 3G y --CPU 4. Si el ensamblaje ya existe de una corrida
# anterior, se reutiliza y no se vuelve a correr. El registro completo
# queda en ensamblaje/trinity_SRR5314137.log. Con --full_cleanup,
# Trinity borra sus archivos intermedios y deja solo el FASTA de
# transcritos (Trinity.fasta) dentro de la carpeta de salida.

dir_salida_trinity <- file.path(
  "ensamblaje", paste0("trinity_", corrida_demo)
)
fasta_trinity <- file.path(dir_salida_trinity, "Trinity.fasta")
registro_trinity <- file.path(
  "ensamblaje", paste0("trinity_", corrida_demo, ".log")
)

correr_trinity <- !file.exists(fasta_trinity)

if (correr_trinity) {
  cat("Corriendo Trinity sobre la submuestra (100,000 pares)...\n")
  cat("Esto tarda varios minutos; el registro queda en ", registro_trinity, "\n")
  # Sin jellyfish, Trinity no puede normalizar las lecturas y se le
  # pasa --no_normalize_reads (ver sección 6). Con jellyfish
  # instalado, la normalización in silico corre por defecto.
  argumentos_trinity <- c(
    "--seqType", "fq",
    "--left", sub_r1,
    "--right", sub_r2,
    "--CPU", "4",
    "--max_memory", "3G",
    "--full_cleanup",
    "--output", dir_salida_trinity
  )
  if (!hay_jellyfish) {
    argumentos_trinity <- c(argumentos_trinity, "--no_normalize_reads")
  }
  estado_trinity <- system2(
    ruta_trinity,
    args = argumentos_trinity,
    stdout = registro_trinity
  )
  if (estado_trinity != 0 || !file.exists(fasta_trinity)) {
    stop("Trinity falló o no escribió ", fasta_trinity,
         ". Revisa ", registro_trinity, ".")
  }
  cat("Trinity terminó. Transcritos en:", fasta_trinity, "\n")
} else {
  cat("Ensamblaje ya existente; se reutiliza:", fasta_trinity, "\n")
}

# --------------------------------------------------------------------------
# 9. Leer el ensamblaje y calcular sus estadísticas con Biostrings
# --------------------------------------------------------------------------
# El archivo Trinity.fasta tiene un transcrito por registro, con su
# nombre (TRINITY_DN*_c*_g*_i*) y su longitud. Con Biostrings se leen
# las secuencias y se calculan, desde los transcritos mismos:
#   - Número de transcritos y número de genes Trinity: el gen es la
#     parte del nombre hasta _i*; cada gen puede tener varias
#     isoformas (_i1, _i2, ...).
#   - Longitud total del ensamblaje (suma de las longitudes, en nt) y
#     longitud media, mediana, mínima y máxima.
#   - N50: se ordenan los transcritos de mayor a menor longitud y se
#     van sumando sus longitudes; el N50 es la longitud del transcrito
#     en el que la suma acumulada llega a la mitad de la longitud
#     total. Es la medida estándar de contigüidad de un ensamblaje:
#     la mitad de las bases del ensamblaje está en transcritos de al
#     menos N50 nt. L50 es cuántos transcritos hacen falta para llegar
#     a esa mitad.
#   - Contenido GC: porcentaje de bases G o C en todo el ensamblaje.
#     En Chlorella vulgaris el GC ronda el 61-62 % (ver Módulo 03).
#   - Transcritos de más de 1,000 nt: los transcritos largos son los
#     que más probablemente corresponden a genes completos.

transcritos <- readDNAStringSet(fasta_trinity)
nombres_transcritos <- names(transcritos)
longitudes <- width(transcritos)

genes_trinity <- sub("_i[0-9]+$", "", nombres_transcritos)
n_genes <- length(unique(genes_trinity))
n_transcritos <- length(transcritos)

longitud_total <- sum(longitudes)
longitud_media <- mean(longitudes)
longitud_mediana <- median(longitudes)
longitud_minima <- min(longitudes)
longitud_maxima <- max(longitudes)

longitudes_ordenadas <- sort(longitudes, decreasing = TRUE)
acumulado <- cumsum(longitudes_ordenadas)
n50 <- longitudes_ordenadas[which(acumulado >= longitud_total / 2)[1]]
l50 <- which(acumulado >= longitud_total / 2)[1]

n_mayor_1000 <- sum(longitudes >= 1000)
pct_mayor_1000 <- 100 * n_mayor_1000 / n_transcritos

# Contenido GC: se cuentan las G y las C de todos los transcritos y
# se dividen entre la longitud total. alphabetFrequency cuenta las
# letras de cada secuencia; colSums las suma para todo el ensamblaje.
conteo_bases <- colSums(alphabetFrequency(
  transcritos, baseOnly = TRUE, as.prob = FALSE
))
pct_gc <- 100 * (conteo_bases["G"] + conteo_bases["C"]) / longitud_total

cat("Transcritos ensamblados:", n_transcritos, "\n")
cat("Genes Trinity (componentes):", n_genes, "\n")
cat("Longitud total del ensamblaje:", format(longitud_total, big.mark = ","),
    "nt\n")
cat("Longitud media:", round(longitud_media, 1), "nt; mediana:",
    longitud_mediana, "nt\n")
cat("Longitud mínima:", longitud_minima, "nt; máxima:",
    format(longitud_maxima, big.mark = ","), "nt\n")
cat("N50:", n50, "nt; L50:", l50, "transcritos\n")
cat("Transcritos de 1,000 nt o más:", n_mayor_1000, " (",
    round(pct_mayor_1000, 1), " %)\n", sep = "")
cat("Contenido GC del ensamblaje:", round(pct_gc, 1), " %\n")

# Los 10 transcritos más largos, con su gen Trinity.
transcritos_largos <- tibble(
  transcrito = nombres_transcritos,
  gen = genes_trinity,
  longitud_nt = longitudes
) %>%
  arrange(desc(longitud_nt)) %>%
  slice_head(n = 10)
cat("Los 10 transcritos más largos:\n")
print(transcritos_largos)

# Transcritos por clase de longitud: cuántos hay en cada clase y qué
# parte de las bases totales aporta cada clase.
clases_longitud <- tibble(
  transcrito = nombres_transcritos,
  longitud_nt = longitudes
) %>%
  mutate(
    clase = case_when(
      longitud_nt < 500 ~ "< 500 nt",
      longitud_nt < 1000 ~ "500-999 nt",
      longitud_nt < 2000 ~ "1,000-1,999 nt",
      TRUE ~ "2,000 nt o más"
    ),
    clase = factor(
      clase,
      levels = c("< 500 nt", "500-999 nt", "1,000-1,999 nt", "2,000 nt o más")
    )
  ) %>%
  group_by(clase) %>%
  summarise(
    transcritos = n(),
    bases = sum(longitud_nt),
    .groups = "drop"
  ) %>%
  mutate(
    pct_transcritos = 100 * transcritos / sum(transcritos),
    pct_bases = 100 * bases / sum(bases)
  )
cat("Transcritos por clase de longitud:\n")
print(clases_longitud)

# --------------------------------------------------------------------------
# 10. Calidad del ensamblaje: alinear las lecturas de vuelta con bowtie2
# --------------------------------------------------------------------------
# Un ensamblaje bueno debe explicar las lecturas de las que salió:
# si se alinean las lecturas de la submuestra sobre los transcritos,
# la mayoría debe encontrar su lugar. bowtie2 (Langmead y Salzberg,
# 2012) alinea lecturas cortas contra una referencia; aquí la
# referencia es el propio archivo Trinity.fasta. Primero se construye
# el índice con bowtie2-build y después se alinean los dos extremos
# de la submuestra. bowtie2 escribe en su salida el porcentaje de
# lecturas alineadas ("overall alignment rate"), que es el número que
# se guarda. En un ensamblaje completo, un valor por encima del 70 %
# se considera bueno; en la submuestra de demostración el valor es
# menor porque el ensamblaje es parcial por diseño: con 100,000 pares
# no se reconstruyen los transcritos poco expresados y las lecturas
# que vienen de ellos no tienen dónde alinearse. Por eso el criterio
# de la demostración es al menos 50 % de lecturas mapeadas, y el del
# ensamblaje completo de las 6 corridas es al menos 70 %. Si bowtie2
# no está disponible, esta sección se omite y el porcentaje queda sin
# medir.

pct_mapeadas <- NA_real_

if (nzchar(bowtie2$bowtie2) && nzchar(bowtie2$bowtie2_build)) {
  indice_bowtie2 <- file.path("ensamblaje", paste0(corrida_demo, "_trinity"))
  cat("Construyendo el índice bowtie2 del ensamblaje...\n")
  estado_indice <- system2(
    bowtie2$bowtie2_build,
    args = c(fasta_trinity, indice_bowtie2),
    stdout = TRUE, stderr = TRUE
  )
  cat(paste(tail(estado_indice, 3), collapse = "\n"), "\n")
  cat("Alineando la submuestra sobre el ensamblaje con bowtie2...\n")
  sam_temporal <- file.path(
    "ensamblaje", paste0(corrida_demo, "_sub100k_contra_trinity.sam")
  )
  salida_bowtie2 <- system2(
    bowtie2$bowtie2,
    args = c(
      "-p", "4", "-q",
      "-x", indice_bowtie2,
      "-1", sub_r1, "-2", sub_r2,
      "-S", sam_temporal
    ),
    stdout = TRUE, stderr = TRUE
  )
  # bowtie2 escribe el resumen en la salida de error, que system2
  # devuelve mezclada con la estándar. La línea que interesa dice:
  # "88.76% overall alignment rate".
  linea_mapeo <- salida_bowtie2[grepl("overall alignment rate", salida_bowtie2)]
  if (length(linea_mapeo) > 0) {
    pct_mapeadas <- as.numeric(sub(
      "^([0-9.]+)% overall alignment rate.*", "\\1",
      linea_mapeo[length(linea_mapeo)]
    ))
    cat("Lecturas de la submuestra alineadas sobre el ensamblaje:",
        pct_mapeadas, " %\n")
  } else {
    cat("No se encontró la línea de overall alignment rate en la salida de bowtie2.\n")
    cat(paste(salida_bowtie2, collapse = "\n"), "\n")
  }
  # El SAM solo sirvió para medir el porcentaje; se borra para no
  # ocupar disco.
  if (file.exists(sam_temporal)) {
    file.remove(sam_temporal)
  }
} else {
  cat("bowtie2 no disponible: el porcentaje de lecturas mapeadas queda sin medir.\n")
}

# --------------------------------------------------------------------------
# 11. Tabla resumen del ensamblaje
# --------------------------------------------------------------------------
# La tabla 05_resumen_ensamblaje.csv junta las estadísticas del
# ensamblaje y el porcentaje de lecturas mapeadas. Es el registro del
# ensamblaje del proyecto y la usarán los módulos siguientes: el
# Módulo 06 cuantifica los transcritos de este FASTA y el Módulo 09
# los anota.

resumen_ensamblaje <- tibble(
  corrida = corrida_demo,
  pares_submuestra = n_pares_submuestra,
  transcritos = n_transcritos,
  genes_trinity = n_genes,
  bases_totales_nt = longitud_total,
  longitud_media_nt = round(longitud_media, 1),
  longitud_mediana_nt = longitud_mediana,
  longitud_minima_nt = longitud_minima,
  longitud_maxima_nt = longitud_maxima,
  n50_nt = n50,
  l50 = l50,
  transcritos_mayor_1000_nt = n_mayor_1000,
  pct_transcritos_mayor_1000_nt = round(pct_mayor_1000, 1),
  pct_gc = round(pct_gc, 1),
  pct_lecturas_mapeadas = round(pct_mapeadas, 2),
  version_trinity = paste(version_trinity, collapse = " ")
)
print(resumen_ensamblaje)
write_csv(resumen_ensamblaje, "05_resumen_ensamblaje.csv")

# --------------------------------------------------------------------------
# 12. Gráfico 1: distribución de longitudes de los transcritos
# --------------------------------------------------------------------------
# El histograma cuenta cuántos transcritos hay de cada longitud. El
# eje horizontal está en escala logarítmica (log10) porque las
# longitudes van de 200 nt a varios miles de nt y en escala lineal
# todo quedaría aplastado a la izquierda. La línea vertical marca el
# N50. La mayoría de los transcritos de un ensamblaje de novo son
# cortos: son fragmentos parciales o transcritos poco expresados que
# la submuestra no alcanza a reconstruir completos.

grafico_longitudes <- ggplot(
  tibble(longitud_nt = longitudes),
  aes(x = longitud_nt)
) +
  geom_histogram(bins = 60, fill = "#4E79A7", color = "white") +
  geom_vline(
    xintercept = n50, linetype = "dashed", color = "#E15759", linewidth = 1
  ) +
  annotate(
    "text", x = n50 * 1.3, y = 0, label = paste0("N50 = ", n50, " nt"),
    color = "#E15759", size = 4, hjust = 0, vjust = -1
  ) +
  scale_x_log10() +
  labs(
    title = "Distribución de longitudes de los transcritos (SRR5314137)",
    subtitle = paste0(
      "Ensamblaje de novo con Trinity, submuestra de ",
      formatC(n_pares_submuestra, format = "d", big.mark = ","), " pares de lecturas"
    ),
    x = "Longitud del transcrito (nucleótidos, nt; escala logarítmica)",
    y = "Número de transcritos"
  ) +
  theme_minimal(base_size = 12)

print(grafico_longitudes)

# Interpretación:
# El eje horizontal es la longitud del transcrito en nucleótidos
# (nt), en escala logarítmica (cada marca multiplica por 10: 200,
# 1,000 y 10,000 nt), y el eje vertical cuenta cuántos transcritos
# tienen esa longitud. La forma típica de un ensamblaje de novo es
# una montaña a la izquierda: la mayoría de los transcritos mide
# entre 200 y 500 nt, que es la longitud mínima que Trinity reporta
# (200 nt por defecto), y la cola de la derecha, mucho más baja,
# llega hasta los transcritos más largos, de varios miles de nt. La
# línea roja marca el N50: la mitad de las bases del ensamblaje está
# en transcritos de al menos esa longitud. Si la montaña de la
# izquierda fuera casi todo el ensamblaje y el N50 quedara por debajo
# de 300 nt, el ensamblaje estaría fragmentado: habría que ensamblar
# con más lecturas o revisar la limpieza del Módulo 04.

ggsave(
  "05_distribucion_longitudes.png",
  plot = grafico_longitudes,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 13. Gráfico 2: curva acumulada de bases y N50
# --------------------------------------------------------------------------
# Los transcritos se ordenan de mayor a menor longitud y se va
# sumando su longitud: el eje horizontal es el puesto del transcrito
# en ese orden (1 es el más largo) y el eje vertical es el porcentaje
# acumulado de las bases totales del ensamblaje. La curva sube rápido
# al principio (los transcritos largos aportan muchas bases) y se
# aplana al final (los miles de transcritos cortos aportan poco). El
# punto donde la curva cruza el 50 % es el N50; su puesto en el eje
# horizontal es el L50.

curva_acumulada <- tibble(
  puesto = seq_along(longitudes_ordenadas),
  longitud_nt = longitudes_ordenadas,
  pct_bases_acumulado = 100 * acumulado / longitud_total
)

grafico_curva_n50 <- ggplot(
  curva_acumulada,
  aes(x = puesto, y = pct_bases_acumulado)
) +
  geom_line(color = "#4E79A7", linewidth = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#59A14F") +
  geom_vline(xintercept = l50, linetype = "dashed", color = "#E15759") +
  geom_point(
    data = curva_acumulada %>% filter(puesto == l50),
    aes(x = puesto, y = pct_bases_acumulado),
    color = "#E15759", size = 3
  ) +
  annotate(
    "text", x = l50, y = 58,
    label = paste0("N50 = ", n50, " nt (L50 = ", l50, ")"),
    color = "#E15759", size = 4, hjust = 1
  ) +
  scale_x_log10() +
  labs(
    title = "Curva acumulada de bases del ensamblaje (SRR5314137)",
    subtitle = "Transcritos ordenados de mayor a menor longitud",
    x = "Puesto del transcrito, del más largo al más corto (escala logarítmica)",
    y = "Bases acumuladas (% de la longitud total del ensamblaje)"
  ) +
  theme_minimal(base_size = 12)

print(grafico_curva_n50)

# Interpretación:
# El eje horizontal es el puesto del transcrito cuando se ordenan
# de mayor a menor longitud (1 es el transcrito más largo), en escala
# logarítmica, y el eje vertical es el porcentaje acumulado de las
# bases totales del ensamblaje: al llegar al transcrito número 100
# ya se lleva acumulado cierto porcentaje de todas las bases. La
# curva sube muy rápido al principio porque los primeros transcritos
# son largos y cada uno aporta muchas bases; después se aplana,
# porque los miles de transcritos cortos que quedan aportan muy poco
# cada uno. La línea verde marca el 50 % de las bases y la línea roja
# el L50: el número de transcritos que hace falta para llegar a esa
# mitad. El punto rojo es el N50: la longitud del transcrito que
# completa el 50 %. Un ensamblaje contiguo tiene una curva que sube
# rápido y un L50 pequeño (pocos transcritos largos explican la mitad
# del ensamblaje); un ensamblaje fragmentado tiene una curva que sube
# despacio y un L50 grande.

ggsave(
  "05_curva_acumulada_n50.png",
  plot = grafico_curva_n50,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 14. Gráfico 3: transcritos por clase de longitud
# --------------------------------------------------------------------------
# Las barras cuentan cuántos transcritos hay en cada clase de
# longitud y qué porcentaje de las bases totales aporta cada clase.
# Las dos cantidades cuentan cosas distintas: hay muchos transcritos
# cortos, pero las bases están en los transcritos largos.

grafico_clases <- clases_longitud %>%
  pivot_longer(
    cols = c(pct_transcritos, pct_bases),
    names_to = "medida", values_to = "porcentaje"
  ) %>%
  mutate(
    medida = recode(
      medida,
      pct_transcritos = "% de los transcritos",
      pct_bases = "% de las bases"
    )
  ) %>%
  ggplot(aes(x = clase, y = porcentaje, fill = medida)) +
  geom_col(position = "dodge") +
  geom_text(
    aes(label = paste0(round(porcentaje, 1), " %")),
    position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5
  ) +
  scale_fill_manual(values = c(
    "% de los transcritos" = "#4E79A7",
    "% de las bases" = "#F28E2B"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Transcritos por clase de longitud (SRR5314137)",
    subtitle = paste0(
      "Ensamblaje de novo con Trinity, submuestra de ",
      formatC(n_pares_submuestra, format = "d", big.mark = ","), " pares de lecturas"
    ),
    x = "Clase de longitud del transcrito",
    y = "Porcentaje (%)",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(grafico_clases)

# Interpretación:
# El eje horizontal agrupa los transcritos en cuatro clases según su
# longitud (menos de 500 nt, 500 a 999 nt, 1,000 a 1,999 nt, y 2,000
# nt o más) y el eje vertical es un porcentaje (%). Cada clase tiene
# dos barras: la azul dice qué porcentaje de todos los transcritos
# cae en esa clase y la naranja dice qué porcentaje de todas las
# bases del ensamblaje aportan esos transcritos. Lo habitual es que
# la barra azul más alta esté en la clase de menos de 500 nt (la
# mayoría de los transcritos son cortos) y que la barra naranja más
# alta esté en las clases largas (las bases están en pocos
# transcritos largos). Esa diferencia es normal en un ensamblaje de
# novo y es justo lo que mide el N50: las bases se concentran en los
# transcritos largos aunque los transcritos sean en su mayoría
# cortos. Si la clase de menos de 500 nt aportara también la mayoría
# de las bases, el ensamblaje estaría fragmentado.

ggsave(
  "05_transcritos_por_clase_longitud.png",
  plot = grafico_clases,
  width = 9, height = 5.5, dpi = 300, bg = "white"
)

# --------------------------------------------------------------------------
# 15. Veredicto: ¿el ensamblaje sirve?
# --------------------------------------------------------------------------
# Se decide con tres números de la tabla resumen: que Trinity haya
# escrito transcritos (más de 0), el N50 (al menos 300 nt para una
# submuestra de 100,000 pares) y el porcentaje de lecturas que se
# alinean de vuelta sobre el ensamblaje (al menos 50 % en la
# demostración con submuestra; en el ensamblaje completo de las 6
# corridas el criterio es al menos 70 %, ver sección 10). La
# submuestra es parcial por diseño, así que el N50 y el porcentaje de
# mapeo son menores que los del ensamblaje completo: el veredicto
# dice si esta demostración sirve para seguir, no si el
# transcriptoma está completo.

cat("Transcritos ensamblados:", n_transcritos, "\n")
cat("N50:", n50, "nt\n")
if (!is.na(pct_mapeadas)) {
  cat("Lecturas mapeadas de vuelta:", pct_mapeadas, " %\n")
} else {
  cat("Lecturas mapeadas de vuelta: sin medir (bowtie2 no disponible).\n")
}

pasa_ensamblaje <- n_transcritos > 0 && n50 >= 300 &&
  (is.na(pct_mapeadas) || pct_mapeadas >= 50)

if (pasa_ensamblaje) {
  cat("Veredicto: el ensamblaje SIRVE para seguir.\n")
  cat("El archivo Trinity.fasta pasa al Módulo 06 (cuantificación con Salmon o Kallisto).\n")
} else {
  cat("Veredicto: el ensamblaje NO cumple los criterios del proyecto.\n")
  cat("Revisa los gráficos y la tabla 05_resumen_ensamblaje.csv antes de seguir.\n")
  cat("Con la submuestra de 100,000 pares, lo primero es probar el ensamblaje\n")
  cat("completo de las 6 corridas con el guion 05_ensamblar_trinity.sh.\n")
}

# --------------------------------------------------------------------------
# 16. Resumen del módulo
# --------------------------------------------------------------------------
# Archivos generados (guárdalos: los usaremos en los módulos siguientes):
#   - 05_resumen_ensamblaje.csv            : tabla con las estadísticas
#     del ensamblaje (transcritos, genes, longitudes, N50, GC) y el
#     porcentaje de lecturas mapeadas de vuelta.
#   - 05_ensamblar_trinity.sh              : guion con los comandos
#     Trinity para las 6 corridas completas del proyecto.
#   - 05_distribucion_longitudes.png       : gráfico de la distribución
#     de longitudes de los transcritos, con el N50 marcado.
#   - 05_curva_acumulada_n50.png           : gráfico de la curva
#     acumulada de bases, con el N50 y el L50 marcados.
#   - 05_transcritos_por_clase_longitud.png: gráfico de transcritos por
#     clase de longitud, en % de transcritos y % de bases.
#   - ensamblaje/SRR5314137_sub100k_1.fq.gz y
#     ensamblaje/SRR5314137_sub100k_2.fq.gz : submuestra de 100,000
#     pares de lecturas limpias usada para ensamblar.
#   - ensamblaje/trinity_SRR5314137/Trinity.fasta : transcritos
#     ensamblados por Trinity (la referencia del Módulo 06).
#   - ensamblaje/trinity_SRR5314137.log    : registro completo de la
#     corrida de Trinity.
#   - trinity/trinityrnaseq-2.15.2/        : Trinity 2.15.2 compilado
#     por el módulo (programa, bowtie2, jellyfish, salmon), solo si
#     Trinity no estaba instalado, con sus registros en
#     trinity/compilar-trinity.log y trinity/compilar-trinity-plugins.log.
#
# Qué sigue (Módulo 06):
#   Cuantificación de transcritos con Salmon o Kallisto: contar
#   cuántas lecturas provienen de cada transcrito del ensamblaje, en
#   cada muestra del proyecto, para armar la matriz de conteos que
#   usa la expresión diferencial del Módulo 07 (DESeq2).
#
# --------------------------------------------------------------------------
# Qué aprendimos
# --------------------------------------------------------------------------
# - El ensamblaje de novo reconstruye los transcritos a partir de las
#   lecturas, sin genoma de referencia: es el camino cuando no hay
#   genoma, cuando está incompleto o cuando la cepa es distinta de la
#   de referencia, como la OW-01 de Chlorella vulgaris de este
#   proyecto.
# - Trinity trabaja en tres pasos: Inchworm rompe las lecturas en
#   k-meros de 25 nt y los une en contigs; Chrysalis agrupa los
#   contigs que comparten k-meros en un grafo de de Bruijn por
#   componente; Butterfly recorre cada grafo y escribe los transcritos
#   e isoformas de cada gen.
# - En una computadora normal se ensambla una submuestra (aquí
#   100,000 pares, semilla fija 20261113) y el ensamblaje completo de
#   las 6 corridas queda en el guion 05_ensamblar_trinity.sh, para un
#   equipo con al menos 10 GB de RAM libres.
# - El nombre de un transcrito Trinity dice de qué gen viene:
#   TRINITY_DN*_c*_g*_i* es la isoforma i del gen g del componente c;
#   por eso hay más transcritos que genes.
# - El N50 es la longitud del transcrito en el que la suma acumulada
#   de las longitudes, ordenadas de mayor a menor, llega a la mitad
#   de la longitud total: la mitad de las bases del ensamblaje está
#   en transcritos de al menos N50 nt. El L50 es cuántos transcritos
#   hacen falta para llegar a esa mitad.
# - La calidad del ensamblaje se mide alineando las lecturas de
#   vuelta sobre los transcritos con bowtie2: si la mayoría de las
#   lecturas encuentra su lugar, el ensamblaje explica los datos de
#   los que salió. En la submuestra el valor es menor que en el
#   ensamblaje completo, porque hay transcritos que la submuestra no
#   alcanza a reconstruir.
# - En R, el paquete Biostrings de Bioconductor lee el FASTA del
#   ensamblaje y permite calcular todas sus estadísticas desde los
#   transcritos mismos: longitudes con width, contenido GC con
#   letterFrequency y alphabetFrequency, y el N50 ordenando las
#   longitudes y acumulando.
# - El archivo Trinity.fasta es la referencia que usa el Módulo 06
#   para cuantificar los transcritos en cada muestra del proyecto.
