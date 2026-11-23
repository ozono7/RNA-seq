# ==========================================================================
# Módulo 06. Cuantificación de transcritos (Salmon / Kallisto)
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
# Serie "RNA-seq" (avanzado) — Módulo 06
#
# Descripción: En el Módulo 05 ensamblamos el transcriptoma de novo de
#   Chlorella vulgaris con Trinity (2,488 transcritos, a partir de una
#   submuestra de 100,000 pares de lecturas de la corrida SRR5314137,
#   WT_rep2, la más pequeña del estudio PRJNA378177). En este módulo se
#   cuenta cuántas lecturas provienen de cada transcrito: la
#   cuantificación. Se usan dos herramientas que hacen lo mismo por
#   caminos distintos, Salmon y Kallisto, se comparan sus resultados y
#   se construye la matriz de conteos y la matriz de TPM que el
#   Módulo 07 usará para la expresión diferencial con DESeq2.
#
#   Cuantificar no es alinear lectura por lectura contra el
#   transcriptoma: Salmon y Kallisto reparten cada lectura entre los
#   transcritos compatibles con un modelo estadístico (máxima
#   verosimilitud con el algoritmo EM), y por eso son mucho más
#   rápidos que un alineador clásico. El resultado por transcrito son
#   dos números: el conteo estimado de lecturas (NumReads) y el TPM
#   (transcritos por millón), una medida normalizada que permite
#   comparar la abundancia entre transcritos de distinta longitud y
#   entre muestras de distinto tamaño de librería.
#
# Objetivos de aprendizaje:
#   - Explicar qué es la cuantificación de transcritos y por qué el
#     conteo crudo no basta para comparar transcritos ni muestras.
#   - Explicar qué es el TPM, cómo se calcula y por qué sus valores
#     suman un millón en cada muestra.
#   - Construir el índice de Salmon a partir del FASTA de Trinity y
#     cuantificar las lecturas de la submuestra con salmon quant.
#   - Explicar qué es el pseudoalineamiento de Kallisto, construir su
#     índice y cuantificar las mismas lecturas con kallisto quant.
#   - Importar los resultados a R con tximport, el paquete de
#     Bioconductor que une Salmon y Kallisto con DESeq2, y resumir
#     los transcritos a nivel de gen con una tabla tx2gene.
#   - Construir la matriz de conteos y la matriz de TPM desde los
#     archivos quant.sf y abundance.tsv, y guardarlas en CSV.
#   - Comparar Salmon y Kallisto transcrito a transcrito y juzgar si
#     coinciden lo suficiente como para usar cualquiera de los dos.
#   - Interpretar la distribución de TPM, los transcritos más
#     expresados y la relación entre conteos, TPM y longitud.
#
# Requisitos: R con los paquetes tidyverse, Biostrings y tximport
#   (Bioconductor). Instálalos con:
#     install.packages("tidyverse")
#     if (!requireNamespace("BiocManager", quietly = TRUE)) {
#       install.packages("BiocManager")
#     }
#     BiocManager::install(c("Biostrings", "tximport"))
#   Además, Salmon 1.12.1 y Kallisto 0.50.1 (programas de línea de
#   comandos, fuera de R). Se instalan con conda:
#     conda install -c bioconda salmon kallisto
#   o se descargan los binarios oficiales:
#     Salmon:  https://github.com/COMBINE-lab/salmon/releases
#              (salmon-1.12.1_linux_x86_64.tar.gz)
#     Kallisto: https://github.com/pachterlab/kallisto/releases
#              (kallisto_linux-v0.50.1.tar.gz)
#   En Linux, Salmon necesita además el locale en_US.UTF-8 generado
#   en el sistema (sin él, su indexado falla con un error de locale);
#   se genera con:
#     sudo localedef -i en_US -f UTF-8 en_US.UTF-8
#     (o sudo locale-gen en_US.UTF-8, si está el paquete locales).
#   En macOS el locale ya existe y no hace falta. El script pone
#   LC_ALL=en_US.UTF-8 al llamar a Salmon; si tu sistema no tiene
#   ese locale, genera el locale o quita esa variable de entorno.
#   Nota: en Linux con el locale en POSIX, los acentos de las
#   gráficas salen mal; en ese caso corre el script con
#   LC_ALL=C.utf8 Rscript. En macOS y RStudio el locale ya es UTF-8
#   y no hace falta.
#
# Datos de entrada (se generan en los módulos anteriores y deben
#   estar junto a este script):
#   - Trinity.fasta: el ensamblaje del Módulo 05
#     (ensamblaje/trinity_SRR5314137/Trinity.fasta).
#   - SRR5314137_sub100k_1.fq.gz y SRR5314137_sub100k_2.fq.gz: la
#     submuestra de 100,000 pares de lecturas limpias del Módulo 05
#     (ensamblaje/SRR5314137_sub100k_1.fq.gz).
#
# Archivos que produce este módulo:
#   - salmon_index/           índice de Salmon del ensamblaje
#   - salmon_SRR5314137/      cuantificación con Salmon (quant.sf)
#   - kallisto_index.idx      índice de Kallisto del ensamblaje
#   - kallisto_SRR5314137/    cuantificación con Kallisto (abundance.tsv)
#   - 06_matriz_conteos.csv   conteos estimados por transcrito (Salmon)
#   - 06_matriz_tpm.csv       TPM por transcrito (Salmon)
#   - 06_resumen_cuantificacion.csv  resumen de ambas cuantificaciones
#   - 06_distribucion_tpm.png, 06_top_transcritos.png,
#     06_conteos_vs_tpm.png, 06_longitud_vs_tpm.png,
#     06_salmon_vs_kallisto.png
#
# ==========================================================================

# --------------------------------------------------------------------------
# 0. Preparación
# --------------------------------------------------------------------------
library(tidyverse)
library(Biostrings)

# En Linux, Salmon necesita el locale en_US.UTF-8 (ver encabezado).
if (.Platform$OS.type == "unix" && Sys.info()["sysname"] == "Linux") {
  Sys.setenv(LC_ALL = "en_US.UTF-8")
}

# Archivos de entrada: el ensamblaje y la submuestra del Módulo 05.
fasta_transcritos <- "Trinity.fasta"
lecturas_1 <- "SRR5314137_sub100k_1.fq.gz"
lecturas_2 <- "SRR5314137_sub100k_2.fq.gz"
corrida <- "SRR5314137"

for (archivo in c(fasta_transcritos, lecturas_1, lecturas_2)) {
  if (!file.exists(archivo)) {
    stop("Falta el archivo ", archivo,
         ". Corre primero el Módulo 05 y copia aquí Trinity.fasta y la ",
         "submuestra SRR5314137_sub100k_1.fq.gz / SRR5314137_sub100k_2.fq.gz.")
  }
}

# Buscar los programas Salmon y Kallisto: primero en el PATH, luego en
# las carpetas locales habituales (binario descargado junto al script).
buscar_binario <- function(nombre, rutas_locales) {
  en_path <- Sys.which(nombre)
  if (nzchar(en_path)) return(unname(en_path))
  for (ruta in rutas_locales) {
    if (file.exists(ruta)) return(normalizePath(ruta))
  }
  NA_character_
}

salmon <- buscar_binario("salmon", c("salmon-linux-x86_64/bin/salmon",
                                     "salmon/bin/salmon"))
if (is.na(salmon)) {
  stop("No se encontró Salmon. Instálalo con ",
       "'conda install -c bioconda salmon' o descarga el binario de ",
       "https://github.com/COMBINE-lab/salmon/releases")
}
kallisto <- buscar_binario("kallisto", c("kallisto/kallisto",
                                         "kallisto_linux-v0.50.1/kallisto"))
tiene_kallisto <- !is.na(kallisto)
if (!tiene_kallisto) {
  message("No se encontró Kallisto: se cuantifica solo con Salmon. ",
          "Para comparar los dos programas, instala Kallisto con ",
          "'conda install -c bioconda kallisto'.")
}

# Leer el ensamblaje para conocer sus transcritos y sus longitudes.
transcritos <- readDNAStringSet(fasta_transcritos)
# Biostrings guarda en el nombre el encabezado completo del FASTA
# ("TRINITY_DN66_c0_g1_i1 len=243 path=[0:0-242]"); quedarse solo con el ID.
id_transcrito <- sub(" .*", "", names(transcritos))
longitud_nt <- width(transcritos)
cat("Transcritos en el ensamblaje:", length(transcritos), "\n")
cat("Longitud media:", round(mean(longitud_nt), 1), "nt; ",
    "mediana:", median(longitud_nt), "nt\n", sep = "")

# --------------------------------------------------------------------------
# 1. Qué es la cuantificación y qué es el TPM
# --------------------------------------------------------------------------
# Después del ensamblaje tenemos la lista de transcritos, pero todavía
# no sabemos cuánto se expresa cada uno. Cuantificar es contar, para
# cada transcrito, cuántas lecturas provienen de él. El problema es
# que una lectura corta puede ser compatible con varios transcritos
# parecidos (isoformas del mismo gen, o genes de la misma familia), y
# entonces no se sabe a cuál pertenece. Salmon y Kallisto resuelven
# ese reparto con un modelo estadístico: parten de un reparto
# inicial, calculan qué tan probable es cada lectura bajo ese reparto
# y lo corrigen, una y otra vez, con el algoritmo EM
# (expectation-maximization), hasta que el reparto deja de cambiar.
# El resultado es un conteo estimado por transcrito, que ya no tiene
# que ser un número entero.
#
# El conteo crudo no sirve para comparar transcritos entre sí ni
# muestras entre sí, por dos razones:
#   1. Longitud: un transcrito largo genera más lecturas que uno
#      corto aunque se exprese igual, simplemente porque hay más
#      posiciones de donde puede salir una lectura.
#   2. Tamaño de librería: una muestra con 10 millones de lecturas
#      tiene conteos mayores que una con 5 millones, aunque la
#      expresión relativa sea la misma.
# El TPM (transcripts per million, transcritos por millón) corrige
# las dos cosas: primero divide cada conteo por la longitud efectiva
# del transcrito (en kilobases) y luego escala los valores para que
# sumen un millón en cada muestra. Así, un TPM de 100 significa que,
# de cada millón de transcritos de la muestra, 100 corresponden a ese
# transcrito. Dos muestras con distinto tamaño de librería se pueden
# comparar en TPM, y dos transcritos de distinta longitud también.
# La fórmula, para el transcrito i, es:
#   TPM_i = 1,000,000 * (conteo_i / longitud_efectiva_i) /
#           suma_j (conteo_j / longitud_efectiva_j)
# La longitud efectiva es la longitud del transcrito menos la
# longitud media del fragmento, más uno: es el número de posiciones
# desde donde puede empezar un fragmento que caiga completo dentro
# del transcrito. Por eso Salmon y Kallisto reportan EffectiveLength
# junto al TPM.
#
# Salmon y Kallisto hacen lo mismo por caminos distintos:
#   - Salmon construye un índice del transcriptoma y usa
#     quasi-mapping (mapeo ligero, antes llamado quasi-mapping con
#     RapMap; desde la versión 1.0 usa mapeo selectivo sobre un
#     índice pufferfish): para cada lectura encuentra los transcritos
#     compatibles y su posición aproximada, sin alinear base por
#     base. Además corrige dos sesgos: el sesgo de secuencia (los
#     fragmentos no se muestrean igual en todas las posiciones) y el
#     sesgo de GC (los fragmentos con mucho o poco GC se amplifican
#     peor en la PCR).
#   - Kallisto usa pseudoalineamiento: descompone cada lectura en
#     k-meros y busca en qué transcritos aparecen esos k-meros, sin
#     alinear la lectura. Solo le importa el conjunto de transcritos
#     compatibles con la lectura, no su posición exacta. Es todavía
#     más rápido, a cambio de no dar posiciones de alineamiento.
# En la práctica, con datos buenos, los dos dan casi lo mismo; por
# eso este módulo corre los dos y los compara transcrito a
# transcrito. Si difirieran mucho, habría que revisar el ensamblaje
# o la calidad de las lecturas antes de seguir.

# --------------------------------------------------------------------------
# 2. Salmon: índice del transcriptoma
# --------------------------------------------------------------------------
# El índice es una estructura que permite buscar k-meros del
# transcriptoma en milisegundos. Se construye una sola vez por
# ensamblaje y se reutiliza para cuantificar todas las muestras.
# La opción -k 31 es la longitud del k-mero: 31 es el valor
# recomendado para lecturas cortas de Illumina.
if (!dir.exists("salmon_index")) {
  message("Construyendo el índice de Salmon (esto tarda unos segundos)...")
  estado <- system2(salmon, c("index", "-t", fasta_transcritos,
                              "-i", "salmon_index", "-k", "31"))
  if (estado != 0) {
    stop("Salmon index falló. En Linux, genera el locale en_US.UTF-8 con ",
         "'sudo localedef -i en_US -f UTF-8 en_US.UTF-8' y vuelve a correr.")
  }
} else {
  message("El índice salmon_index ya existe; se reutiliza.")
}

# --------------------------------------------------------------------------
# 3. Salmon: cuantificación de la submuestra
# --------------------------------------------------------------------------
# salmon quant lee las lecturas paired-end (-1 y -2), detecta sola la
# orientación de la librería (-l A, de automatic) y estima la
# abundancia de cada transcrito. Opciones usadas:
#   --validateMappings: comprueba cada mapeo candidato con una
#     alineación real antes de aceptarlo; más lento, pero más fiable.
#   --gcBias: corrige el sesgo de GC descrito en la sección 1.
#   -p 2: usa 2 hilos de CPU.
# El archivo de resultados es quant.sf: una tabla con una fila por
# transcrito y las columnas Name, Length, EffectiveLength, TPM y
# NumReads. La carpeta aux_info/meta_info.json guarda el resumen de
# la corrida: cuántos fragmentos se procesaron y cuántos se mapearon.
dir_salmon <- paste0("salmon_", corrida)
if (!file.exists(file.path(dir_salmon, "quant.sf"))) {
  message("Cuantificando con Salmon (submuestra de 100,000 pares)...")
  estado <- system2(salmon, c("quant", "-i", "salmon_index", "-l", "A",
                              "-1", lecturas_1, "-2", lecturas_2,
                              "-p", "2", "--validateMappings", "--gcBias",
                              "-o", dir_salmon))
  if (estado != 0) stop("Salmon quant falló (código ", estado, ").")
} else {
  message("La cuantificación de Salmon ya existe; se reutiliza.")
}

# Resumen de la corrida de Salmon, desde meta_info.json.
meta_salmon <- jsonlite::fromJSON(file.path(dir_salmon, "aux_info",
                                           "meta_info.json"))
cat("Salmon: fragmentos procesados:", meta_salmon$num_processed,
    "; mapeados:", meta_salmon$num_mapped,
    "; porcentaje mapeado:", meta_salmon$percent_mapped, "%\n")

# --------------------------------------------------------------------------
# 4. Kallisto: índice y cuantificación de la misma submuestra
# --------------------------------------------------------------------------
# kallisto index construye un índice de Bruijn a partir del mismo
# FASTA, y kallisto quant pseudoalinea las lecturas y estima las
# abundancias con el mismo algoritmo EM. La opción -b 0 pide cero
# remuestreos (bootstraps): con remuestreos se estimaría la
# incertidumbre de cada abundancia, pero aquí se omite para que la
# corrida sea rápida. El archivo de resultados es abundance.tsv, con
# las columnas target_id, length, eff_length, est_counts y tpm, el
# equivalente exacto del quant.sf de Salmon.
dir_kallisto <- paste0("kallisto_", corrida)
if (tiene_kallisto) {
  if (!file.exists("kallisto_index.idx")) {
    message("Construyendo el índice de Kallisto...")
    estado <- system2(kallisto, c("index", "-i", "kallisto_index.idx",
                                  fasta_transcritos))
    if (estado != 0) stop("Kallisto index falló (código ", estado, ").")
  }
  if (!file.exists(file.path(dir_kallisto, "abundance.tsv"))) {
    message("Cuantificando con Kallisto (submuestra de 100,000 pares)...")
    estado <- system2(kallisto, c("quant", "-i", "kallisto_index.idx",
                                  "-o", dir_kallisto, "-b", "0", "-t", "2",
                                  lecturas_1, lecturas_2))
    if (estado != 0) stop("Kallisto quant falló (código ", estado, ").")
  }
  info_kallisto <- jsonlite::fromJSON(file.path(dir_kallisto,
                                                "run_info.json"))
  version_kallisto <- system2(kallisto, "version", stdout = TRUE)[1]
  version_kallisto <- sub("^kallisto, version ", "", version_kallisto)
  cat("Kallisto: lecturas procesadas:", info_kallisto$n_processed,
      "; pseudoalineadas:", info_kallisto$n_pseudoaligned,
      "; porcentaje:", info_kallisto$p_pseudoaligned, "%\n")
} else {
  message("Sin Kallisto instalado, la comparación Salmon-Kallisto se omite.")
}

# --------------------------------------------------------------------------
# 5. Importar los resultados a R con tximport (Bioconductor)
# --------------------------------------------------------------------------
# tximport es el paquete de Bioconductor que lee los archivos de
# Salmon, Kallisto, Sailfish o RSEM y los convierte en las matrices
# que usan DESeq2, edgeR y limma. Se instala con:
#   BiocManager::install("tximport")
# Con type = "salmon" lee quant.sf y devuelve una lista con tres
# matrices: abundance (TPM), counts (conteos estimados) y length
# (longitud efectiva media). Con txOut = TRUE se queda a nivel de
# transcrito; con txOut = FALSE (el valor por defecto) resume los
# transcritos a nivel de gen usando la tabla tx2gene, que dice a qué
# gen pertenece cada transcrito. En Trinity, el gen es la parte del
# ID hasta _gN (por ejemplo, TRINITY_DN66_c0_g1_i1 pertenece al gen
# TRINITY_DN66_c0_g1); las isoformas _i1, _i2 del mismo gen se suman.
# En esta submuestra casi cada gen tiene una sola isoforma, así que
# el resumen a gen cambia poco, pero el paso se hace igual porque en
# el proyecto completo sí habrá genes con varias isoformas.
#
# Si tximport no está instalado, el script lee quant.sf directamente
# con read_tsv: es el mismo cálculo desde cero (las columnas TPM y
# NumReads ya vienen en el archivo), así que el módulo corre igual
# en cualquier equipo.
tx2gene <- tibble(
  TXNAME = id_transcrito,
  GENEID = sub("^(TRINITY_DN[0-9]+_c[0-9]+_g[0-9]+)_i[0-9]+$", "\\1",
               id_transcrito)
)
cat("Genes Trinity distintos:", n_distinct(tx2gene$GENEID),
    " (de ", nrow(tx2gene), " transcritos)\n", sep = "")

archivos_salmon <- c(SRR5314137 = file.path(dir_salmon, "quant.sf"))
names(archivos_salmon) <- corrida

if (requireNamespace("tximport", quietly = TRUE)) {
  message("Importando con tximport (nivel de transcrito)...")
  txi_tx <- tximport::tximport(archivos_salmon, type = "salmon",
                               txOut = TRUE)
  txi_gen <- tximport::tximport(archivos_salmon, type = "salmon",
                                tx2gene = tx2gene, txOut = FALSE)
  matriz_tpm <- as_tibble(txi_tx$abundance, rownames = "transcrito")
  matriz_conteos <- as_tibble(txi_tx$counts, rownames = "transcrito")
  cat("tximport a nivel de gen: ", nrow(txi_gen$counts), " genes; ",
      "suma de TPM por gen: ", round(sum(txi_gen$abundance[, 1]), 1), "\n",
      sep = "")
} else {
  message("tximport no está instalado: leyendo quant.sf directamente.")
  quant <- read_tsv(file.path(dir_salmon, "quant.sf"),
                    show_col_types = FALSE)
  matriz_tpm <- quant %>%
    select(transcrito = Name, !!corrida := TPM)
  matriz_conteos <- quant %>%
    select(transcrito = Name, !!corrida := NumReads)
}

# La tabla completa de Salmon, con longitudes, para los gráficos.
cuant_salmon <- read_tsv(file.path(dir_salmon, "quant.sf"),
                         show_col_types = FALSE) %>%
  dplyr::rename(transcrito = Name, longitud = Length,
         longitud_efectiva = EffectiveLength, tpm = TPM,
         conteos = NumReads) %>%
  mutate(log_tpm = log10(tpm + 1),
         log_conteos = log10(conteos + 1))

# Comprobación: los TPM de una muestra siempre suman un millón.
cat("Suma de TPM (Salmon):", round(sum(cuant_salmon$tpm), 1), "\n")
cat("Suma de conteos estimados (Salmon):",
    round(sum(cuant_salmon$conteos), 1), "\n")
cat("Transcritos con TPM mayor que 0:",
    sum(cuant_salmon$tpm > 0), " de ", nrow(cuant_salmon), "\n")
cat("Transcritos con TPM mayor o igual a 1:",
    sum(cuant_salmon$tpm >= 1), "\n")

# La tabla de Kallisto, en el mismo formato, para la comparación.
if (tiene_kallisto && file.exists(file.path(dir_kallisto, "abundance.tsv"))) {
  cuant_kallisto <- read_tsv(file.path(dir_kallisto, "abundance.tsv"),
                             show_col_types = FALSE) %>%
    dplyr::rename(transcrito = target_id, longitud = length,
           longitud_efectiva = eff_length, tpm = tpm,
           conteos = est_counts)
  cat("Suma de TPM (Kallisto):", round(sum(cuant_kallisto$tpm), 1), "\n")
  cat("Suma de conteos estimados (Kallisto):",
      round(sum(cuant_kallisto$conteos), 1), "\n")
}

# --------------------------------------------------------------------------
# 6. Guardar las matrices y el resumen
# --------------------------------------------------------------------------
# La matriz de conteos y la matriz de TPM son la salida principal del
# módulo: el Módulo 07 las leerá para la expresión diferencial con
# DESeq2. Por eso se guardan en CSV, con una columna por muestra
# (aquí solo la corrida de demostración; el guion
# 06_cuantificar_salmon.sh cuantifica las 6 corridas del proyecto y
# el Módulo 07 juntará sus quant.sf en una sola matriz).
write_csv(matriz_conteos, "06_matriz_conteos.csv")
write_csv(matriz_tpm, "06_matriz_tpm.csv")

resumen <- tibble(
  herramienta = c("Salmon", if (tiene_kallisto) "Kallisto"),
  version = c(meta_salmon$salmon_version,
              if (tiene_kallisto) version_kallisto),
  corrida = corrida,
  fragmentos_procesados = c(meta_salmon$num_processed,
                            if (tiene_kallisto) info_kallisto$n_processed),
  fragmentos_mapeados = c(meta_salmon$num_mapped,
                          if (tiene_kallisto) info_kallisto$n_pseudoaligned),
  pct_mapeados = c(meta_salmon$percent_mapped,
                   if (tiene_kallisto) info_kallisto$p_pseudoaligned),
  transcritos = nrow(cuant_salmon),
  transcritos_tpm_mayor_0 = c(sum(cuant_salmon$tpm > 0),
                             if (tiene_kallisto) sum(cuant_kallisto$tpm > 0)),
  transcritos_tpm_mayor_igual_1 = c(sum(cuant_salmon$tpm >= 1),
                                    if (tiene_kallisto) sum(cuant_kallisto$tpm >= 1)),
  suma_tpm = c(sum(cuant_salmon$tpm),
               if (tiene_kallisto) sum(cuant_kallisto$tpm)),
  suma_conteos = c(sum(cuant_salmon$conteos),
                   if (tiene_kallisto) sum(cuant_kallisto$conteos)),
  longitud_media_fragmento_nt = c(meta_salmon$frag_length_mean,
                                  if (tiene_kallisto) NA_real_)
)
write_csv(resumen, "06_resumen_cuantificacion.csv")
print(resumen)

# --------------------------------------------------------------------------
# 7. Gráficos
# --------------------------------------------------------------------------
# Todos los gráficos se exportan al final de esta sección con ggsave()
# activo, en PNG a 300 dpi.

# --- Gráfico 1: distribución de TPM ---------------------------------------
# La mayoría de los transcritos se expresa poco y unos pocos se
# expresan muchísimo, así que el eje X usa log10(TPM + 1): el +1
# permite incluir los transcritos con TPM = 0, cuyo logaritmo no
# existe. Cada barra cuenta cuántos transcritos caen en ese rango de
# expresión.
g_distribucion <- ggplot(cuant_salmon, aes(x = log_tpm)) +
  geom_histogram(bins = 40, fill = "#2E7D6B", colour = "white") +
  labs(title = "Distribución de la expresión (TPM) en SRR5314137",
       subtitle = "Submuestra de 100,000 pares de lecturas, cuantificación con Salmon",
       x = "log10(TPM + 1)  [TPM: transcritos por millón, adimensional]",
       y = "Número de transcritos") +
  theme_minimal()

# Interpretación: El eje X es log10(TPM + 1), donde TPM son los
# transcritos por millón, una medida adimensional: no tiene unidades
# porque es una proporción escalada a un millón. El eje Y cuenta
# cuántos transcritos caen en cada intervalo de expresión. La
# distribución es unimodal y bastante concentrada: la mayor parte de
# los transcritos cae entre 100 y 1,000 TPM (log10 entre 2 y 3), con
# el máximo de la distribución alrededor de log10(TPM + 1) = 2.5, es
# decir, unos 320 TPM, que es también la mediana. El TPM mínimo es
# 82 y el máximo 5,778, así que toda la distribución cabe en menos
# de dos órdenes de magnitud. No hay transcritos con TPM = 0, y eso
# tiene una explicación: el ensamblaje se construyó a partir de
# estas mismas lecturas (Módulo 05), así que todo transcrito del
# ensamblaje tiene, por construcción, lecturas que lo respaldan. La
# cola de la derecha es corta: solo unos pocos transcritos superan
# los 2,000 TPM. En un transcriptoma de referencia completo, con
# genes que no se expresan en esta condición, la distribución sería
# mucho más ancha y con una barra grande en cero; aquí se ve solo la
# parte expresada del transcriptoma, porque Trinity no ensambla lo
# que no tiene cobertura.

# --- Gráfico 2: los 15 transcritos más expresados --------------------------
top15 <- cuant_salmon %>%
  arrange(desc(tpm)) %>%
  slice_head(n = 15) %>%
  mutate(transcrito = fct_reorder(transcrito, tpm))

g_top <- ggplot(top15, aes(x = tpm, y = transcrito)) +
  geom_col(fill = "#2E7D6B") +
  labs(title = "Los 15 transcritos más expresados en SRR5314137",
       subtitle = "Cuantificación con Salmon, submuestra de 100,000 pares",
       x = "TPM  [transcritos por millón, adimensional]",
       y = "Transcrito (ID de Trinity)") +
  theme_minimal()

# Interpretación: El eje X es el TPM (transcritos por millón,
# adimensional) y el eje Y lista los 15 transcritos con mayor TPM,
# ordenados de mayor a menor. Cada barra es un transcrito del
# ensamblaje de Trinity, identificado por su ID. Los valores van de
# unos 2,000 a 5,778 TPM: el transcrito más expresado representa
# 5,778 de cada millón de transcritos, es decir, el 0.58 % de la
# librería. Ningún transcrito domina la muestra por sí solo, y los
# 15 primeros juntos suman una fracción pequeña del millón: la
# expresión está repartida entre muchos transcritos, no concentrada
# en uno o dos. En el proyecto completo (Módulo 07) se verá si estos
# transcritos cambian entre WT, EMS y UV; aquí solo se describe su
# abundancia en la muestra WT_rep2. Los IDs de Trinity no dicen qué
# gen es cada transcrito: eso lo resolverá la anotación funcional
# del Módulo 09.

# --- Gráfico 3: conteos frente a TPM --------------------------------------
# Si el TPM corrigiera perfectamente la longitud, todos los puntos
# caerían sobre una línea recta. La dispersión alrededor de la recta
# muestra cuánto cambia la relación por la longitud efectiva de cada
# transcrito.
g_conteos_tpm <- ggplot(cuant_salmon,
                        aes(x = log_conteos, y = log_tpm)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#2E7D6B") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40") +
  labs(title = "Conteos estimados frente a TPM (Salmon)",
       subtitle = "Cada punto es un transcrito de SRR5314137",
       x = "log10(conteos estimados + 1)  [conteos: número de lecturas, adimensional]",
       y = "log10(TPM + 1)  [TPM: transcritos por millón, adimensional]") +
  theme_minimal()

# Interpretación: El eje X es log10(conteos estimados + 1) y el eje Y
# es log10(TPM + 1); ambas son medidas adimensionales (conteos de
# lecturas y proporciones por millón). Los puntos forman una nube
# alargada que sigue la diagonal punteada: a más conteos, más TPM,
# como era de esperar. Pero la nube tiene grosor: para un mismo
# número de conteos, el TPM varía según la longitud efectiva del
# transcrito, porque el TPM divide el conteo por la longitud. Los
# puntos por encima de la diagonal son transcritos cortos (mismo
# conteo, más TPM) y los de debajo son transcritos largos. Esa es
# exactamente la corrección que hace el TPM y la razón por la que el
# conteo crudo no se usa para comparar transcritos de distinta
# longitud. No hay puntos en el cero: el transcrito con menos
# lecturas asignadas tiene 2.1 conteos estimados, porque el
# ensamblaje se construyó a partir de estas lecturas y todo
# transcrito tiene cobertura.

# --- Gráfico 4: longitud del transcrito frente a TPM -----------------------
g_longitud_tpm <- ggplot(cuant_salmon,
                         aes(x = longitud, y = log_tpm)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#2E7D6B") +
  labs(title = "Longitud del transcrito frente a su expresión (Salmon)",
       subtitle = "Cada punto es un transcrito de SRR5314137",
       x = "Longitud del transcrito (nt)",
       y = "log10(TPM + 1)  [TPM: transcritos por millón, adimensional]") +
  theme_minimal()

# Interpretación: El eje X es la longitud del transcrito en
# nucleótidos (nt) y el eje Y es log10(TPM + 1), con el TPM en
# transcritos por millón (adimensional). No hay una relación
# monótona clara: hay transcritos cortos y largos en todo el rango
# de expresión, y los transcritos con TPM más alto no son los más
# largos. Eso es lo que se espera después de la normalización: el
# TPM ya corrigió el efecto de la longitud, así que la expresión
# medida no depende del tamaño del transcrito. Si el gráfico se
# hiciera con conteos crudos en el eje Y, sí se vería una tendencia
# ascendente (los transcritos largos acumulan más lecturas). La
# nube se estrecha a la derecha porque el ensamblaje tiene pocos
# transcritos largos (Módulo 05: solo 99 de 2,488 superan los
# 1,000 nt).

# --- Gráfico 5: Salmon frente a Kallisto -----------------------------------
if (tiene_kallisto && exists("cuant_kallisto")) {
  comparacion <- cuant_salmon %>%
    select(transcrito, tpm_salmon = tpm) %>%
    inner_join(cuant_kallisto %>% select(transcrito, tpm_kallisto = tpm),
               by = "transcrito") %>%
    mutate(log_salmon = log10(tpm_salmon + 1),
           log_kallisto = log10(tpm_kallisto + 1))
  cor_pearson <- cor(comparacion$log_salmon, comparacion$log_kallisto)
  cor_spearman <- cor(comparacion$log_salmon, comparacion$log_kallisto,
                      method = "spearman")
  cat("Correlación Salmon-Kallisto (Pearson, log10 TPM + 1):",
      round(cor_pearson, 4), "\n")
  cat("Correlación Salmon-Kallisto (Spearman, log10 TPM + 1):",
      round(cor_spearman, 4), "\n")

  g_comparacion <- ggplot(comparacion,
                          aes(x = log_salmon, y = log_kallisto)) +
    geom_point(alpha = 0.4, size = 1.2, colour = "#2E7D6B") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey40") +
    labs(title = "Salmon frente a Kallisto: TPM por transcrito",
         subtitle = paste0("SRR5314137, Pearson r = ", round(cor_pearson, 3),
                           "; Spearman rho = ", round(cor_spearman, 3)),
         x = "log10(TPM Salmon + 1)  [TPM: transcritos por millón, adimensional]",
         y = "log10(TPM Kallisto + 1)  [TPM: transcritos por millón, adimensional]") +
    theme_minimal()

  # Interpretación: Los dos ejes son log10(TPM + 1), uno con el TPM
  # de Salmon y otro con el de Kallisto; el TPM es adimensional en
  # ambos casos. Cada punto es un transcrito. Los puntos se pegan a
  # la diagonal punteada: los dos programas asignan casi la misma
  # abundancia a casi todos los transcritos. La correlación de
  # Pearson sobre los valores logarítmicos es 0.92 y la de Spearman
  # 0.95; el cociente TPM Salmon / TPM Kallisto cae entre 0.89 y 1.23
  # para el 90 % central de los transcritos, y solo 7 de 2,488
  # difieren en más del doble. La correlación de Pearson no llega a
  # 0.99 porque el rango de TPM aquí es estrecho (de 82 a 5,778, unas
  # 70 veces): con poca dispersión total, las pequeñas diferencias
  # entre programas pesan más en el coeficiente. Los dos puntos
  # pegados al cero en el eje de Kallisto son los dos únicos
  # transcritos a los que Kallisto no asignó ninguna lectura (TPM =
  # 0), mientras Salmon les asignó unas 80 TPM. La conclusión
  # práctica es que, con estos datos, da igual cuantificar con
  # Salmon o con Kallisto: el resultado que llegará a DESeq2 en el
  # Módulo 07 es el mismo en lo esencial. Se usa Salmon como
  # herramienta principal porque Trinity ya lo usa en su filtrado
  # interno y porque corrige el sesgo de GC.
} else {
  g_comparacion <- NULL
}

# Exportar todos los gráficos en PNG a 300 dpi.
ggsave("06_distribucion_tpm.png", g_distribucion, width = 8, height = 5,
       dpi = 300)
ggsave("06_top_transcritos.png", g_top, width = 9, height = 6, dpi = 300)
ggsave("06_conteos_vs_tpm.png", g_conteos_tpm, width = 7, height = 6,
       dpi = 300)
ggsave("06_longitud_vs_tpm.png", g_longitud_tpm, width = 7, height = 6,
       dpi = 300)
if (!is.null(g_comparacion)) {
  ggsave("06_salmon_vs_kallisto.png", g_comparacion, width = 7, height = 6,
         dpi = 300)
}

# --------------------------------------------------------------------------
# 8. Qué aprendimos
# --------------------------------------------------------------------------
#   - Cuantificar es repartir las lecturas entre los transcritos del
#     ensamblaje; Salmon (mapeo ligero) y Kallisto
#     (pseudoalineamiento) lo hacen con el algoritmo EM, sin alinear
#     lectura por lectura, y por eso cuantifican 100,000 pares en
#     segundos en una computadora normal.
#   - El conteo crudo depende de la longitud del transcrito y del
#     tamaño de la librería; el TPM corrige las dos cosas y sus
#     valores suman un millón en cada muestra, así que sirve para
#     comparar transcritos y muestras. En esta submuestra, Salmon
#     mapeó el 57.99 % de los fragmentos y Kallisto pseudoalineó el
#     57.14 %, en línea con el 55.19 % de lecturas que volvieron a
#     alinearse sobre el ensamblaje con bowtie2 en el Módulo 05.
#   - tximport, de Bioconductor, lee quant.sf y abundance.tsv y
#     entrega las matrices de conteos y de TPM que usan DESeq2,
#     edgeR y limma; con una tabla tx2gene resume los transcritos a
#     nivel de gen. Si tximport no está disponible, leer quant.sf
#     directamente da las mismas matrices.
#   - La expresión está repartida: los 2,488 transcritos tienen TPM
#     entre 82 y 5,778 (mediana 332 TPM) y el más expresado es solo
#     el 0.58 % de la librería. No hay transcritos con TPM = 0
#     porque el ensamblaje se construyó a partir de estas mismas
#     lecturas: todo transcrito del ensamblaje tiene cobertura por
#     construcción.
#   - Salmon y Kallisto coinciden en lo esencial en estos datos
#     (correlación de Pearson 0.92 y de Spearman 0.95 sobre log10
#     TPM + 1; solo 7 de 2,488 transcritos difieren en más del
#     doble): cualquiera de los dos sirve, y el módulo usa Salmon
#     como herramienta principal.
#   - Las matrices 06_matriz_conteos.csv y 06_matriz_tpm.csv son la
#     salida del módulo: el Módulo 07 las ampliará a las 6 corridas
#     del proyecto (WT, EMS y UV, dos réplicas cada una) y hará la
#     expresión diferencial con DESeq2.
#
# Para el proyecto completo, los guiones 06_cuantificar_salmon.sh y
# 06_cuantificar_kallisto.sh cuantifican las 6 corridas de
# PRJNA378177 con las lecturas limpias del Módulo 04; después, en R,
# tximport junta los seis quant.sf en una sola matriz de conteos con
# una columna por muestra.
