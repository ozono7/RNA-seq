# ==========================================================================
# Módulo 01. Búsqueda de información pública: PubMed, GenBank y UniProt
# ==========================================================================
# Autor: Dr. José Juan Torres Martínez
# Docente de Ciencias Químicas
# Universidad Autónoma de San Luis Potosí (UASLP)
# México
#
# Material académico para estudiantes
# Serie "RNA-seq" (avanzado) — Módulo 01
#
# Descripción: Primer módulo del proyecto de RNA-seq. Antes de analizar un
#   solo dato, hay que saber qué se ha publicado y qué secuencias existen.
#   En este módulo se busca información pública en tres bases de datos:
#   PubMed (literatura), GenBank (secuencias de nucleótidos) y UniProt
#   (proteínas), todo desde R y sin supercomputadora. Los resultados se
#   guardan en archivos para usarlos en los módulos siguientes.
#
# Objetivos de aprendizaje:
#   - Buscar literatura en PubMed con rentrez (paquete de R para NCBI).
#   - Buscar y descargar secuencias de GenBank.
#   - Consultar UniProt con su API REST.
#   - Guardar los resultados en archivos (CSV y FASTA) para el proyecto.
#
# Requisitos: R con los paquetes tidyverse, rentrez, httr y jsonlite.
#   Instálalos con: install.packages(c("tidyverse", "rentrez", "httr", "jsonlite"))
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

# --------------------------------------------------------------------------
# 2. PubMed: ¿qué se ha publicado sobre Chlorella y lípidos?
# --------------------------------------------------------------------------
# Buscamos artículos sobre acumulación de lípidos en Chlorella vulgaris
# estudiada con RNA-seq. Guardamos título, autores, revista y año.

busqueda_pubmed <- entrez_search(
  db = "pubmed",
  term = "Chlorella vulgaris lipid accumulation RNA-seq",
  retmax = 10
)

# Resumen de cada artículo encontrado
resumen_pubmed <- entrez_summary(
  db = "pubmed",
  id = busqueda_pubmed$ids
)

# Si solo hay un artículo, entrez_summary devuelve un solo resumen y no
# una lista; lo convertimos a lista de uno para que el resto funcione igual.
if (length(busqueda_pubmed$ids) == 1) {
  resumen_pubmed <- list(resumen_pubmed)
}

# Extraer los campos con la función de rentrez (extract_from_esummary),
# que maneja la estructura de los resúmenes correctamente.
tabla_pubmed <- tibble(
  pmid = busqueda_pubmed$ids,
  titulo = extract_from_esummary(resumen_pubmed, "title"),
  revista = extract_from_esummary(resumen_pubmed, "fulljournalname"),
  anio = extract_from_esummary(resumen_pubmed, "pubdate"),
  autores = vapply(resumen_pubmed, function(x) {
    a <- x$authors
    if (is.null(a) || length(a) == 0) return(NA_character_)
    nombres <- vapply(a, function(au) {
      nm <- au$name
      if (is.null(nm)) NA_character_ else as.character(nm)
    }, character(1))
    paste(nombres, collapse = "; ")
  }, character(1))
)

print(tabla_pubmed)
write_csv(tabla_pubmed, "01_pubmed_chlorella_lipidos.csv")

# --------------------------------------------------------------------------
# 3. GenBank: secuencias del gen rbcL de Chlorella vulgaris
# --------------------------------------------------------------------------
# rbcL codifica la subunidad grande de la Rubisco (fotosíntesis). Buscamos
# secuencias de este gen en Chlorella vulgaris y descargamos las primeras
# cinco en formato FASTA.

busqueda_genbank <- entrez_search(
  db = "nucleotide",
  term = "Chlorella vulgaris[Organism] AND rbcL[Gene]",
  retmax = 5
)

# Descargar las secuencias en formato FASTA
secuencias_fasta <- entrez_fetch(
  db = "nucleotide",
  id = busqueda_genbank$ids,
  rettype = "fasta"
)

# Guardar el archivo FASTA para los módulos siguientes
writeLines(secuencias_fasta, "01_rbcL_chlorella.fasta")

# Tabla con los accesos encontrados
tabla_genbank <- tibble(
  acceso = busqueda_genbank$ids
)

print(tabla_genbank)
write_csv(tabla_genbank, "01_genbank_rbcL.csv")

# --------------------------------------------------------------------------
# 4. UniProt: proteínas de Chlorella vulgaris relacionadas con lípidos
# --------------------------------------------------------------------------
# Consultamos la API REST de UniProt. TaxID 3077 = Chlorella vulgaris.
# Pedimos proteínas anotadas con el término "lipid" (10 resultados).

respuesta <- GET(
  "https://rest.uniprot.org/uniprotkb/search",
  query = list(
    query = "taxonomy_id:3077 AND lipid",
    format = "tsv",
    fields = "accession,id,protein_name,organism_name,length",
    size = 10
  )
)

# Leer la respuesta como tabla. I() indica que el texto es el dato mismo
# y no el nombre de un archivo (así evitamos la advertencia de readr).
tabla_uniprot <- content(respuesta, as = "text", encoding = "UTF-8") %>%
  read_tsv(I(.), show_col_types = FALSE)

print(tabla_uniprot)
write_csv(tabla_uniprot, "01_uniprot_chlorella_lipidos.csv")

# --------------------------------------------------------------------------
# 5. Resumen del módulo
# --------------------------------------------------------------------------
# Archivos generados (guárdalos: los usaremos en los módulos siguientes):
#   - 01_pubmed_chlorella_lipidos.csv : artículos sobre Chlorella y lípidos.
#   - 01_rbcL_chlorella.fasta        : secuencias rbcL de Chlorella vulgaris.
#   - 01_genbank_rbcL.csv            : accesos de GenBank encontrados.
#   - 01_uniprot_chlorella_lipidos.csv : proteínas de Chlorella y lípidos.
#
# Qué sigue (Módulo 02):
#   Descargar datos reales de RNA-seq del SRA (Sequence Read Archive)
#   para Chlorella vulgaris, y revisarlos antes de analizarlos.
