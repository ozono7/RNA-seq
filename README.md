# RNA-seq — Análisis de RNA-seq sin supercomputadora

Serie avanzada, en R con Bioconductor, para pasar de la búsqueda en bases
de datos públicas a un análisis de RNA-seq completo, a nivel de
publicación, sin supercomputadora.

**Autor:** Dr. José Juan Torres Martínez
**Docente de Ciencias Químicas, UASLP, México**

Material académico para estudiantes.

## Organismo de trabajo

*Chlorella vulgaris* (microalga verde, TaxID 3077), con énfasis en la
acumulación de lípidos y la biorremediación.

## Herramientas

- **Lenguaje:** R, con paquetes de Bioconductor: DESeq2, edgeR, limma,
  SummarizedExperiment, GenomicRanges y Biostrings.
- **Datos de ejemplo:** conjuntos de ExperimentData de Bioconductor, en
  especial `airway` (RNA-seq, conjunto clásico para DESeq2). Cuando un
  módulo necesite datos de ejemplo, se usa uno de Bioconductor
  ExperimentData o del SRA.
- **Bases de datos:** PubMed, GenBank y UniProt.

Para instalar Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "edgeR", "limma", "SummarizedExperiment",
                       "GenomicRanges", "Biostrings", "airway"))
```

## Módulos

La Lección 00 es la introducción a la transcriptómica y se hace antes
del Módulo 01. El análisis diferencial completo con DESeq2 es el
Módulo 07.

| # | Archivo | Tema |
|---|---------|------|
| 00 | `00-introduccion-a-la-transcriptomica.R` | Introducción a la transcriptómica: qué mide, la matriz de conteos y un primer análisis de exploración con el conjunto `airway` (tamaño de librería, normalización básica y PCA) |
| 01 | `01-busqueda-pubmed-genbank-uniprot.R` | Búsqueda de información pública: PubMed, GenBank y UniProt |

Próximo módulo: 02 — Descargar datos reales de RNA-seq del SRA
(Sequence Read Archive) para *Chlorella vulgaris*, y revisarlos antes
de analizarlos.

## Reglas

- Todo calculado **de cero**: los datos van primero y cada valor se
  deriva de ellos; nunca parámetros pre-cocinados.
- Los resultados se guardan en archivos (CSV y FASTA) para usarlos en
  los módulos siguientes.
- Gráficos claros y bien presentados, bien documentados en español.

## Licencia

Material académico de uso libre para fines educativos.
