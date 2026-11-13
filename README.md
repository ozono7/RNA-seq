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
| 02 | `02-descargar-datos-sra.R` | Descargar datos de RNA-seq del SRA con SRA Toolkit (prefetch y fasterq-dump): estudio PRJNA378177, 6 muestras de *Chlorella vulgaris* (WT, EMS y UV) |
| 03 | `03-control-calidad-fastqc-multiqc.R` | Control de calidad de lecturas con FastQC y MultiQC: calidad Phred por base y por lectura, contenido GC, composición de bases y longitud, con el paquete ShortRead de Bioconductor (corrida de demostración SRR5314137) |
| 04 | `04-limpieza-lecturas-trimmomatic-fastp.R` | Limpieza de lecturas con Trimmomatic y fastp: quitar adaptadores, recortar los extremos de baja calidad y filtrar por longitud mínima, con verificación del antes y el después (corrida de demostración SRR5314137) |
| 05 | `05-ensamblaje-de-novo-trinity.R` | Ensamblaje de novo del transcriptoma con Trinity: reconstruir los transcritos sin genoma de referencia, estadísticas del ensamblaje (N50, GC) y porcentaje de lecturas mapeadas de vuelta con bowtie2 (corrida de demostración SRR5314137, submuestra de 100,000 pares) |

Próximo módulo: 06 — Cuantificación de transcritos con Salmon o Kallisto.

## Reglas

- Todo calculado **de cero**: los datos van primero y cada valor se
  deriva de ellos; nunca parámetros pre-cocinados.
- Los resultados se guardan en archivos (CSV y FASTA) para usarlos en
  los módulos siguientes.
- Gráficos claros y bien presentados, bien documentados en español.

## Licencia

Material académico de uso libre para fines educativos.

