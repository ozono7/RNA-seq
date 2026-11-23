#!/bin/bash
# Cuantificación con Kallisto de las 6 corridas de RNA-seq de
# Chlorella vulgaris (PRJNA378177), con las lecturas limpias del
# Módulo 04 y el ensamblaje de Trinity del Módulo 05.
# Generado por el Módulo 06 de la serie RNA-seq.
# Uso: bash 06_cuantificar_kallisto.sh
# Requisitos: Kallisto en el PATH (conda install -c bioconda kallisto)
#   o en kallisto/kallisto. El ensamblaje Trinity.fasta debe estar en
#   la carpeta de trabajo (ver 06_cuantificar_salmon.sh).
# Cada corrida produce una carpeta kallisto_<corrida>/ con abundance.tsv.
# Después, en R, tximport junta los seis abundance.tsv en una sola
# matriz de conteos (ver Módulo 06, sección 5, y Módulo 07).
# La opción -b 100 pide 100 remuestreos (bootstraps) por muestra para
# estimar la incertidumbre de cada abundancia; con -b 0 la corrida es
# más rápida y no se estima esa incertidumbre.

set -e

if command -v kallisto >/dev/null 2>&1; then
  KALLISTO=kallisto
elif [ -x kallisto/kallisto ]; then
  KALLISTO=kallisto/kallisto
else
  echo "No se encontró Kallisto. Instálalo con 'conda install -c bioconda kallisto'." >&2
  exit 1
fi

if [ ! -f kallisto_index.idx ]; then
  "$KALLISTO" index -i kallisto_index.idx Trinity.fasta
fi

for corrida in SRR5314138 SRR5314139 SRR5314140 SRR5314141 SRR5314136 SRR5314137; do
  if [ ! -f "kallisto_${corrida}/abundance.tsv" ]; then
    "$KALLISTO" quant -i kallisto_index.idx -o "kallisto_${corrida}" \
      -b 100 -t 4 \
      "limpieza/fastp/${corrida}_1.clean.fq.gz" \
      "limpieza/fastp/${corrida}_2.clean.fq.gz"
  fi
done

echo "Cuantificación con Kallisto terminada: 6 carpetas kallisto_SRR* con abundance.tsv"
