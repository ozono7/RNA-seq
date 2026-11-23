#!/bin/bash
# Cuantificación con Salmon de las 6 corridas de RNA-seq de
# Chlorella vulgaris (PRJNA378177), con las lecturas limpias del
# Módulo 04 y el ensamblaje de Trinity del Módulo 05.
# Generado por el Módulo 06 de la serie RNA-seq.
# Uso: bash 06_cuantificar_salmon.sh
# Requisitos: Salmon en el PATH (conda install -c bioconda salmon)
#   o en salmon-linux-x86_64/bin/salmon, y el locale en_US.UTF-8
#   generado en Linux (sudo localedef -i en_US -f UTF-8 en_US.UTF-8).
#   El ensamblaje Trinity.fasta debe estar en la carpeta de trabajo
#   (es ensamblaje/trinity_SRR5314137/Trinity.fasta del Módulo 05;
#   en el proyecto completo se usa el ensamblaje de las 6 corridas).
# Cada corrida produce una carpeta salmon_<corrida>/ con quant.sf.
# Después, en R, tximport junta los seis quant.sf en una sola matriz
# de conteos (ver Módulo 06, sección 5, y Módulo 07).

set -e
export LC_ALL=en_US.UTF-8

if command -v salmon >/dev/null 2>&1; then
  SALMON=salmon
elif [ -x salmon-linux-x86_64/bin/salmon ]; then
  SALMON=salmon-linux-x86_64/bin/salmon
else
  echo "No se encontró Salmon. Instálalo con 'conda install -c bioconda salmon'." >&2
  exit 1
fi

if [ ! -d salmon_index ]; then
  "$SALMON" index -t Trinity.fasta -i salmon_index -k 31
fi

for corrida in SRR5314138 SRR5314139 SRR5314140 SRR5314141 SRR5314136 SRR5314137; do
  if [ ! -f "salmon_${corrida}/quant.sf" ]; then
    "$SALMON" quant -i salmon_index -l A \
      -1 "limpieza/fastp/${corrida}_1.clean.fq.gz" \
      -2 "limpieza/fastp/${corrida}_2.clean.fq.gz" \
      -p 4 --validateMappings --gcBias -o "salmon_${corrida}"
  fi
done

echo "Cuantificación con Salmon terminada: 6 carpetas salmon_SRR* con quant.sf"
