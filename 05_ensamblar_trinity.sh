#!/bin/bash
# Ensamblaje de novo con Trinity de las 6 corridas de RNA-seq de
# Chlorella vulgaris (PRJNA378177), con las lecturas limpias del
# Módulo 04. Generado por el Módulo 05 de la serie RNA-seq.
# Uso: bash 05_ensamblar_trinity.sh
# Requisitos: Trinity 2.15.2 en el PATH (conda install -c bioconda
# trinity) o en trinity/trinityrnaseq-2.15.2/Trinity, Java,
# jellyfish (conda install -c bioconda jellyfish) para la
# normalización in silico que Trinity hace por defecto, samtools
# (conda install -c bioconda samtools), y al menos 10 GB de RAM
# libres. Sin autoconf, Trinity se compila con
# make no_bamsifter (ver Módulo 05).

mkdir -p ensamblaje

trinity/trinityrnaseq-v2.15.2/Trinity --seqType fq --left limpieza/fastp/SRR5314138_1.clean.fq.gz --right limpieza/fastp/SRR5314138_2.clean.fq.gz --CPU 4 --max_memory 10G --full_cleanup --output ensamblaje/trinity_SRR5314138

trinity/trinityrnaseq-v2.15.2/Trinity --seqType fq --left limpieza/fastp/SRR5314139_1.clean.fq.gz --right limpieza/fastp/SRR5314139_2.clean.fq.gz --CPU 4 --max_memory 10G --full_cleanup --output ensamblaje/trinity_SRR5314139

trinity/trinityrnaseq-v2.15.2/Trinity --seqType fq --left limpieza/fastp/SRR5314140_1.clean.fq.gz --right limpieza/fastp/SRR5314140_2.clean.fq.gz --CPU 4 --max_memory 10G --full_cleanup --output ensamblaje/trinity_SRR5314140

trinity/trinityrnaseq-v2.15.2/Trinity --seqType fq --left limpieza/fastp/SRR5314141_1.clean.fq.gz --right limpieza/fastp/SRR5314141_2.clean.fq.gz --CPU 4 --max_memory 10G --full_cleanup --output ensamblaje/trinity_SRR5314141

trinity/trinityrnaseq-v2.15.2/Trinity --seqType fq --left limpieza/fastp/SRR5314136_1.clean.fq.gz --right limpieza/fastp/SRR5314136_2.clean.fq.gz --CPU 4 --max_memory 10G --full_cleanup --output ensamblaje/trinity_SRR5314136

trinity/trinityrnaseq-v2.15.2/Trinity --seqType fq --left limpieza/fastp/SRR5314137_1.clean.fq.gz --right limpieza/fastp/SRR5314137_2.clean.fq.gz --CPU 4 --max_memory 10G --full_cleanup --output ensamblaje/trinity_SRR5314137

