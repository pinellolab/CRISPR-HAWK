#!/bin/bash

# reproduce crisprhawk analysis 1000G (Cas9)
crisprhawk search \
    -f data/genomes/hg38/ \
    -r data/regions/cas9.bed \
    -p NGG \
    -g 20 \
    -o results/CAS9/1000G \
    -v data/variants/1000G \
    --annotation data/annotations/encode.bed.gz data/annotations/dhs.bed.gz data/annotations/cosmic.bed.gz \
    --annotation-colnames ENCODE DHS COSMIC \
    --gene-annotation data/annotations/gencode.protein_coding.bed.gz \
    --gene-annotation-colnames GENCODE \
    --haplotype-table \
    --threads 32 \
    --verbosity 2 

# reproduce crisprhawk analysis 1000G (Cpf1)
crisprhawk search \
    -f data/genomes/hg38/ \
    -r data/regions/cpf1.bed \
    -p TTTV \
    -g 23 \
    -o results/CPF1/1000G \
    -v data/variants/1000G \
    --annotation data/annotations/encode.bed.gz data/annotations/dhs.bed.gz data/annotations/cosmic.bed.gz \
    --annotation-colnames ENCODE DHS COSMIC \
    --gene-annotation data/annotations/gencode.protein_coding.bed.gz \
    --gene-annotation-colnames GENCODE \
    --right \
    --haplotype-table \
    --threads 32 \
    --verbosity 2

# reproduce crisprhawk analysis HGDP (Cas9)
crisprhawk search \
    -f data/genomes/hg38/ \
    -r data/regions/cas9.bed \
    -p NGG \
    -g 20 \
    -o results/CAS9/HGDP \
    -v data/variants/HGDP \
    --annotation data/annotations/encode.bed.gz data/annotations/dhs.bed.gz data/annotations/cosmic.bed.gz \
    --annotation-colnames ENCODE DHS COSMIC \
    --gene-annotation data/annotations/gencode.protein_coding.bed.gz \
    --gene-annotation-colnames GENCODE \
    --haplotype-table \
    --threads 32 \
    --verbosity 2

# reproduce crisprhawk analysis HGDP (Cpf1)
crisprhawk search \
    -f data/genomes/hg38/ \
    -r data/regions/cpf1.bed \
    -p TTTV \
    -g 23 \
    -o results/CPF1/HGDP \
    -v data/variants/HGDP \
    --annotation data/annotations/encode.bed.gz data/annotations/dhs.bed.gz data/annotations/cosmic.bed.gz \
    --annotation-colnames ENCODE DHS COSMIC \
    --gene-annotation data/annotations/gencode.protein_coding.bed.gz \
    --gene-annotation-colnames GENCODE \
    --right \
    --haplotype-table \
    --threads 32 \
    --verbosity 2

# reproduce crisprhawk analysis GNOMAD (Cas9)
crisprhawk search \
    -f data/genomes/hg38/ \
    -r data/regions/cas9.bed \
    -p NGG \
    -g 20 \
    -o results/CAS9/GNOMAD \
    -v data/variants/GNOMAD \
    --annotation data/annotations/encode.bed.gz data/annotations/dhs.bed.gz data/annotations/cosmic.bed.gz \
    --annotation-colnames ENCODE DHS COSMIC \
    --gene-annotation data/annotations/gencode.protein_coding.bed.gz \
    --gene-annotation-colnames GENCODE \
    --haplotype-table \
    --threads 32 \
    --verbosity 2

# reproduce crisprhawk analysis GNOMAD (Cpf1)
crisprhawk search \
    -f data/genomes/hg38/ \
    -r data/regions/cpf1.bed \
    -p TTTV \
    -g 23 \
    -o results/CPF1/GNOMAD \
    -v data/variants/GNOMAD \
    --annotation data/annotations/encode.bed.gz data/annotations/dhs.bed.gz data/annotations/cosmic.bed.gz \
    --annotation-colnames ENCODE DHS COSMIC \
    --gene-annotation data/annotations/gencode.protein_coding.bed.gz \
    --gene-annotation-colnames GENCODE \
    --right \
    --haplotype-table \
    --threads 32 \
    --verbosity 2

