#!/bin/bash
set -euxo pipefail

VERSION=$1
ASSEMBLIES=$2

########################################################################################################################
# GRCh38
########################################################################################################################

wget -O "${ASSEMBLIES}/GRCh38/GRCh38.primary_assembly.genome.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/GRCh38.primary_assembly.genome.fa.gz"

wget -O "${ASSEMBLIES}/GRCh38/gencode/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz"

wget -O "${ASSEMBLIES}/GRCh38/gencode/gencode.v${VERSION}.primary_assembly.annotation.gtf.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${VERSION}/gencode.v${VERSION}.primary_assembly.annotation.gtf.gz"

########################################################################################################################
# Postprocessing
########################################################################################################################

for file in "${ASSEMBLIES}"/*/*.fa.gz; do
  saveto="${file/.fa.gz/.fa.bgz}"
  gunzip "${file}" --stdout | bgzip --threads "$(nproc)" -l 6 -o "${saveto}" /dev/stdin
  rm "${file}"
  echo "Indexing" "${saveto}"
  samtools faidx "${saveto}"
done
