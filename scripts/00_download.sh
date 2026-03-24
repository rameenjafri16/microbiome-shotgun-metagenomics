#!/bin/bash
# 00_download.sh
# Downloads raw paired-end shotgun metagenomics FASTQ files from NCBI SRA
# Dataset: De Filippis et al. (2019), NCBI SRA accession SRP126540
# 3 vegan samples: SRR8146968, SRR8146978, SRR8146982
# 3 omnivore samples: SRR8146938, SRR8146935, SRR8146936
# Tools: SRA Toolkit v3.2.1, GNU Parallel
# Usage: bash 00_download.sh

# Download 3 samples in parallel using GNU Parallel
# prefetch: downloads SRA normalized format file
# fasterq-dump: converts to paired-end FASTQ (--split-files)
# --threads 4: use 4 threads for fasterq-dump conversion
# gzip: compress output FASTQ files to save disk space

parallel -j 3 '
  echo "Downloading {}"
  prefetch {}
  fasterq-dump --split-files --progress --threads 4 {}
  gzip {}_1.fastq {}_2.fastq
  echo "Done: {}"
' ::: SRR8146968 SRR8146978 SRR8146982 SRR8146938 SRR8146935 SRR8146936
