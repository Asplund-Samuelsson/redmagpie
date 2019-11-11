#!/usr/bin/env bash

# Create list of examples genomes accession IDs
(
  cut -f 1 intermediate/example_genomes.tab | tail -n +2
  cut -f 2 intermediate/example_genomes.tab | tail -n +2
) > intermediate/example_genomes.txt

# Create temporary directories for DeepEC
mkdir intermediate/orf intermediate/deepec

# Perform DeepEC analysis on the example genomes
cat intermediate/example_genomes.txt | parallel --no-notice --jobs 16 '
  gunzip -c data/orf/{}.fasta.gz > intermediate/orf/{}.fasta;
  mkdir intermediate/deepec/{};
  deepec -i intermediate/orf/{}.fasta -o intermediate/deepec/{} > /dev/null;
  tail -n +2 intermediate/deepec/{}/DeepEC_Result.txt | sed -e "s/^/{}\t/";
  rm -rf intermediate/deepec/{};
  rm intermediate/orf/{}.fasta
' > intermediate/deepec.tab 2> intermediate/deepec.error

# Zip output
pigz intermediate/deepec.tab
pigz intermediate/deepec.error

# Remove temporary directories
rmdir intermediate/orf intermediate/deepec

# Perform Pfam analysis on the example genomes
cat intermediate/example_genomes.txt | parallel --no-notice --jobs 16 '
  gunzip -c data/orf/{}.fasta.gz > intermediate/{}.fasta;
  hmmsearch --cut_tc --noali -o /dev/null \
  --tblout intermediate/{}.tblout data/Pfam-A.hmm.gz intermediate/{}.fasta;
  cat intermediate/{}.tblout | grep -v "^#" | sed -e "s/ \+/\t/g" | \
  cut -f 1,4 | sed -e "s/^/{}\t/";
  rm intermediate/{}.tblout;
  rm intermediate/{}.fasta
' > intermediate/pfam.tab 2> intermediate/pfam.error

# Zip output
pigz intermediate/pfam.tab
