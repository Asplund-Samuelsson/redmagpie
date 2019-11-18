#!/usr/bin/env bash

# Create alignment output directory
mkdir intermediate/aligned_enzymes

# Align enzymes using mafft
ls intermediate/enzymes_to_align/ | parallel --no-notice --jobs 16 '
Outfile=`echo "intermediate/aligned_enzymes/{}" | sed -e "s/.fasta/.ali.fasta/"`
cat intermediate/enzymes_to_align/{} | mafft - > ${Outfile}
'

# Create tab-delimited table with alignments
ls intermediate/aligned_enzymes |
  cut -f 1-4 -d . |
  while read EC
    do
      seqmagick convert --line-wrap 0 intermediate/aligned_enzymes/${EC}.ali.fasta - |
      cut -f 1 -d \  |
      tr -d ">" |
      paste - - |
      sed -e "s/^/${EC}\t/"
    done > intermediate/enzyme_alignment.tab

# Gzip fasta and alignment files
pigz intermediate/enzymes_to_align/*
pigz intermediate/aligned_enzymes/*
pigz intermediate/enzyme_alignment.tab

# Create Accession and Sequence ID table
paste <(
  gunzip -c intermediate/enzyme_alignment.tab.gz | cut -f 2 | cut -f 1 -d \|
  ) <(
  gunzip -c intermediate/enzyme_alignment.tab.gz | cut -f 2
  ) | sort | uniq |
  pigz > intermediate/enzyme_alignment.accession_sequence.tab.gz
