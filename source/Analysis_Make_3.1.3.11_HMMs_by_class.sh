#!/usr/bin/env bash

# Extract negative 3.1.3.11 sequences
paste \
intermediate/alignment_features.accession_ids.txt \
intermediate/alignment_features.y.txt |
grep "0$" | cut -f 1 | parallel --no-notice --jobs 64 '
  seqmagick convert \
  --pattern-include "^{}\|lcl" \
  --output-format fasta intermediate/aligned_enzymes/3.1.3.11.ali.fasta.gz -
' > intermediate/3.1.3.11_negative.ali.fasta

# Extract positive 3.1.3.11 sequences
paste \
intermediate/alignment_features.accession_ids.txt \
intermediate/alignment_features.y.txt |
grep "1$" | cut -f 1 | parallel --no-notice --jobs 64 '
  seqmagick convert \
  --pattern-include "^{}\|lcl" \
  --output-format fasta intermediate/aligned_enzymes/3.1.3.11.ali.fasta.gz -
' > intermediate/3.1.3.11_positive.ali.fasta

# Create HMMs
hmmbuild \
results/3.1.3.11_negative.hmm intermediate/3.1.3.11_negative.ali.fasta

hmmbuild \
results/3.1.3.11_positive.hmm intermediate/3.1.3.11_positive.ali.fasta
