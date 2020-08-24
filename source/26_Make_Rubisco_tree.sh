#!/usr/bin/env bash

# Create alignment of Tabita examples and GTDB Rubiscos
mafft --thread 24 \
<(cat data/2007_Tabita_Rubisco_examples.fasta  data/rubisco.fasta) > \
intermediate/rubisco_with_Tabita_examples.ali.fasta

# Make tree of Tabita exampels and GTDB Rubiscos
FastTreeMP intermediate/rubisco_with_Tabita_examples.ali.fasta > \
intermediate/rubisco_with_Tabita_examples.tree
