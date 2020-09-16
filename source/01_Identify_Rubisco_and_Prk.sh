#!/usr/bin/env bash

# Obtain lists of UniProt IDs for Rubisco and Prk from KEGG
wget -qO - "https://www.genome.jp/dbget-bin/get_linkdb?-t+uniprot+ko:K01601" | \
grep "\[" | tr "[]" "\t" | cut -f 2 | sort | uniq > \
data/kegg_uniprot_ids.K01601_rubisco.txt

wget -qO - "https://www.genome.jp/dbget-bin/get_linkdb?-t+uniprot+ko:K00855" | \
grep "\[" | tr "[]" "\t" | cut -f 2 | sort | uniq > \
data/kegg_uniprot_ids.K00855_prk.txt

# Manually download FASTA files from UniProt based on IDs obtained from KEGG
# data/kegg_uniprot.K01601_rubisco.fasta
# data/kegg_uniprot.K00855_prk.fasta

# Cluster Rubisco and Prk sequences at 70% identity
cd-hit -d 0 -c 0.7 -i data/kegg_uniprot.K01601_rubisco.fasta \
-o intermediate/kegg_uniprot.K01601_rubisco.70.fasta

cd-hit -d 0 -c 0.7 -i data/kegg_uniprot.K00855_prk.fasta \
-o intermediate/kegg_uniprot.K00855_prk.70.fasta

# Align Rubisco and Prk sequences
mafft --thread 16 intermediate/kegg_uniprot.K01601_rubisco.70.fasta > \
intermediate/kegg_uniprot.K01601_rubisco.70.ali.fasta

mafft --thread 16 intermediate/kegg_uniprot.K00855_prk.70.fasta > \
intermediate/kegg_uniprot.K00855_prk.70.ali.fasta

# Remove positions with more than 50% gaps
seqmagick convert --squeeze-threshold 0.5 \
intermediate/kegg_uniprot.K01601_rubisco.70.ali.fasta \
intermediate/kegg_uniprot.K01601_rubisco.70.ali.squeezed.fasta

seqmagick convert --squeeze-threshold 0.5 \
intermediate/kegg_uniprot.K00855_prk.70.ali.fasta \
intermediate/kegg_uniprot.K00855_prk.70.ali.squeezed.fasta

# Check length of the alignments
seqmagick info intermediate/kegg_uniprot.*.squeezed.fasta
# Length is 401 for Rubisco and 265 for Prk

# Remove sequences that have gaps in more than 50% of the positions
seqmagick convert --min-ungapped-length 200 \
intermediate/kegg_uniprot.K01601_rubisco.70.ali.squeezed.fasta \
intermediate/kegg_uniprot.K01601_rubisco.70.ali.squeezed.purged.fasta

seqmagick convert --min-ungapped-length 132 \
intermediate/kegg_uniprot.K00855_prk.70.ali.squeezed.fasta \
intermediate/kegg_uniprot.K00855_prk.70.ali.squeezed.purged.fasta

# Remove the gaps and align the sequences again
seqmagick convert --output-format fasta --ungap \
intermediate/kegg_uniprot.K01601_rubisco.70.ali.squeezed.purged.fasta - | \
mafft - > intermediate/kegg_uniprot.K01601_rubisco.70.clean.ali.fasta

seqmagick convert --output-format fasta --ungap \
intermediate/kegg_uniprot.K00855_prk.70.ali.squeezed.purged.fasta - | \
mafft - > intermediate/kegg_uniprot.K00855_prk.70.clean.ali.fasta

# Create Rubisco and Prk HMMs from the clean alignments
hmmbuild -n K01601_rubisco data/K01601_rubisco.hmm \
intermediate/kegg_uniprot.K01601_rubisco.70.clean.ali.fasta

hmmbuild -n K00855_prk data/K00855_prk.hmm \
intermediate/kegg_uniprot.K00855_prk.70.clean.ali.fasta

# Run hmmsearch versus the original files from Uniprot to check quality
hmmsearch --noali --cpu 8 --tblout results/K01601_rubisco_self_test.tblout \
data/K01601_rubisco.hmm data/kegg_uniprot.K01601_rubisco.fasta

hmmsearch --noali --cpu 8 --tblout results/K00855_prk_self_test.tblout \
data/K00855_prk.hmm data/kegg_uniprot.K00855_prk.fasta

# Create FASTA file with all ORFs
zcat data/orf/* > intermediate/gtdb_r89_ncbi_orf.fasta

# Run hmmsearch with Rubisco and Prk HMMs on all ORFs
hmmsearch --noali --cpu 8 --tblout \
intermediate/gtdb_r89_ncbi_orf.rubisco.tblout data/K01601_rubisco.hmm \
intermediate/gtdb_r89_ncbi_orf.fasta

hmmsearch --noali --cpu 8 --tblout \
intermediate/gtdb_r89_ncbi_orf.prk.tblout data/K00855_prk.hmm \
intermediate/gtdb_r89_ncbi_orf.fasta

# Extract significant Rubisco and Prk hits at 0.01 sequence E-value cutoff
source/Filter_hmmsearch_output_by_seqE.py \
intermediate/gtdb_r89_ncbi_orf.rubisco.tblout 0.01 \
intermediate/gtdb_r89_ncbi_orf.prk.1e-2_seqids.txt

source/Filter_hmmsearch_output_by_seqE.py \
intermediate/gtdb_r89_ncbi_orf.prk.tblout 0.01 \
intermediate/gtdb_r89_ncbi_orf.rubisco.1e-2_seqids.txt

# Extract the identified Rubisco and Prk sequences from genomes
cat intermediate/gtdb_r89_ncbi_orf.rubisco.1e-2_seqids.txt | \
parallel --no-notice --jobs 32 '
  Accession=`echo {} | cut -f 1 -d \|`;
  seqmagick convert --include-from-file <(echo {}) --output-format fasta \
  data/orf/${Accession}* -
' > intermediate/gtdb_r89_ncbi_orf.rubisco.fasta

cat intermediate/gtdb_r89_ncbi_orf.prk.1e-2_seqids.txt | \
parallel --no-notice --jobs 32 '
  Accession=`echo {} | cut -f 1 -d \|`;
  seqmagick convert --include-from-file <(echo {}) --output-format fasta \
  /tmp/orf/${Accession}* -
' > intermediate/gtdb_r89_ncbi_orf.prk.fasta

# Send Rubisco and Prk FASTA sequences to KEGG BlastKOALA for KO annotation
# data/gtdb_r89_ncbi_orf.rubisco.BlastKOALA.txt
# data/gtdb_r89_ncbi_orf.prk.BlastKOALA.chunk_1.txt
# data/gtdb_r89_ncbi_orf.prk.BlastKOALA.chunk_2.txt
# data/gtdb_r89_ncbi_orf.prk.BlastKOALA.chunk_3.txt

# Run DeepEC on the putative Rubisco and Prk sequences
python /path/to/deepec/deepec.py \
-i intermediate/gtdb_r89_ncbi_orf.rubisco.fasta \
-o intermediate/deepec_rubisco/

python /path/to/deepec/deepec.py \
-i intermediate/gtdb_r89_ncbi_orf.prk.fasta \
-o intermediate/deepec_prk/

# Determine which sequences only have an incorrect EC annotation
source/Find_annotation_mismatch.py \
intermediate/deepec_rubisco/DeepEC_Result.txt EC:4.1.1.39 \
intermediate/deepec_rejected_rubisco.txt

source/Find_annotation_mismatch.py \
intermediate/deepec_prk/DeepEC_Result.txt EC:2.7.1.19 \
intermediate/deepec_rejected_prk.txt

# Filter sequence files to those that do not only have an incorrect EC
seqmagick convert --exclude-from-file intermediate/deepec_rejected_rubisco.txt \
intermediate/gtdb_r89_ncbi_orf.rubisco.fasta \
intermediate/gtdb_r89_ncbi_orf.rubisco.deepec_filtered.fasta

seqmagick convert --exclude-from-file intermediate/deepec_rejected_prk.txt \
intermediate/gtdb_r89_ncbi_orf.prk.fasta \
intermediate/gtdb_r89_ncbi_orf.prk.deepec_filtered.fasta

# Align Rubisco sequences to identify RLPs
mafft --thread 16 \
intermediate/gtdb_r89_ncbi_orf.rubisco.deepec_filtered.fasta > \
intermediate/gtdb_r89_ncbi_orf.rubisco.deepec_filtered.ali.fasta

# Determine the lysine status of position 5265 (K174 in Hanson and Tabita 2001)
seqmagick convert --cut 5265:5265 --output-format fasta \
intermediate/gtdb_r89_ncbi_orf.rubisco.deepec_filtered.ali.fasta - | \
tr -d ">" | cut -f 1 -d \  | paste - - > \
intermediate/gtdb_r89_ncbi_orf.rubisco.deepec_filtered.K_status.5265.tab

# Create list of K174 (K5265) and DeepEC filtered Rubisco
grep -P "\tK$" \
intermediate/gtdb_r89_ncbi_orf.rubisco.deepec_filtered.K_status.5265.tab | \
cut -f 1 > \
intermediate/gtdb_r89_ncbi_orf.rubisco.deepec_filtered.K_filtered.5265.txt

# Create list of DeepEC-filtered Prk sequence IDs
grep ">" intermediate/gtdb_r89_ncbi_orf.prk.deepec_filtered.fasta | \
cut -f 1 -d \  | tr -d ">" > \
intermediate/gtdb_r89_ncbi_orf.prk.deepec_filtered.txt

# Determine which sequences only have an incorrect BlastKOALA KO annotation
source/Find_annotation_mismatch.py \
<(cat intermediate/gtdb_r89_ncbi_orf.rubisco.BlastKOALA.txt | grep -P "\t") \
K01601 intermediate/kegg_rejected_rubisco.txt

source/Find_annotation_mismatch.py \
<(cat intermediate/gtdb_r89_ncbi_orf.prk.BlastKOALA.chunk_* | grep -P "\t") \
K00855 intermediate/kegg_rejected_prk.txt

# Create final accepted Rubisco and Prk sequence ID lists
diff --new-line-format="" --unchanged-line-format="" \
<(sort intermediate/gtdb_r89_ncbi_orf.rubisco.deepec_filtered.K_filtered.5265.txt) \
<(sort intermediate/kegg_rejected_rubisco.txt) > data/rubisco.txt

diff --new-line-format="" --unchanged-line-format="" \
<(sort intermediate/gtdb_r89_ncbi_orf.prk.deepec_filtered.txt) \
<(sort intermediate/kegg_rejected_prk.txt) > data/prk.txt

# Create list of positive genomes and remove Cyanobacteria
diff \
<(
  (cat data/rubisco.txt | cut -f 1 -d \| | sort | uniq;
   cat data/prk.txt | cut -f 1 -d \| | sort | uniq) | \
   sort | uniq -d | sort
) \
<(
  cut -f 1,5 intermediate/accession_taxonomy.tab | \
  grep -P "\tCyanobacteria$" | cut -f 1 | sort
) | grep "<" | cut -f 2 -d \  > data/positive_genomes.txt
