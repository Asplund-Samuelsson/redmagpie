#!/usr/bin/python3
import os, gzip

# Read input arguments
orfdir = "data/orf"
table = "intermediate/enzymes_to_align.tab"
outdir = "intermediate/enzymes_to_align"

# Create output directory
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Load table of ECs to extract and the corresponding Accession and ORF IDs
accession_orf_ec = {}
opened_ecs = {}

for line in open(table):
    # Skip first line
    if line == "Accession\tORF\tEC\n":
        continue
    # Split line into Accession, ORF, and EC
    accession, orf, ec = line.rstrip().split('\t')
    # Open EC output fasta if not open already
    if ec not in opened_ecs:
        opened_ecs[ec] = open(os.path.join(outdir, ec + ".fasta"), 'w')
    # Add Accession, ORF, and EC to dictionary
    try:
        # If Accession and ORF has been added, append EC
        accession_orf_ec[accession][orf].append(ec)
    except KeyError:
        try:
            # If ORF has not been added to Accession, create new list with EC
            accession_orf_ec[accession][orf] = [ec]
        except KeyError:
            # If Accession has not been added, create new dictionary with ORF+EC
            accession_orf_ec[accession] = {orf : [ec]}

# Iterate over Accessions, writing ORFs to the correct EC outfiles
for accession in accession_orf_ec:
    # Open fasta for the Accession
    for line in gzip.open(os.path.join(orfdir, accession + '.fasta.gz'), 'rt'):
        if line.startswith(">"):
            seqid = line.lstrip(">").split(" ")[0]
            line = ">" + accession + "|" + seqid + "\n"
        if seqid in accession_orf_ec[accession]:
            for ec in accession_orf_ec[accession][seqid]:
                junk = opened_ecs[ec].write(line)

# Close all the file handles
for ec in opened_ecs:
    opened_ecs[ec].close()
