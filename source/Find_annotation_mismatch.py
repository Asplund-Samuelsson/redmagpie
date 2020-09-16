#!/usr/bin/python3
import sys

# Read infile and outfile from command line
infile = sys.argv[1] # A DeepEC output file
expected = sys.argv[2] # Expected EC (e.g. "EC:4.1.1.39")
outfile = sys.argv[3] # A list of rejected sequence IDs

# Store all and also accepted sequence IDs
seqids = set()
accepted = set()

# Open infile
infile = open(infile, 'r')

# Iterate over lines of DeepEC file
for line in infile:
    # Skip first line
    if not line == "Query ID\tPredicted EC number\n":
        # Split line into two column entries
        line = line.strip().split()
        # Add all sequence IDs to first set
        seqids.add(line[0])
        # Add expected EC sequence IDs to second set
        if expected == line[1]:
            accepted.add(line[0])

# Close infile
infile.close()

# Open outfile
outfile = open(outfile, 'w')

# Write rejected sequence IDs to outfile
for seqid in seqids:
    if seqid not in accepted:
        junk = outfile.write(seqid + "\n")

# Close outfile
outfile.close()
