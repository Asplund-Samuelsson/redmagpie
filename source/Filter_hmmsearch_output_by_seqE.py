#!/usr/bin/python3
import sys

# Read infile and outfile from command line
infile = sys.argv[1] # A hmmer tblout file
cutoff = float(sys.argv[2]) # Sequence E-value cutoff
outfile = sys.argv[3] # A list of sequence IDs

if infile.endswith(".gz"):
    import gzip
    infile = gzip.open(infile, 'rt')
else:
    infile = open(infile, 'r')

# Iterate over HMMER tblout file, writing sequence IDs to outfile
written = set([])
outfile = open(outfile, 'w')
for line in infile:
    if not line.startswith("#"):
        line = line.split()
        if float(line[4]) < cutoff and line[0] not in written:
            outfile.write(line[0] + "\n")
            written.add(line[0])

# Close file handles
outfile.close()
infile.close()
