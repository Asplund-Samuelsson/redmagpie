options(width=100)
library(tidyverse)

# Load data
meta = read_tsv("data/gtdb_metadata.tab.gz")
ncya = scan("data/cyanobacteria_without_CBB.txt", character())

# Combine data
comb = meta %>%
  # Select only cyanobacteria
  filter(grepl("p__Cyanobacteria", gtdb_taxonomy)) %>%
  # Reduce to relevant data
  select(accession) %>%
  # Add CBB status
  mutate(CBB = ifelse(accession %in% ncya, 0, 1))

# Select Cyanobacteria with CBB cycle and save IDs
writeLines(filter(comb, CBB == 1)$accession, "intermediate/cyano_with_CBB.txt")
