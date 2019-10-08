options(width=150)
library(tidyverse)

# Define infiles
meta_file = "data/gtdb_metadata.tab.gz"

# Load data
meta = read_tsv(meta_file)

# Select relevant data
meta = meta %>% select(accession, gtdb_taxonomy)

# Get phylogenetic group
meta = meta %>%
  mutate(
    Phylum = unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 2)),
    Class = unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 3)),
    Order = unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 4))
  ) %>%
  select(-gtdb_taxonomy) %>%
  mutate(
    Phylum = str_replace(Phylum, "p__", ""),
    Class = str_replace(Class, "c__", ""),
    Order = str_replace(Order, "o__", "")
  )

meta = mutate(meta, Group = ifelse(Phylum == "Proteobacteria", Class, Phylum))
meta = mutate(meta, Group = unlist(lapply(str_split(Group, "_"), "[[", 1)))

# Save final table
write_tsv(
  rename(meta, Accession = accession),
  "intermediate/accession_taxonomy.tab"
)
