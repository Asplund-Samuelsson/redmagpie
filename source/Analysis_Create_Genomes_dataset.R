options(width=150)
library(tidyverse)

# Define infiles
meta_file = "data/gtdb_metadata.tab.gz"
exgn_file = "intermediate/example_genomes.tab"

# Load data
meta = read_tsv(meta_file)
exgn = read_tsv(exgn_file)

# Load features
dpec = read_csv(
  "intermediate/EC_count_features.X.csv.gz",
  col_names = scan(
    "intermediate/EC_count_features.feature_names.txt", character()
    )
  ) %>%
  # Add Accession IDs
  mutate(
    Accession = scan(
      "intermediate/EC_count_features.accession_ids.txt", character()
    )
  ) %>%
  # Gather into long format
  gather(Feature, Count, -Accession) %>%
  # Store origin of Data
  mutate(Data = "EC")

pfam = read_csv(
    "intermediate/pfam_features.X.csv.gz",
    col_names = scan(
      "intermediate/pfam_features.feature_names.txt", character()
    )
  ) %>%
  # Add Accession IDs
  mutate(
    Accession = scan("intermediate/pfam_features.accession_ids.txt", character())
  ) %>%
  # Gather into long format
  gather(Feature, Count, -Accession) %>%
  mutate(Data = "Pfam")

# Combine data
dset = bind_rows(
    exgn %>%
      rename(Accession = Genome) %>%
      mutate(CBB_status = "CBB-positive"),
    exgn %>%
      rename(Accession = Relative, Relative = Genome) %>%
      mutate(CBB_status = "CBB-negative")
  ) %>%
  inner_join(
    meta %>%
      select(accession, checkm_completeness, gtdb_taxonomy) %>%
      rename(Accession = accession)
  ) %>%
  inner_join(dpec %>% select(-Data) %>% spread(Feature, Count)) %>%
  inner_join(pfam %>% select(-Data) %>% spread(Feature, Count))

dset = dset %>% select(-Domain)

# Save supplementary dataset
write_tsv(dset, "results/Supplementary_Example_genomes.tab")
system("gzip results/Supplementary_Example_genomes.tab")
