options(width=150)
library(tidyverse)

# Define infiles
cbbx_file = "intermediate/rpp_examples.no_cyano.tab"
dpec_file = "intermediate/rpp_examples.ec.tab.gz"
meta_file = "data/gtdb_metadata.tab.gz"

# Load data
cbbx = read_tsv(cbbx_file)
dpec = read_tsv(dpec_file, col_names=c("Sequence", "EC"))
meta = read_tsv(meta_file)

# Prepare taxonomic information
taxo = meta %>%
  # Select relevant data
  select(accession, gtdb_taxonomy) %>%
  # Extract Domain, Phylum, and Class
  mutate(
    Domain = str_replace(
      unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 1)), "d__", ""
    ),
    Phylum = str_replace(
      unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 2)), "p__", ""
    ),
    Class = str_replace(
      unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 3)), "c__", ""
    )
  ) %>%
  # Drop the nasty taxonomy string and rename accession
  select(-gtdb_taxonomy) %>%
  rename(Genome = accession) %>%
  # Create Organism classification
  mutate(Organism = ifelse(Phylum == "Proteobacteria", Class, Phylum)) %>%
  select(Genome, Organism)

# Calculate average distance to relative per taxonomic group
inner_join(cbbx, taxo) %>%
  group_by(Organism) %>%
  summarise(Distance = mean(Distance), Genomes = length(Genome)) %>%
  arrange(-Distance) %>%
  as.data.frame()

# Prepare examples and features matrix
preX = dpec %>%
  # Parse Accession
  mutate(
    Accession = unlist(lapply(str_split(Sequence, fixed("|")), "[[", 1))
  ) %>%
  # Remove unwanted genomes
  filter(
    Accession %in% c(cbbx$Genome, cbbx$Relative)
  ) %>%
  # Tidy up EC
  mutate(EC = str_replace(EC, ":", "")) %>%
  # Exclude PRK and rubisco
  filter(!(EC %in% c("EC2.7.1.19", "EC4.1.1.39"))) %>%
  # Count occurrences of each EC
  group_by(Accession, EC) %>%
  summarise(Count = length(EC)) %>%
  ungroup() %>%
  # Spread ECs into columns
  spread(EC, Count)

# Create X matrix
X = preX %>% select(-Accession) %>% as.matrix()

# Set row names
rownames(X) = preX$Accession

# Set NA to 0 occurrences of the EC in question
X[is.na(X)] = 0

# Write X to file
write_csv(as_tibble(X), "intermediate/rpp_examples_X.no_cyano.csv", col_names=F)

# Create y class vector
y = ifelse(preX$Accession %in% cbbx$Genome, 1, 0)

# Write y to file
write_lines(y, "intermediate/rpp_examples_y.no_cyano.txt")

# Write feature name vector
write_lines(colnames(X), "intermediate/rpp_examples_feature_names.no_cyano.txt")

# Write data points accession ID vector
write_lines(rownames(X), "intermediate/rpp_examples_accession_ids.no_cyano.txt")
