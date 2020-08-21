options(width=150)
library(tidyverse)

# Define infiles
dpec_file = "intermediate/deepec.tab.gz"
pfam_file = "intermediate/pfam.tab.gz"
exgn_file = "intermediate/example_genomes.tab"

# Load data
dpec = read_tsv(dpec_file, col_names = c("Accession", "ORF", "EC"))
pfam = read_tsv(pfam_file, col_names = c("Accession", "ORF", "Pfam"))
exgn = read_tsv(exgn_file)


# 1. Features from enzyme counts

# Prepare examples and features matrix
preX = dpec %>%
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

# Create y class vector
y = ifelse(preX$Accession %in% exgn$Genome, 1, 0)

# Write X to file
write_csv(
  as_tibble(X),
  gzfile("intermediate/EC_count_features.X.csv.gz"),
  col_names=F
)

# Write y to file
write_lines(y, "intermediate/EC_count_features.y.txt")

# Write feature name vector
write_lines(colnames(X), "intermediate/EC_count_features.feature_names.txt")

# Write data points accession ID vector
write_lines(rownames(X), "intermediate/EC_count_features.accession_ids.txt")


# 2. Features from Pfams

# Prepare examples and features matrix
preX = pfam %>%
  # Exclude PRK and rubisco, including small subunit
  filter(
    !(Pfam %in% c("PF02788.16", "PF00016.20", "PF00485.18", "PF00101.20"))
  ) %>%
  # Count occurrences of each Pfam
  group_by(Accession, Pfam) %>%
  summarise(Count = length(Pfam)) %>%
  ungroup() %>%
  # Spread Pfams into columns
  spread(Pfam, Count)

# Create X matrix
X = preX %>% select(-Accession) %>% as.matrix()

# Set row names
rownames(X) = preX$Accession

# Set NA to 0 occurrences of the EC in question
X[is.na(X)] = 0

# Create y class vector
y = ifelse(preX$Accession %in% exgn$Genome, 1, 0)

# Write X to file
write_csv(
  as_tibble(X),
  gzfile("intermediate/pfam_features.X.csv.gz"),
  col_names=F
)

# Write y to file
write_lines(y, "intermediate/pfam_features.y.txt")

# Write feature name vector
write_lines(colnames(X), "intermediate/pfam_features.feature_names.txt")

# Write data points accession ID vector
write_lines(rownames(X), "intermediate/pfam_features.accession_ids.txt")
