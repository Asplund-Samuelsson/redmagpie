options(width=150)
library(tidyverse)

# Load metadata
exgn = read_tsv("intermediate/example_genomes.tab")
comp = read_tsv("data/gtdb_metadata.tab.gz") %>%
  select(accession, checkm_completeness) %>%
  rename(Accession = accession, Completeness = checkm_completeness) %>%
  filter(Accession %in% c(exgn$Genome, exgn$Relative)) %>%
  mutate(Completeness = Completeness / 100)

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

# Consider only genomes with at least 95% completeness
feat = comp %>%
  filter(Completeness >= 0.95) %>%
  select(Accession) %>%
  mutate(CBB = ifelse(Accession %in% exgn$Genome, "Positive", "Negative")) %>%
  inner_join(bind_rows(dpec, pfam))

fwil = feat %>%
  group_by(CBB, Feature) %>%
  summarise(Count = list(Count)) %>%
  spread(CBB, Count) %>%
  group_by(Feature) %>%
  mutate(
    p = wilcox.test(unlist(Negative), unlist(Positive))$p.value,
    NCV = sd(unlist(Negative)) / mean(unlist(Negative)),
    Negative = mean(unlist(Negative)),
    PCV = sd(unlist(Positive)) / mean(unlist(Positive)),
    Positive = mean(unlist(Positive))
  ) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p, method="BH"))

# Save enrichment table
write_tsv(fwil, "intermediate/feature_enrichment.tab")
