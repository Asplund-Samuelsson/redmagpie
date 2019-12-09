options(width=150)
library(tidyverse)

# Define infiles
alec_file = "data/alignment_enzymes.txt"
deep_file = "intermediate/deepec.tab.gz"
exgn_file = "intermediate/example_genomes.tab"

# Load data
alec = scan(alec_file, character())
deep = read_tsv(deep_file, col_names=c("Accession", "ORF", "EC"))
exgn = read_tsv(exgn_file)

# Tidy up EC
deep = mutate(deep, EC = str_replace(EC, "EC:", ""))

# Determine ECs to keep
ecfr = deep %>%
  # Remove ORF information and make rows unique
  select(-ORF) %>%
  distinct() %>%
  # Add Class label
  inner_join(
    exgn %>%
      select(Genome, Relative) %>%
      gather(Class, Accession)
  ) %>%
  mutate(
    Class = recode(
      Class,
      "Genome" = "Positive",
      "Relative" = "Negative"
    )
  ) %>%
  # Count number of genomes per EC and Class
  group_by(EC, Class) %>%
  summarise(Genomes = length(Accession)) %>%
  # Calculate Fraction of genomes per Class that contain each EC
  ungroup() %>%
  mutate(Fraction = Genomes / nrow(exgn)) %>%
  # Select ECs found in at least 50% of Positive and Negative genomes
  select(-Genomes) %>%
  spread(Class, Fraction) %>%
  filter(Negative >= 0.5 & Positive >= 0.5 & EC %in% alec)

# Write desired DeepEC annotations to new file
write_tsv(filter(deep, EC %in% ecfr$EC), "intermediate/enzymes_to_align.tab")
