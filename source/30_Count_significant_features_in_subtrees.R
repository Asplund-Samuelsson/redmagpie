options(width=150)
library(tidyverse)

# Define infile
supA_file = "results/Supplementary_ACE.tab"

# Load data
supA = read_tsv(supA_file)

# Count number of significant features per subtree
supA %>%
  select(Feature, Subtree, Significant) %>%
  distinct() %>%
  group_by(Subtree) %>%
  summarise(
    Insignificant = sum(Significant == 0),
    Significant = sum(Significant)
  ) %>%
  mutate(Features = Significant + Insignificant, Ratio = Significant / Features)
