options(width=150)
library(tidyverse)

# Load data
meta = read_tsv("data/gtdb_metadata.tab.gz")
ncya = scan("data/cyanobacteria_without_CBB.txt", character())

# Combine data
comb = meta %>%
  # Select only cyanobacteria
  filter(grepl("p__Cyanobacteria", gtdb_taxonomy)) %>%
  # Reduce to relevant data
  select(accession, checkm_completeness) %>%
  rename(Accession = accession, Completeness = checkm_completeness) %>%
  # Add CBB status
  mutate(CBB = ifelse(Accession %in% ncya, 0, 1))

# Test difference
wilcox.test(Completeness~CBB, comb)$p.value
# [1] 5.849147e-33

# Calculate median completeness
comb %>% group_by(CBB) %>% summarise(Completeness = median(Completeness))
# # A tibble: 2 x 2
#     CBB Completeness
#   <dbl>        <dbl>
# 1     0         82.8
# 2     1         98.8

# Check how many complete genomes are CBB-negative
comb %>%
  mutate(Completeness = ifelse(Completeness < 99, "<99", "≥99")) %>%
  group_by(Completeness) %>% summarise(CBB = sum(CBB) / length(CBB))
# # A tibble: 2 x 2
#   Completeness   CBB
#   <chr>        <dbl>
# 1 <99          0.595
# 2 ≥99          0.937
