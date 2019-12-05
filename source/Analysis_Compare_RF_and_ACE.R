options(width=150)
library(tidyverse)

# Define infiles
ftpa_file = "intermediate/feature_history_correlation.archaeal_tree.tab.gz"
ftpb_file = "intermediate/feature_history_correlation.bacterial_subtrees.tab.gz"
ecim_file = "intermediate/EC_count_features.importance.tab.gz"
pfim_file = "intermediate/pfam_features.importance.tab.gz"

# Load data
ftpa = read_tsv(ftpa_file) %>% mutate(Subtree = 0)
ftpb = read_tsv(ftpb_file) %>% mutate(Domain = "Bacteria")
pfim = read_tsv(pfim_file)
ecim = read_tsv(ecim_file)

# Calculate adjusted p-values and filter ftpc file
topf = bind_rows(ftpa, ftpb) %>%
  select(Feature, Domain, Subtree, pWilcox, R, pCor) %>%
  distinct() %>%
  filter(is.finite(R)) %>%
  mutate(
    pW = p.adjust(pWilcox, method="BH"),
    pC = p.adjust(pCor, method="BH")
  ) %>%
  filter(pC < 0.001 & pW < 0.001) %>%
  arrange(-abs(R))

# Compare to random forest EC and Pfam feature importances
pfim = pfim %>%
  group_by(Feature) %>%
  summarise(
    CV = sd(Importance) / mean(Importance),
    Importance = mean(Importance)
  ) %>%
  arrange(-Importance)

ecim = ecim %>%
  group_by(Feature) %>%
  summarise(
    CV = sd(Importance) / mean(Importance),
    Importance = mean(Importance)
  ) %>%
  arrange(-Importance)

pfco = topf %>%
  filter(startsWith(Feature, "PF"))

ecco = topf %>%
  filter(startsWith(Feature, "EC"))

# Correlate mean R-squared to importance
ecco %>%
  group_by(Feature) %>%
  summarise(R2 = mean(R^2)) %>%
  inner_join(ecim) %>%
  summarise(Correlation = cor(R2, Importance, method="spearman"))

pfco %>%
  group_by(Feature) %>%
  summarise(R2 = mean(R^2)) %>%
  inner_join(pfim) %>%
  summarise(Correlation = cor(R2, Importance, method="spearman"))
