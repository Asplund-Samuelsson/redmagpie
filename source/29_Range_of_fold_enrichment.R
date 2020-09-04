options(width=150)
library(tidyverse)

# Load data
enri = read_tsv("results/Supplementary_Enrichment.tab")

# Check distribution of ratios
filter(enri, q < 0.05) %>%
  select(Feature, mean_Negative, mean_Positive) %>%
  distinct() %>%
  filter(mean_Positive > 0, mean_Negative > 0) %>%
  mutate(Ratio = mean_Positive / mean_Negative) %>%
  pull(Ratio) %>%
  log2() %>%
  hist(50)

# Check range of ratios
filter(enri, q < 0.05) %>%
  select(Feature, mean_Negative, mean_Positive) %>%
  distinct() %>%
  filter(mean_Positive > 0, mean_Negative > 0) %>%
  mutate(Ratio = mean_Positive / mean_Negative) %>%
  pull(Ratio) %>%
  log2() %>%
  range()

# Check range of 95% middle values
filter(enri, q < 0.05) %>%
  select(Feature, mean_Negative, mean_Positive) %>%
  distinct() %>%
  filter(mean_Positive > 0, mean_Negative > 0) %>%
  mutate(Ratio = mean_Positive / mean_Negative) %>%
  arrange(Ratio) %>%
  slice(ceiling(0.025*nrow(.)):floor(0.975*nrow(.))) %>%
  pull(Ratio) %>%
  log2() %>%
  range()

# Find smallest significant ratio
filter(enri, q < 0.05) %>%
  select(Feature, mean_Negative, mean_Positive) %>%
  distinct() %>%
  filter(mean_Positive > 0, mean_Negative > 0) %>%
  mutate(
    Ratio = mean_Positive / mean_Negative,
    Ratio = log2(Ratio),
    absRatio = abs(Ratio)
  ) %>%
  arrange(absRatio)
