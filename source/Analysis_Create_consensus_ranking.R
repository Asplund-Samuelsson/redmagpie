options(width=150)
library(tidyverse)

# Load data
supp = bind_rows(
  # Load enrichment
  read_tsv("results/Supplementary_Enrichment.tab") %>%
    select(Rank, Feature_Type, Feature, mean_Negative, mean_Positive) %>%
    distinct() %>%
    mutate(Method = "E") %>%
    gather(Data, Value, -Method, -Feature_Type, -Feature),
  # Load ACE
  read_tsv("results/Supplementary_ACE.tab") %>%
  select(Feature_Type, Feature, Subtree, r, q_Correlation, q_Wilcox) %>%
  distinct() %>%
  # Calculate weighted absolute correlation
  mutate(
    X = abs(r) *
    log(q_Correlation + median(q_Correlation))/log(median(q_Correlation)) *
    log(q_Wilcox + median(q_Wilcox))/log(median(q_Correlation))
  ) %>%
  # Calculate mean weighted r per feature
  group_by(Feature_Type, Feature) %>%
  summarise(weighted_r = mean(X)) %>%
  # Add back rank
  inner_join(
    read_tsv("results/Supplementary_ACE.tab") %>%
    select(Rank, Feature) %>%
    distinct()
  ) %>%
  mutate(Method = "A") %>%
  gather(Data, Value, -Method, -Feature_Type, -Feature),
  # Load random forest
  read_tsv("results/Supplementary_Random_forest.tab") %>%
    select(Rank, Feature_Type, Feature, Importance) %>%
    distinct() %>%
    mutate(Method = "R") %>%
    gather(Data, Value, -Method, -Feature_Type, -Feature)
)

# Spread data into columns and calculate consensus rank
supp = supp %>%
  mutate(Data = ifelse(Data == "Rank", Method, Data)) %>%
  select(-Method) %>%
  spread(Data, Value) %>%
  # Set rank to one above max if NA
  mutate(
    A = ifelse(is.na(A), max(A, na.rm=T)+1, A),
    E = ifelse(is.na(E), max(E, na.rm=T)+1, E),
    R = ifelse(is.na(R), max(R, na.rm=T)+1, R)
  ) %>%
  # Calculate consensus rank
  mutate(Rank = rank(A+E+R)) %>%
  arrange(Rank) %>%
  select(-E, -A, -R)

# Load annotations
anno = distinct(bind_rows(
  read_tsv("results/Supplementary_Enrichment.tab") %>%
    select(Feature, Name, Description, KEGG_EC),
  read_tsv("results/Supplementary_ACE.tab") %>%
    select(Feature, Name, Description, KEGG_EC),
  read_tsv("results/Supplementary_Random_forest.tab") %>%
    select(Feature, Name, Description, KEGG_EC)
))

# Add back annotations and organize
supp = supp %>%
  inner_join(anno) %>%
  select(
    Rank, Feature_Type, Feature, Name,
    mean_Negative, mean_Positive,
    weighted_r,
    Importance,
    Description, KEGG_EC
  )

# Save consensus rank table
write_tsv(supp, "results/Supplementary_Consensus_rank.tab")

# Save small consensus table with top 20 ECs and top 20 Pfams
write_tsv(
  filter(supp, Feature_Type == "DeepEC") %>%
    head(20) %>%
    select(-KEGG_EC) %>%
    mutate(Description = Name) %>%
    select(-Name, -Feature_Type),
  "results/Table_Top_20_consensus_ECs.tab"
)

write_tsv(
  filter(supp, Feature_Type == "Pfam") %>%
    head(20) %>%
    select(-KEGG_EC) %>%
    select(-Name, -Feature_Type) %>%
    mutate(Feature = unlist(lapply(str_split(Feature, "\\."), "[[", 1))),
  "results/Table_Top_20_consensus_Pfams.tab"
)
