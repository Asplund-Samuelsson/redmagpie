options(width=150)
library(tidyverse)

# Define infiles
ftpa_file = "intermediate/feature_history_correlation.archaeal_tree.tab.gz"
ftpb_file = "intermediate/feature_history_correlation.bacterial_subtrees.tab.gz"
fwil_file = "intermediate/feature_enrichment.tab"
ecim_file = "intermediate/EC_count_features.importance.tab.gz"
pfim_file = "intermediate/pfam_features.importance.tab.gz"
kgon_file = "data/kegg_enzyme.old_new.tab"
ecan_file = "data/kegg_enzyme.tab"
pfan_file = "data/Pfam-A_32.0.NAME_ACC_DESC.tab"

# Load data
ftpa = read_tsv(ftpa_file) %>% mutate(Subtree = 0)
ftpb = read_tsv(ftpb_file) %>% mutate(Domain = "Bacteria")
fwil = read_tsv(fwil_file)
pfim = read_tsv(pfim_file)
ecim = read_tsv(ecim_file)
kgon = read_tsv(kgon_file, col_names=c("OldEC", "EC"))
ecan = read_tsv(ecan_file, col_names=c("ec", "Description"))
pfan = read_tsv(pfan_file, col_names=c("Name", "Feature", "Description"))

# Calculate adjusted p-values
acco = bind_rows(ftpa, ftpb) %>%
  select(Feature, Domain, Subtree, pWilcox, R, pCor) %>%
  distinct() %>%
  group_by(is.finite(R)) %>%
  mutate(
    qWilcox = p.adjust(pWilcox, method="BH"),
    qCor = p.adjust(pCor, method="BH")
  ) %>%
  ungroup() %>%
  select(-"is.finite(R)")

# Modify data for combination
ecan = ecan %>%
  mutate(
    Feature = str_replace(ec, "ec:", "EC"),
    Name = unlist(lapply(str_split(Description, ";"), "[[", 1))
  ) %>%
  select(-ec) %>%
  mutate(
    EC = str_replace(Feature, "EC", ""),
    Description = ifelse(
      Feature == "EC3.6.3.17",
      "monosaccharide-transporting ATPase",
      Description
    ),
    Name = ifelse(
      Feature == "EC3.6.3.17",
      "monosaccharide-transporting ATPase",
      Name
    )
  ) %>%
  left_join(kgon) %>%
  mutate(
    Feature = ifelse(
      is.na(OldEC), Feature,
      paste("EC", OldEC, sep="")
    )
  ) %>%
  bind_rows(
    filter(., !is.na(OldEC)) %>%
      mutate(Feature = paste("EC", EC, sep=""))
  ) %>%
  filter(
    !grepl("Transferred to", Description)
  )

feim = bind_rows(ecim, pfim) %>%
  group_by(Feature) %>%
  summarise(
    CV = sd(Importance) / mean(Importance),
    Importance = mean(Importance)
  )

anno = bind_rows(
    pfan %>% mutate(KEGG_EC = NA),
    ecan %>% select(-OldEC) %>% rename(KEGG_EC = EC)
  )

# There are now three results tables (fwil, acco, feim) + annotations (anno)

# Sort results and add annotations
fwil = fwil %>% arrange(padj) %>% left_join(anno)
acco = acco %>% arrange(-R^2) %>% left_join(anno)
feim = feim %>% arrange(-Importance) %>% left_join(anno)

# Rename and reorder columns for consistency
fwil = fwil %>%
  rename(
    mean_Negative = Negative,
    mean_Positive = Positive,
    CV_Negative = NCV,
    CV_Positive = PCV,
    q = padj
  ) %>%
  mutate(
    Feature_Type = ifelse(startsWith(Feature, "EC"), "DeepEC", "Pfam")
  ) %>%
  select(
    Feature_Type, Feature, Name,
    mean_Negative, mean_Positive,
    CV_Negative, CV_Positive,
    p, q,
    Description, KEGG_EC
  )

acco = acco %>%
  rename(
    p_Wilcox = pWilcox,
    r = R,
    p_Correlation = pCor,
    q_Wilcox = qWilcox,
    q_Correlation = qCor
  ) %>%
  mutate(
    Feature_Type = ifelse(startsWith(Feature, "EC"), "DeepEC", "Pfam")
  ) %>%
  select(
    Feature_Type, Feature, Name,
    Domain, Subtree,
    r, p_Correlation, q_Correlation,
    p_Wilcox, q_Wilcox,
    Description, KEGG_EC
  )

feim = feim %>%
  rename(CV_Importance = CV) %>%
  mutate(
    Feature_Type = ifelse(startsWith(Feature, "EC"), "DeepEC", "Pfam")
  ) %>%
  select(
    Feature_Type, Feature, Name,
    Importance, CV_Importance,
    Description, KEGG_EC
  )

# Also create table with ACE mean (only significant)
accm = acco %>%
  filter(q_Correlation < 0.05, q_Wilcox < 0.05) %>%
  select(Feature_Type, Feature, Subtree, r) %>%
  distinct() %>%
  group_by(Feature_Type, Feature) %>%
  summarise(
    n_Subtrees = length(Subtree),
    Agreement = sum(rep((mean(r) > 0), n_Subtrees) == (r > 0)) / n_Subtrees,
    r = mean(r),
    Subtrees = paste(unique(Subtree), collapse=" ")
  ) %>%
  left_join(anno) %>%
  arrange(-Agreement, -n_Subtrees, -r^2) %>%
  select(
    Feature_Type, Feature, Name,
    r, n_Subtrees, Agreement, Subtrees,
    Description, KEGG_EC
  ) %>%
  ungroup()

# Add Rank
fwil = fwil %>% mutate(Rank = 1:nrow(.)) %>% select(Rank, everything())
acco = acco %>% mutate(Rank = 1:nrow(.)) %>% select(Rank, everything())
accm = accm %>% mutate(Rank = 1:nrow(.)) %>% select(Rank, everything())
feim = feim %>% mutate(Rank = 1:nrow(.)) %>% select(Rank, everything())

# Save tables
write_tsv(fwil, "results/Supplementary_Enrichment.tab")
write_tsv(acco, "results/Supplementary_ACE.tab")
write_tsv(accm, "results/Supplementary_ACE_mean.tab")
write_tsv(feim, "results/Supplementary_Random_forest.tab")
