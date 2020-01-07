options(width=150)
library(tidyverse)

# Define infiles
ftpa_file = "intermediate/feature_history_correlation.archaeal_tree.tab.gz"
ftpb_file = "intermediate/feature_history_correlation.bacterial_subtrees.tab.gz"
fwil_file = "intermediate/feature_enrichment.tab"
ecim_file = "intermediate/EC_count_features.importance.tab.gz"
pfim_file = "intermediate/pfam_features.importance.tab.gz"

# Load data
ftpa = read_tsv(ftpa_file) %>% mutate(Subtree = 0)
ftpb = read_tsv(ftpb_file) %>% mutate(Domain = "Bacteria")
fwil = read_tsv(fwil_file)
pfim = read_tsv(pfim_file)
ecim = read_tsv(ecim_file)

# Calculate adjusted p-values and filter ancestral character estimation data
topf = bind_rows(ftpa, ftpb) %>%
  select(Feature, Domain, Subtree, pWilcox, R, pCor) %>%
  distinct() %>%
  filter(is.finite(R)) %>%
  mutate(
    pW = p.adjust(pWilcox, method="BH"),
    pC = p.adjust(pCor, method="BH")
  ) %>%
  filter(pC < 0.05 & pW < 0.05) %>%
  arrange(-abs(R))

# Filter Wilcoxon enrichment data
topw = fwil %>%
  filter(padj < 0.05) %>%
  arrange(padj)

# Compare methods
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

pfwi = topw %>%
  filter(startsWith(Feature, "PF"))

ecwi = topw %>%
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

# Perform comparison
comp = bind_rows(
  topf %>%
    group_by(Feature) %>%
    summarise(Value = mean(R^2), Data = "R2"),
  topw %>%
    mutate(Value = -log10(padj), Data = "nLog10padj") %>%
    select(Feature, Value, Data),
  bind_rows(pfim, ecim) %>%
    select(-CV) %>%
    rename(Value = Importance) %>%
    mutate(Data = "Importance")
)

# Prepare normalized rank function
normrank = function(x){
  nRank = -1 + order(x) %>% order()
  nRank = ifelse(is.na(x), NA, nRank)
  nRank = nRank / max(nRank, na.rm=T)
  return(abs(1 - nRank))
}

comp = comp %>%
  # Add data Type
  mutate(Type = ifelse(startsWith(Feature, "EC"), "EC", "Pfam")) %>%
  # For each feature Type and Data method, calculate a normalized Rank
  group_by(Data) %>%
  mutate(Rank = normrank(Value)) %>%
  # Reduce to Features that appear for more than one Data type
  filter(
    Feature %in% (
      group_by(., Feature) %>%
        summarise(Count = length(Feature)) %>%
        filter(Count > 1) %>%
        pull(Feature)
    )
  )

# Calculate correlations
corr = comp %>%
  select(-Rank) %>%
  spread(Data, Value) %>%
  group_by(Type) %>%
  summarise(
    EN_AC = cor(nLog10padj, R2, use="complete.obs", method="spearman"),
    EN_RF = cor(nLog10padj, Importance, use="complete.obs", method="spearman"),
    AC_RF = cor(R2, Importance, use="complete.obs", method="spearman"),
    pEN_AC = cor.test(nLog10padj, R2, method="spearman")$p.value,
    pEN_RF = cor.test(nLog10padj, Importance, method="spearman")$p.value,
    pAC_RF = cor.test(R2, Importance, method="spearman")$p.value
  ) %>%
  gather(Pair, R, -Type) %>%
  mutate(
    Value = ifelse(startsWith(Pair, "p"), "p", "R"),
    Pair = str_remove(Pair, "p")
  ) %>%
  spread(Value, R)

write_tsv(corr, "results/method_correlation.tab")

# Plot it
copp = comp %>%
  mutate(
    Method = recode(
      Data,
      "nLog10padj" = "Enrichment",
      "R2" = "ACE",
      "Importance" = "Random forest"
    ),
    Method = factor(Method, levels = c("Random forest", "ACE", "Enrichment")),
    Feature = factor(
      Feature,
      levels = ungroup(.) %>%
        group_by(Feature) %>%
        summarise(Rank = mean(Rank)) %>%
        arrange(Rank) %>%
        pull(Feature)
    )
  )

gp = ggplot(copp, aes(x=Feature, fill=Rank, y=Method))
gp = gp + geom_tile()
gp = gp + theme_bw()
gp = gp + facet_grid(~Type, space="free_x", scales="free_x")
gp = gp + scale_fill_viridis_c(direction = -1)
gp = gp + theme(
  strip.background = element_blank(),
  axis.text.y = element_text(colour="black"),
  axis.ticks.y = element_line(colour="black"),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

ggsave("results/method_feature_rank_comparison.pdf", gp, h=4, w=18, units="cm")
