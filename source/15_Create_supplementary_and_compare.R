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
    Feature_Type, Feature,
    mean_Negative, mean_Positive,
    CV_Negative, CV_Positive,
    p, q
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
    Feature_Type = ifelse(startsWith(Feature, "EC"), "DeepEC", "Pfam"),
    Significant = as.numeric(q_Wilcox < 0.001 & q_Correlation < 0.001)
  ) %>%
  select(
    Feature_Type, Feature,
    Domain, Subtree,
    r, Significant,
    p_Correlation, q_Correlation,
    p_Wilcox, q_Wilcox
  ) %>%
  filter(is.finite(r))

feim = feim %>%
  rename(CV_Importance = CV) %>%
  mutate(
    Feature_Type = ifelse(startsWith(Feature, "EC"), "DeepEC", "Pfam")
  ) %>%
  select(
    Feature_Type, Feature,
    Importance, CV_Importance
  )

# Add Rank unified for ECs and Pfams
fwil = fwil %>%
  # First rank within feature type
  group_by(Feature_Type) %>%
  mutate(Rank = rank(q)) %>%
  # Then rank based on rank and the sorting value
  ungroup() %>%
  arrange(Rank, q) %>%
  mutate(Rank = 1:length(Rank)) %>%
  arrange(Rank) %>%
  # Sort columns as wanted
  select(Rank, everything())

acca = acco %>%
  # Check number of subtrees in agreement...
  filter(Significant == 1) %>%
  group_by(Feature) %>%
  summarise(
    r_pos = sum(r >= 0),
    n = length(Subtree),
    r = mean(r),
    r2 = mean(r^2)
  ) %>%
  filter(r_pos == 0 | r_pos == n, n > 1) %>%
  arrange(-n, -r2) %>%
  select(-r_pos)

accr = acco %>%
  # Calculate the max r2 value per feature
  select(Feature_Type, Feature, Subtree, r) %>%
  group_by(Feature_Type, Feature) %>%
  summarise(max_r2 = max(r^2)) %>%
  # Create rank based on maximum r2
  mutate(Rank = rank(-max_r2)) %>%
  # Then create final rank regardless of feature type
  ungroup() %>%
  arrange(Rank, -max_r2) %>%
  mutate(Rank = 1:length(Rank)) %>%
  arrange(Rank)

acco = acco %>%
  inner_join(select(accr, Feature, Rank)) %>%
  select(Rank, everything()) %>%
  arrange(Rank, -Significant, -r^2)

feim = feim %>%
  # First rank within feature type
  group_by(Feature_Type) %>%
  mutate(Rank = rank(-Importance)) %>%
  # Then rank based on rank and the sorting value
  ungroup() %>%
  arrange(Rank, -Importance) %>%
  mutate(Rank = 1:length(Rank)) %>%
  arrange(Rank) %>%
  # Sort columns as wanted
  select(Rank, everything())

# Save tables with annotations
write_tsv(left_join(fwil, anno), "results/Supplementary_Enrichment.tab")
write_tsv(left_join(acco, anno), "results/Supplementary_ACE.tab")
write_tsv(left_join(feim, anno), "results/Supplementary_Random_forest.tab")

# Correlate methods
rnks =
  bind_rows(
    select(fwil, Feature, Rank) %>% mutate(Method="E"),
    select(accr, Feature, Rank) %>% mutate(Method="A"),
    select(feim, Feature, Rank) %>% mutate(Method="R")
  ) %>%
  spread(Method, Rank)

cors = select(rnks, -Feature) %>%
  as.matrix() %>%
  cor(use="complete.obs", method="spearman")

cors = cors[lower.tri(cors)]

corp = function(a, b){
  cor.test(a, b, use="complete.obs", method="spearman")$p.value
}

cors = tibble(
  Comparison = c("A:E", "A:R", "E:R"),
  Correlation = cors,
  p = c(
    corp(rnks$A, rnks$E),
    corp(rnks$A, rnks$R),
    corp(rnks$E, rnks$R)
  )
)

write_tsv(cors, "results/method_correlation.tab")

# Determine features that are significant in each method
sigE = filter(fwil, q < 0.05)$Feature
sigA = unique(filter(acco, Significant == 1)$Feature)
sigR = feim$Feature
sigf = tibble(
  Feature = c(sigE, sigA, sigR),
  Method = c(
    rep("E", length(sigE)), rep("A", length(sigA)), rep("R", length(sigR))
  )
)

# Plot ranks
normrank = function(x){
  nRank = -1 + rank(x)
  nRank = ifelse(is.na(x), NA, nRank)
  nRank = nRank / max(nRank, na.rm=T)
  return(abs(1 - nRank))
}

copp = rnks %>%
  mutate(A = normrank(A), E = normrank(E), R = normrank(R)) %>%
  gather(Method, Rank, -Feature) %>%
  # Keep only features that are significant
  inner_join(sigf) %>%
  # Remove features that occur in only one method
  spread(Method, Rank) %>%
  filter((is.na(A) + is.na(E) + is.na(R)) < 2) %>%
  gather(Method, Rank, -Feature) %>%
  filter(!is.na(Rank)) %>%
  mutate(
    Method = recode(
      Method,
      "E" = "Enrichment",
      "A" = "ACE",
      "R" = "Random forest"
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

gp = ggplot(copp, aes(x=Feature, color=Rank, y=Method))
gp = gp + geom_jitter(size=1)

# Add correlations
gp = gp + scale_x_discrete(expand=expand_scale(add=c(0,300)))

gp = gp + geom_segment(x=1280, y=0.75, xend=1280, yend=1.75, color="black")
cAR = round(filter(cors, Comparison == "A:R")$Correlation, 2)
gp = gp + geom_text(x=1340, y=1.25, color="black", label=cAR, size=3)

gp = gp + geom_segment(x=1280, y=2.25, xend=1280, yend=3.25, color="black")
cAE = round(filter(cors, Comparison == "A:E")$Correlation, 2)
gp = gp + geom_text(x=1340, y=2.75, color="black", label=cAE, size=3)

gp = gp + geom_segment(x=1440, y=1, xend=1440, yend=3, color="black")
cER = round(filter(cors, Comparison == "E:R")$Correlation, 2)
gp = gp + geom_text(x=1500, y=2, color="black", label=cER, size=3)

gp = gp + theme_bw()
gp = gp + scale_color_viridis_c(direction = 1)
gp = gp + theme(
  axis.title.y = element_blank(),
  axis.text.y = element_text(colour="black"),
  axis.ticks.y = element_line(colour="black"),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  panel.grid = element_blank()
)

ggsave("results/method_feature_rank_comparison.pdf", gp, h=5, w=18, units="cm")
