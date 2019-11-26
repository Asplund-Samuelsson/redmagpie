options(width=150)
library(tidyverse)
library(phytools)
library(MidpointRooter)
library(ggtree)

# Define infiles
artr_file = "intermediate/archaea.tree"
batr_file = "intermediate/bacteria_midpoint_rooted_tree.Rdata"
exgn_file = "intermediate/example_genomes.tab"
ecim_file = "intermediate/EC_count_features.importance.tab.gz"
pfim_file = "intermediate/pfam_features.importance.tab.gz"

# Load data
artr = read.tree(artr_file)
load(batr_file) # Loads "batr" object; Already midpoint-rooted bacterial tree
exgn = read_tsv(exgn_file)
ecim = read_tsv(ecim_file)
pfim = read_tsv(pfim_file)
posg = exgn$Genome

# Calculate average importances
avim = bind_rows(ecim, pfim) %>%
  group_by(Feature) %>%
  summarise(Importance = mean(Importance)) %>%
  arrange(-Importance) %>%
  mutate(Type = ifelse(str_starts(Feature, "EC"), "EC", "Pfam"))

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

# Combine features
feat = bind_rows(dpec, pfam)

# Root Archaeal tree (Bacterial tree is already rooted)
artr = midpoint.root2(artr)

# Prune trees to example genomes
artr = drop.tip(artr, setdiff(artr$tip.label, c(exgn$Genome, exgn$Relative)))
batr = drop.tip(batr, setdiff(batr$tip.label, c(exgn$Genome, exgn$Relative)))

# Create vectors with pathway status
apos = ifelse(artr$tip.label %in% posg, "Positive", "Negative")
names(apos) = artr$tip.label
bpos = ifelse(batr$tip.label %in% posg, "Positive", "Negative")
names(bpos) = batr$tip.label

# Perform Ancestral Character Estimation of genome positivity
aACE = ace(apos, artr, model="ER", type="discrete")
bACE = ace(bpos, batr, model="ER", type="discrete")


# Get subtrees of bacteria of height corresponding to the archaeal tree
btrs = treeSlice(batr, max(nodeHeights(artr)), orientation="rootwards")
# That doesn't work well

# Compare trees quick and dirty
gp = ggplot(
  bind_rows(
    mutate(ggtree(artr)$data, Domain = "Archaea"),
    mutate(ggtree(batr)$data, Domain = "Bacteria")
  )
)
gp +
  geom_tree(size=0.1) +
  theme_bw() +
  facet_grid(~Domain, scales="free_x", space="free") +
  coord_flip() +
  theme(strip.background=element_blank())

# Get subtrees of bacteria with size of the archaeal tree
btrs = getCladesofSize(batr, artr$Nnode)
# Three trees, not helpful

# Slice the tree at height = artr_height - batr_height
btrs = treeSlice(batr, max(nodeHeights(batr))-max(nodeHeights(artr)))

# Calculate distances between tips
bdis = cophenetic.phylo(batr)

bdtb = bdis %>%
  as_tibble() %>%
  mutate(Accession = colnames(.)) %>%
  gather(Relative, Distance, -Accession)

# For each tip Accession, select n closest tips, where n+1 is number of archaea
bcld = bdtb %>%
  filter(Distance != 0) %>%
  group_by(Accession) %>%
  top_n(length(artr$tip.label) - 1, -Distance)

# Check closest relatives
bcld %>%
  spread(Relative, Distance) %>%
  gather(Relative, Distance, -Accession) %>%
  mutate(Related = ifelse(is.na(Distance), 0, 1)) %>%
  select(-Distance) %>%
  spread(Relative, Related) %>%
  ungroup() %>%
  select(-Accession) %>% distinct()
# All groups of closest relatives are unique

# Use hclust
bhcl = hclust(as.dist(bdis))

# Find optimal number of clusters as close as possible to Arch. tree properties
foreach(K=1:nrow(bdis)) %dopar% {
  clst = cutree(hclust(as.dist(bdis)), k=K)
  clst = tibble(Accession = names(clst), Cluster = clst)
  cltr = midpoint.root2(
    drop.tip(
      batr, setdiff(batr$tip.label, filter(clst, Cluster == 66)$Accession)
    )
  )
}
