options(width=150)
library(tidyverse)
library(phytools)
library(MidpointRooter)
library(ggtree)

# Define infiles
artr_file = "intermediate/archaea.tree"
exgn_file = "intermediate/example_genomes.tab"
ecim_file = "intermediate/EC_count_features.importance.tab.gz"
pfim_file = "intermediate/pfam_features.importance.tab.gz"
taxo_file = "intermediate/accession_taxonomy.tab"

# Load data
artr = read.tree(artr_file)
exgn = read_tsv(exgn_file)
ecim = read_tsv(ecim_file)
pfim = read_tsv(pfim_file)
posg = exgn$Genome
taxo = read_tsv(taxo_file)

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

# Create vectors with pathway status
apos = ifelse(artr$tip.label %in% posg, "Positive", "Negative")
names(apos) = artr$tip.label

# Perform Ancestral Character Estimation of genome positivity
aACE = ace(apos, artr, model="ER", type="discrete")

# Perform analysis for each feature
library(foreach)
library(doMC)
registerDoMC(20)

ftpc = foreach(f=unique(feat$Feature)) %dopar% {
  # Create feature vectors
  arft = feat %>%
    filter(Feature == f) %>%
    filter(Accession %in% artr$tip.label) %>%
    arrange(match(Accession, artr$tip.label)) %>%
    pull(Count)

  names(arft) = artr$tip.label

  # Perform ACE on FBP using Brownian Motion model (default)
  arftACE = fastAnc(artr, arft)

  # Calculate correlations and p-values
  aaes = tibble(
    Likelihood = aACE$lik.anc[,2],
    Calvin = ifelse(aACE$lik.anc[,2] > 0.5, 1, 0),
    Enzyme = arftACE
  )

  aaes %>%
    group_by(Calvin) %>%
    summarise(
      Median = median(Enzyme), MAD = mad(Enzyme),
      Mean = mean(Enzyme), SD = sd(Enzyme)
    ) %>%
    mutate(
      pWilcox = wilcox.test(Enzyme~Calvin, aaes)$p.value,
      R = cor(aaes$Likelihood, aaes$Enzyme, method="spearman"),
      pCor = cor.test(aaes$Likelihood, aaes$Enzyme, method="spearman")$p.value,
      Domain = "Archaea",
      Feature = f
    ) %>% select(
    Feature, Domain, Calvin,
    Median, MAD, pWilcox,
    Mean, SD, R, pCor
  )
} %>% bind_rows()

# Save table; took 3 hours to create
write_tsv(
  ftpc, gzfile("intermediate/feature_history_correlation.archaeal_tree.tab.gz")
)

# Select only the p-value data
ftpv = ftpc %>%
  select(-Calvin, -Median, -MAD, -Mean, -SD) %>%
  distinct()

# Use correlation R to select 10 most important Features
topf = ftpv %>%
  filter(is.finite(pWilcox)) %>%
  mutate(
    padjW = p.adjust(pWilcox, method="BH"),
    padjC = p.adjust(pCor, method="BH")
  ) %>%
  filter(padjW < 0.05 & padjC < 0.05) %>%
  top_n(10, abs(R))

# Determine top Orders
min_count = 10

orgs = taxo %>%
  # Select genomes among selected Archaea
  filter(Accession %in% artr$tip.label) %>%
  # Count genomes in each Order
  group_by(Group, Order) %>%
  # Use Group instead of Order if min_count is not achieved
  mutate(
    Organism0 = ifelse(
      length(Accession) >= min_count,
      Order,
      ifelse(
        Group %in% filter(., length(Accession) >= min_count)$Group,
        paste("Other", Group),
        Group
      )
    )
  ) %>%
  # Count genomes for each Organism
  group_by(Organism0) %>%
  # Set to "Other Archaea" unless min_count is achieved
  mutate(
    Organism = ifelse(
      length(Accession) >= min_count,
      Organism0,
      "Other Archaea"
    )
  )

# Count the number of genomes in each Organism set
orgc = orgs %>%
  group_by(Organism) %>%
  summarise(Count = length(Accession)) %>%
  # Add back Group
  left_join(
    orgs %>%
      ungroup() %>%
      select(Group, Organism) %>%
      mutate(Group = ifelse(Organism == "Other Archaea", NA, Group)) %>%
      distinct()
    ) %>%
  arrange(-Count)

# Decide Colour
orgc = orgc %>%
  arrange(Group, -Count) %>%
  mutate(Colour = c("#1b7837", "#5aae61", "#a6dba0", "#d9f0d3", "#e7d4e8"))

#   Organism           Count Group         Colour
#   <chr>              <int> <chr>         <chr>
# 1 Methanomicrobiales    71 Halobacterota #1b7837
# 2 Methanosarcinales     34 Halobacterota #5aae61
# 3 Archaeoglobales       17 Halobacterota #a6dba0
# 4 Methanotrichales      13 Halobacterota #d9f0d3
# 5 Other Archaea         15 NA            #e7d4e8

# Add colour to Organism to Accession association
orgs = inner_join(orgs, select(orgc, Organism, Colour))

library(ggnewscale)

# Create dataframe with group association for heatmap
txfl = data.frame(
  row.names = orgs$Accession,
  Organism = orgs$Colour,
  stringsAsFactors=F
)

# Create dataframe with group association for heatmap
tfls = txfl[artr$tip.label, , drop=F]

# Plot trees for top Features
registerDoMC(20)

garbage = foreach(i=1:nrow(topf)) %dopar% {

  iftp = topf[i,]
  f = iftp$Feature
  D = iftp$Domain
  R = iftp$R
  pW = iftp$padjW
  pC = iftp$padjC

  # Create feature vectors
  arft = feat %>%
    filter(Feature == f) %>%
    filter(Accession %in% artr$tip.label) %>%
    arrange(match(Accession, artr$tip.label)) %>%
    pull(Count)

  names(arft) = artr$tip.label

  # Perform ACE on FBP using Brownian Motion model (default)
  arftACE = fastAnc(artr, arft)

  # Prepare extra variables
  aant = tibble(
    Calvin = aACE$lik.anc %>%
      as_tibble() %>%
      pull(Positive),
    Feature = arftACE
  )  %>%
    mutate(node = 1:nrow(.) + length(artr$tip.label)) %>%
    bind_rows(
      tibble(
        node = 1:length(artr$tip.label),
        Feature = arft,
        Calvin = ifelse(apos == "Positive", 1, 0)
      )
    )

  # Plot tree
  gp = ggtree(artr, aes(colour=Feature), layout="fan")
  gp$data = gp$data %>% inner_join(aant) # Add extra variables
  gp = gp + geom_nodepoint(mapping=aes(fill=Calvin), shape=21)
  gp = gp + geom_tippoint(mapping=aes(fill=Calvin), shape=24)
  gp = gp + scale_colour_viridis_c(option="B")
  gp = gp + scale_fill_gradient2(
    high="#5ab4ac", mid="#f5f5f5", low="#d8b365", midpoint=0.5
  )
  gp = gp + new_scale_fill()
  gp = gheatmap(
    gp, tfls,
    offset = 0.05, width = 0.05, colnames = F, color=NA
  )
  gp = gp + scale_fill_identity()

  ggsave(
    paste(
      "results/ace_archaea_", f,
      "_R", round(R, 3),
      "_Wilcoxlog10padj", round(log10(pW), 1),
      "_Corlog10padj", round(log10(pC), 1),
      ".pdf",
      sep=""
    ), gp,
    w=15, h=12, units="cm"
  )
}
