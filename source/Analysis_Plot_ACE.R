#!/usr/bin/env Rscript

### COMMAND LINE ARGUMENTS #####################################################

library(optparse)

# Parse command line arguments
option_list = list(
  make_option(
    c("-t", "--tree"), type="integer", default=0,
    help="Bacterial subtree (0 for Archaea) to plot."
  ),
  make_option(
    c("-f", "--feature"), type="character", default="",
    help="Feature to plot."
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

Tree = opt$tree
f = opt$feature

# For testing
# Tree = 7
# f = "PF08406.10"

################################################################################

options(width=150)
library(tidyverse)
library(phytools)
library(MidpointRooter)
library(ggtree)
library(ggnewscale)

# Define infiles
artr_file = "intermediate/archaea.tree"
batr_file = "intermediate/bacteria_midpoint_rooted_tree.Rdata"
exgn_file = "intermediate/example_genomes.tab"
orgs_file = "intermediate/accession_organism_colours.tab"
acea_file = "intermediate/feature_history_correlation.archaeal_tree.tab.gz"
aceb_file = "intermediate/feature_history_correlation.bacterial_subtrees.tab.gz"
subt_file = "intermediate/ace_bacterial_subtrees.tab"

# Load data
artr = read.tree(artr_file)
load(batr_file)
exgn = read_tsv(exgn_file)
posg = exgn$Genome
orgs = read_tsv(orgs_file)
acea = read_tsv(acea_file)
aceb = read_tsv(aceb_file)
subt = read_tsv(subt_file)

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

# Relabel bacterial tree nodes with actual node ID in full tree
batr$node.label = 1:batr$Nnode + length(batr$tip.label)

# Check what kind of Subtrees are in there
subs = subt %>%
  group_by(Subtree) %>%
  summarise(
    Count = length(Subtree),
    Calvin = sum(Accession %in% posg)/length(Subtree)
  )

# Get rid of trees smaller than 50 genomes
subs = filter(subs, Count >= 50)
subt = filter(subt, Subtree %in% subs$Subtree)

# Select only the p-value data
ftpc = bind_rows(acea, aceb) %>%
  mutate(Subtree = ifelse(is.na(Subtree), 0, Subtree)) %>%
  filter(is.finite(pWilcox)) %>%
  mutate(
    padjW = p.adjust(pWilcox, method="BH"),
    padjC = p.adjust(pCor, method="BH")
  )

ftpv = ftpc %>%
  select(-Calvin, -Median, -MAD, -Mean, -SD) %>%
  distinct()

# Create dataframe with group association for heatmap
txfl = data.frame(
  row.names = orgs$Accession,
  Organism = orgs$Colour,
  stringsAsFactors=F
)

# Create dataframe with group association for heatmap
tfla = txfl[artr$tip.label, , drop=F]
tflb = txfl[batr$tip.label, , drop=F]

if (Tree > 0){
  # Extract subtree
  sbtr = drop.tip(
    batr,
    setdiff(batr$tip.label, filter(subt, Subtree == Tree)$Accession)
  )
  # Perform Calvin ACE on subtree
  cpos = ifelse(sbtr$tip.label %in% posg, "Positive", "Negative")
  names(cpos) = sbtr$tip.label
  cACE = ace(cpos, sbtr, model="ER", type="discrete")
}

# Create feature vectors
arft = feat %>%
  filter(Feature == f) %>%
  filter(Accession %in% artr$tip.label) %>%
  arrange(match(Accession, artr$tip.label)) %>%
  pull(Count)

names(arft) = artr$tip.label

if (Tree > 0){
  baft = feat %>%
    filter(Feature == f) %>%
    filter(Accession %in% sbtr$tip.label) %>%
    arrange(match(Accession, sbtr$tip.label)) %>%
    pull(Count)

  names(baft) = sbtr$tip.label
}

# Perform ACE using Brownian Motion model (default)
arftACE = fastAnc(artr, arft)
if(Tree>0){baftACE = fastAnc(sbtr, baft)}

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

if(Tree > 0){
  bant = tibble(
    Calvin = cACE$lik.anc %>%
      as_tibble() %>%
      pull(Positive),
    Feature = baftACE
  ) %>%
    mutate(node = 1:nrow(.) + length(sbtr$tip.label)) %>%
    bind_rows(
      tibble(
        node = 1:length(sbtr$tip.label),
        Feature = baft,
        Calvin = ifelse(cpos == "Positive", 1, 0)
      )
    )
}

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
  gp, tfla,
  offset = 0.05, width = 0.05, colnames = F, color=NA
)
gp = gp + scale_fill_identity()
gp = gp + geom_treescale(x=0, y=0, fontsize=3)


iftp = filter(ftpv, Feature == f & Subtree == 0)
R = iftp$R
pW = iftp$pWilcox
pC = iftp$padjC

if (Tree == 0){
  ggsave(
    paste(
      "results/ace_archaea_", f,
      "_R", round(R, 3),
      "_Wilcoxlog10padj", round(log10(pW), 1),
      "_Corlog10padj", round(log10(pC), 1),
      ".on-demand.pdf",
      sep=""
    ), gp,
    w=15, h=12, units="cm"
  )
}

# Plot tree
if (Tree > 0){

  gp = ggtree(sbtr, aes(colour=Feature), layout="fan")
  gp$data = gp$data %>% inner_join(bant) # Add extra variables
  gp = gp + geom_nodepoint(mapping=aes(fill=Calvin), shape=21)
  gp = gp + geom_tippoint(mapping=aes(fill=Calvin), shape=24)
  gp = gp + scale_colour_viridis_c(option="B")
  gp = gp + scale_fill_gradient2(
    high="#5ab4ac", mid="#f5f5f5", low="#d8b365", midpoint=0.5
  )
  gp = gp + new_scale_fill()
  gp = gheatmap(
    gp, tflb,
    offset = 0.05, width = 0.05, colnames = F,
    color = NA
  )
  gp = gp + scale_fill_identity()
  gp = gp + geom_treescale(x=0, y=0, fontsize=3)

  sclf = max(c(sqrt(length(sbtr$tip.label) / length(artr$tip.label)), 1))

  iftp = filter(ftpv, Feature == f & Subtree == Tree)
  R = iftp$R
  pW = iftp$padjW
  pC = iftp$padjC

  ggsave(
    paste(
      "results/ace_bacteria_subtrees_",
      Tree, "_", f,
      "_R", round(R, 3),
      "_Wilcoxlog10padj", round(log10(pW), 1),
      "_Corlog10padj", round(log10(pC), 1),
      ".on-demand.pdf",
      sep=""
    ),
    gp, w=sclf*15, h=sclf*12, units="cm", limitsize=F
  )

}
