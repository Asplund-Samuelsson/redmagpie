options(width=150)
library(tidyverse)
library(phytools)
library(MidpointRooter)
library(ggtree)

# Define infiles
artr_file = "intermediate/archaea.tree"
batr_file = "intermediate/bacteria_midpoint_rooted_tree.Rdata"
# batr_file = "intermediate/bacteria.tree"
exgn_file = "intermediate/example_genomes.tab"
orgs_file = "intermediate/accession_organism_colours.tab"

# Load data
artr = read.tree(artr_file)
load(batr_file) # Loads "batr" object; Already midpoint-rooted bacterial tree
# batr = read.tree(batr_file)
exgn = read_tsv(exgn_file)
posg = exgn$Genome
orgs = read_tsv(orgs_file)

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

# Root trees
artr = midpoint.root2(artr)
# batr = midpoint.root2(batr) # Slow! Alternatively, load already created object

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

# Prepare extra variables
bant = tibble(
  Calvin = bACE$lik.anc %>%
    as_tibble() %>%
    pull(Positive)
) %>%
  mutate(node = 1:nrow(.) + length(batr$tip.label)) %>%
  bind_rows(
    tibble(
      node = 1:length(batr$tip.label),
      Calvin = ifelse(bpos == "Positive", 1, 0)
    )
  ) %>%
  mutate(Domain="Bacteria")

# Prepare extra variables
aant = tibble(
  Calvin = aACE$lik.anc %>%
    as_tibble() %>%
    pull(Positive)
)  %>%
  mutate(node = 1:nrow(.) + length(artr$tip.label)) %>%
  bind_rows(
    tibble(
      node = 1:length(artr$tip.label),
      Calvin = ifelse(apos == "Positive", 1, 0)
    )
  ) %>%
  mutate(Domain="Archaea")


# Compare trees quick and dirty
gp = ggplot(
  bind_rows(
    mutate(ggtree(artr)$data, Domain = "Archaea"),
    mutate(ggtree(batr)$data, Domain = "Bacteria")
  )
)
gp$data = gp$data %>% inner_join(bind_rows(bant, aant)) # Add extra variables
gp = gp +
  geom_tree(size=0.1, mapping=aes(colour=Calvin)) +
  theme_bw() +
  facet_grid(~Domain, scales="free_x", space="free") +
  coord_flip() +
  theme(strip.background=element_blank()) +
  scale_colour_gradient2(
    high="#5ab4ac", mid="#888888", low="#d8b365", midpoint=0.5
  )

ggsave("results/arc_bac_tree_comparison.pdf", gp, h=15/2.54, w=40/2.54)

# Count Positive and Negative among Archaea
acal = tibble(
    Calvin = ifelse(artr$tip.label %in% posg, 1, 0)
  ) %>%
  group_by(Calvin) %>%
  summarise(Reference = length(Calvin))

# Calculate height of Archaeal tree
ahgt = max(nodeHeights(artr))

# Calculate coefficient of variation of Archaeal edge length
aecv = sd(artr$edge.length) / mean(artr$edge.length)

# Load parallel processing libraries
library(foreach)
library(doMC)
registerDoMC(16)

# Function to calculate a rough similarity between a tree and the archaeal tree
archsim = function(sbtr){
  suppressMessages(
    tibble(
      Calvin = ifelse(sbtr$tip.label %in% posg, 1, 0)
    ) %>%
    group_by(Calvin) %>%
    summarise(Count = length(Calvin)) %>%
    inner_join(acal) %>%
    # Compare count of positive and negative genomes to Archaea
    mutate(
      Factor = exp(-abs(log(Count / Reference))),
      Calvin = ifelse(Calvin == 1, "Positive", "Negative")
    ) %>%
    rename(Aspect = Calvin) %>%
    select(Aspect, Factor) %>%
    # Compare height and CV to Archaea
    bind_rows(
      tibble(
        Aspect = c("Height", "CV"),
        Factor = c(
          exp(-abs(log(max(nodeHeights(sbtr))/ahgt))),
          exp(-abs(log((sd(sbtr$edge.length)/mean(sbtr$edge.length)/aecv))))
        )
      )
    )
  )
}

# Relabel bacterial tree nodes with actual node ID in full tree
batr$node.label = 1:batr$Nnode + length(batr$tip.label)

# Precompute Subtree Similarity for each node
sims = bind_rows(foreach(n=1:batr$Nnode + length(batr$tip.label)) %dopar% {
  TREE = extract.clade(batr, n)
  tibble(
    Node = n,
    # Calculate the four Factors and multiply them
    Similarity = archsim(TREE) %>% pull(Factor) %>% prod(),
    Size = TREE$tip.label %>% length()
  )
})

# Calculate midpoints between heights (to accomodate treeSlice behaviour)
mpnt = tibble(High = sort(unique(c(nodeHeights(batr))))) %>%
  mutate(
    Low = lag(High, 1),
    Height = (High + Low)/2,
    Offset = Height - Low
  ) %>%
  filter(!is.na(Height))

# Precompute tree slices
slcs = bind_rows(foreach(h=mpnt$Height) %dopar% {
  # Determine the subtrees
  tryCatch(
    {
      tibble(
        Height = h,
        Node = unlist(lapply(
          treeSlice(batr, h),
          function(sbtr){sbtr$node.label[1]}
        ))
      )
    },
    # If there is no subtree, do nothing
    error = function(e) {}
  )
})

# Add the root node
slcs = bind_rows(tibble(Height = 0, Node = batr$node.label[1]), slcs)

# Combine slices and similarities
slsm = inner_join(slcs, sims) %>%
  # Filter slices to between 50% and 200% of Archaeal tree size
  filter(Size >= 50, Size <= 300) %>%
  # Height is irrelevant
  select(-Height) %>%
  distinct()

# Determine all nodes in all subtrees
nods = bind_rows(lapply(
  1:nrow(slsm),
  function(n){
    tibble(
      Subtree = n,
      Node = slsm$Node[n],
      Members = getDescendants(batr, slsm$Node[n])
    )
  }
))

# Create table with subtree information
sbts = nods %>%
  select(Subtree, Node) %>%
  distinct() %>%
  inner_join(slsm)

# Function for determining valid subtree companions of set of subtrees
valid_subtrees = function(s){
  # Determine invalid subtrees
  ival = filter(nods, Members %in% filter(nods, Subtree %in% s)$Members)$Subtree
  # Determine valid subtrees
  filter(nods, !(Subtree %in% ival))$Subtree %>% unique()
}

# Find optimal subtrees by iteratively picking the best valid tree
optimal_subtrees = function(s=c()) {
  # Find valid subtree partners
  vals = valid_subtrees(s)
  # If there are no valid partners...
  if (length(vals) == 0) {
    # ...return subtrees
    return(c(s))
  } else {
    # ...otherwise pick the best subtree
    b = sbts %>%
      # The best subtree must be valid
      filter(Subtree %in% vals) %>%
      # Take the most similar subtree
      slice_max(Similarity) %>%
      pull(Subtree)
    # Continue finding valid subtrees
    optimal_subtrees(c(s, b))
  }
}

# Get the list of optimal subtrees
subo = optimal_subtrees()

# Create optimal subtree table
subt = bind_rows(lapply(
  subo, function(i){
    tibble(
      Subtree = i,
      Accession = extract.clade(batr, filter(sbts, Subtree == i)$Node)$tip.label
    )
  }
))

# Rename subtrees in order of Size
subt = subt %>%
  group_by(Subtree) %>%
  summarise(Size = length(Subtree)) %>%
  arrange(-Size) %>%
  mutate(SubtreeX = 1:length(Subtree)) %>%
  select(-Size) %>%
  inner_join(subt) %>%
  select(-Subtree) %>%
  rename(Subtree = SubtreeX)

# Check what kind of Subtrees are in there
subs = subt %>%
  group_by(Subtree) %>%
  summarise(
    Count = length(Subtree),
    Calvin = sum(Accession %in% posg)/length(Subtree)
  )

# Save subtree table
write_tsv(subt, "intermediate/ace_bacterial_subtrees.tab")

library(ggnewscale)

# Create dataframe with group association for heatmap
txfl = data.frame(
  row.names = orgs$Accession,
  Organism = orgs$Colour,
  stringsAsFactors=F
)
tfls = txfl[batr$tip.label, , drop=F]

# Prepare extra variables
bant = tibble(
  Calvin = bACE$lik.anc %>%
    as_tibble() %>%
    pull(Positive)
) %>%
  mutate(node = 1:nrow(.) + length(batr$tip.label)) %>%
  bind_rows(
    tibble(
      node = 1:length(batr$tip.label),
      Calvin = ifelse(bpos == "Positive", 1, 0)
    )
  )

# Calculate last common ancestors of subtrees
clrs = c(
  "#bf812d", "#f6e8c3", "#9970ab", "#e7d4e8",
  "#35978f", "#c7eae5", "#1b7837", "#a6dba0",
  "#8c510a", "#dfc27d", "#762a83", "#c2a5cf",
  "#01665e", "#80cdc1", "#5aae61", "#d9f0d3"
)

mrca = bind_rows(lapply(
  unique(subt$Subtree),
  function(x){
    tibble(
      Subtree=x,
      node=findMRCA(batr, filter(subt, Subtree == x)$Accession)
    )
  }
)) %>%
  mutate(Colour = clrs[1:nrow(.)])

# Check the distribution of the Subtrees on the phylogenetic tree
gp = ggtree(batr, layout="fan")
gp$data = left_join(
  gp$data,
  subt %>% mutate(Subtree = as.character(Subtree)) %>% rename(label = Accession)
) %>% inner_join(bant)
for (n in mrca$node){
  gp = gp + geom_hilight(node=n)
  gp = gp + geom_cladelabel(
    node=n, label=filter(mrca, node == n)$Subtree, offset.text=0.05
  )
}
gp = gp + scale_fill_manual(

)
gp = gp + new_scale_fill()
gp = gp + geom_nodepoint(mapping=aes(fill=Calvin), shape=21)
gp = gp + geom_tippoint(mapping=aes(fill=Calvin), shape=24)
gp = gp + scale_colour_viridis_c(option="B")
gp = gp + scale_fill_gradient2(
  high="#5ab4ac", mid="#f5f5f5", low="#d8b365", midpoint=0.5
)
gp = gp + new_scale_fill()
gp = gheatmap(
  gp, tfls,
  offset=.05, width=.05, colnames=F,
  color = NA
)
gp = gp + scale_fill_identity()
gp = gp + geom_treescale(x=1.9, y=0, fontsize=3)

ggsave("results/ace_bacterial_subtrees.orgs.pdf", gp, w=30, h=30, units="cm")

# Perform correlation and testing on Subtrees
ftpc = bind_rows(lapply(subs$Subtree, function(k){
  # Extract subtree
  sbtr = drop.tip(
    batr,
    setdiff(batr$tip.label, filter(subt, Subtree == k)$Accession)
  )

  # Perform Calvin ACE on subtree
  cpos = ifelse(sbtr$tip.label %in% posg, "Positive", "Negative")
  names(cpos) = sbtr$tip.label
  cACE = ace(cpos, sbtr, model="ER", type="discrete")

  # Analyze
  foreach(f=unique(feat$Feature)) %dopar% {
    # Create feature vector
    baft = feat %>%
      filter(Feature == f) %>%
      filter(Accession %in% sbtr$tip.label) %>%
      arrange(match(Accession, sbtr$tip.label)) %>%
      pull(Count)

    names(baft) = sbtr$tip.label

    # Perform ACE on Feature using Brownian Motion model (default)
    baftACE = fastAnc(sbtr, baft)

    # Calculate correlations and p-values
    baes = tibble(
      Likelihood = cACE$lik.anc[,2],
      Calvin = ifelse(cACE$lik.anc[,2] > 0.5, 1, 0),
      Enzyme = baftACE
    )
    baes %>%
      group_by(Calvin) %>%
      summarise(
        Median = median(Enzyme), MAD = mad(Enzyme),
        Mean = mean(Enzyme), SD = sd(Enzyme)
      ) %>%
      mutate(
        pWilcox = wilcox.test(Enzyme~Calvin, baes)$p.value,
        R = cor(baes$Likelihood, baes$Enzyme, method="spear"),
        pCor = cor.test(baes$Likelihood, baes$Enzyme, method="spear")$p.value,
        Subtree = k,
        Feature = f
      ) %>%
      select(
        Feature, Subtree, Calvin,
        Median, MAD, pWilcox,
        Mean, SD, R, pCor
      )
    } %>% bind_rows()
  }
))

# Save table; takes some time to generate
write_tsv(
  ftpc,
  gzfile("intermediate/feature_history_correlation.bacterial_subtrees.tab.gz")
)

# ftpc = read_tsv("intermediate/feature_history_correlation.bacterial_subtrees.tab.gz")

# Select only the p-value data
ftpv = ftpc %>%
  select(-Calvin, -Median, -MAD, -Mean, -SD) %>%
  distinct()

# Use correlation R to select 3 most important Features in each Subtree
topf = ftpv %>%
  filter(is.finite(pWilcox)) %>%
  mutate(
    padjW = p.adjust(pWilcox, method="BH"),
    padjC = p.adjust(pCor, method="BH")
  ) %>%
  filter(padjW < 0.05 & padjC < 0.05) %>%
  group_by(Subtree) %>%
  top_n(5, abs(R))

# Plot trees for top Features
garbage = foreach(i=1:nrow(topf)) %dopar% {

  iftp = topf[i,]
  f = iftp$Feature
  k = iftp$Subtree
  R = iftp$R
  pW = iftp$padjW
  pC = iftp$padjC

  # Extract subtree
  sbtr = drop.tip(
    batr,
    setdiff(batr$tip.label, filter(subt, Subtree == k)$Accession)
  )

  # Perform Calvin ACE on subtree
  cpos = ifelse(sbtr$tip.label %in% posg, "Positive", "Negative")
  names(cpos) = sbtr$tip.label
  cACE = ace(cpos, sbtr, model="ER", type="discrete")

  # Create feature vector
  baft = feat %>%
    filter(Feature == f) %>%
    filter(Accession %in% sbtr$tip.label) %>%
    arrange(match(Accession, sbtr$tip.label)) %>%
    pull(Count)

  names(baft) = sbtr$tip.label

  # Perform ACE on Feature using Brownian Motion model (default)
  baftACE = as.double(fastAnc(sbtr, baft))

  # Prepare extra variables
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

  # Create dataframe with group association for heatmap
  tfls = txfl[sbtr$tip.label, , drop=F]

  # Plot tree
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
    gp, tfls,
    offset = 0.05, width = 0.05, colnames = F,
    color = NA
  )
  gp = gp + scale_fill_identity()
  gp = gp + geom_treescale(x=0, y=0, fontsize=3)

  sclf = max(c(sqrt(length(sbtr$tip.label) / length(artr$tip.label)), 1))

  ggsave(
    paste(
      "results/ace_bacteria_subtrees_",
      k, "_", f,
      "_R", round(R, 3),
      "_Wilcoxlog10padj", round(log10(pW), 1),
      "_Corlog10padj", round(log10(pC), 1),
      ".pdf",
      sep=""
    ),
    gp, w=sclf*15, h=sclf*12, units="cm", limitsize=F
  )
}
