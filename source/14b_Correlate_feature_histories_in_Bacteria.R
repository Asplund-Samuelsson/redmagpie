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
ecim_file = "intermediate/EC_count_features.importance.tab.gz"
pfim_file = "intermediate/pfam_features.importance.tab.gz"
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
gp +
  geom_tree(size=0.1, mapping=aes(colour=Calvin)) +
  theme_bw() +
  facet_grid(~Domain, scales="free_x", space="free") +
  coord_flip() +
  theme(strip.background=element_blank()) +
  scale_colour_gradient2(
    high="#5ab4ac", mid="#888888", low="#d8b365", midpoint=0.5
  )

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
  # Tree must have at least 2 nodes for correlation and testing
  filter(Size > 3)

# Find optimal subtrees as close as possible to Arch. tree properties
optimal_subtrees = function(n) {
  # Filter to descendants of node n
  sbts = slsm %>%
    filter(Node %in% getDescendants(batr, n))

  # Determine the optimal height
  hopt = sbts %>%
    # Calculate the average Similarity for each Height
    group_by(Height) %>%
    summarise(Similarity = mean(Similarity)) %>%
    distinct() %>%
    top_n(1, Similarity) %>%
    # In case of multiple top Similarity Heights, select the lowest Height
    top_n(1, -Height)

  # Fix for when there is no subtree
  if (nrow(hopt) == 0){hopt=tibble(Similarity = 0)}

  # If Similarity of the tree is better than that of the subtrees...
  if (filter(sims, Node == n)$Similarity > hopt$Similarity){
    # ...return the Accessions of the Subtree
    return(extract.clade(batr, n)$tip.label)
  } else {
    # ...otherwise apply analysis recursively to descendants
    lapply(filter(sbts, Height == hopt$Height)$Node, optimal_subtrees)
  }
}

# Get a hierarchical list of lists for the optimal Subtrees
subo = optimal_subtrees(batr$node.label[1])

# List flattening function by Tommy, 2011
# https://stackoverflow.com/a/8139959/6018441
flatten2 <- function(x) {
  len <- sum(rapply(x, function(x) 1L))
  y <- vector('list', len)
  i <- 0L
  rapply(x, function(x) { i <<- i+1L; y[[i]] <<- x })
  y
}

# Flatten the list
subf = flatten2(subo)

# Create Subtree table
subt = bind_rows(lapply(
  1:length(subf), function(i){tibble(Subtree = i, Accession = subf[[i]])}
))

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

# Save subtree table
write_tsv(
  bind_rows(lapply(
    1:length(subf), function(i){tibble(Subtree = i, Accession = subf[[i]])}
  )),
  "intermediate/ace_bacterial_subtrees.tab"
)

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

library(ggnewscale)

# Create dataframe with group association for heatmap
txfl = data.frame(
  row.names = orgs$Accession,
  Organism = orgs$Colour,
  stringsAsFactors=F
)
tfls = txfl[batr$tip.label, , drop=F]

# Check the distribution of the Subtrees on the phylogenetic tree
gp = ggtree(batr, layout="fan")
gp$data = left_join(
  gp$data,
  subt %>% mutate(Subtree = as.character(Subtree)) %>% rename(label = Accession)
)
gp = gp + geom_tippoint(mapping=aes(fill=Subtree), shape=24)
gp = gp + scale_fill_manual(
    values = c(
      "#bf812d", "#f6e8c3", "#9970ab", "#e7d4e8",
      "#35978f", "#c7eae5", "#1b7837", "#a6dba0",
      "#8c510a", "#dfc27d", "#762a83", "#c2a5cf",
      "#01665e", "#80cdc1", "#5aae61", "#d9f0d3"
    )
  )
gp = gp + new_scale_fill()
gp = gheatmap(
  gp, tfls,
  offset=.05, width=.05, colnames=F,
  color = NA
)
gp = gp + scale_fill_identity()

ggsave("results/ace_bacterial_subtrees.orgs.pdf", gp, w=40, h=40, units="cm")


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
  baftACE = fastAnc(sbtr, baft)

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
  gp

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
