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

# Generate all subtrees
ntip = length(batr$tip.label)
nods = batr$Nnode
sbts = lapply(1:nods + ntip, function(n){extract.clade(batr, n)})

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

# Calculate the similarity of each subtree with > 2 tips to Archaeal tree
sbsf = bind_rows(lapply(
  1:length(sbts),
  function(i){
    if (length(sbts[[i]]$tip.label) > 2) {
      archsim(sbts[[i]]) %>%
        mutate(Subtree = i, Size = length(sbts[[i]]$tip.label))
    }
  }
))

sbsm = sbsf %>%
  # For each Subtree..
  group_by(Subtree) %>%
  # Make sure that all Aspects are included (Positive/Negative can be missing)
  filter(setequal(c("Positive", "Negative", "Height", "CV"), Aspect)) %>%
  # Calculate the Similarity
  summarise(Similarity = prod(Factor), Size = unique(Size)) %>%
  # Add the node in the tree
  mutate(
    Node = Subtree + ntip
  )

# Create a table with Subtrees, listing Accessions
sbtt = bind_rows(lapply(
  sbsm$Subtree,
  function(i){tibble(Subtree = i, Accession = sbts[[i]]$tip.label)}
))

# Calculate what other trees each tree is compatible with
tcmp = bind_rows(lapply(
  sbtt$Subtree %>% unique(),
  function(i){
    accs = filter(sbtt, Subtree == i)$Accession
    sbtt %>%
      filter(Subtree != i) %>%
      group_by(Subtree) %>%
      summarise(Compatible = sum(Accession %in% accs) == 0) %>%
      filter(Compatible) %>%
      select(-Compatible) %>%
      rename(Compatible = Subtree) %>%
      mutate(Subtree = i)
  }
))

# Function to find the Subtrees compatible with all trees in a set of Subtrees
get_compatible = function(ssbt){
  tcmp %>%
    filter(Subtree %in% ssbt) %>%
    group_by(Compatible) %>%
    summarise(Count = length(Compatible)) %>%
    filter(Count == length(ssbt)) %>%
    pull(Compatible)
}

# Iteratively select the best subtrees
sbti = sbsm %>% top_n(1, Similarity) %>% pull(Subtree)
sims = sbsm %>% filter(Subtree %in% get_compatible(sbti))

while (nrow(sims) > 0) {
  sbti = c(sbti, sims %>% top_n(1, Similarity) %>% pull(Subtree))
  sims = sims %>% filter(Subtree %in% get_compatible(sbti))
}

# Compare to previous method
lapply(oclf, function(x){prod(archsim(keep.tip(batr, x))$Factor)}) %>%
  unlist() %>% mean()
# Mean similarity was 0.18
unlist(oclf) %>% length
# Number of genomes was 2524
filter(sbtt, Subtree %in% sbti) %>% nrow
# Number of genomes is now 2508
filter(sbsm, Subtree %in% sbti) %>% pull(Similarity) %>% mean
# ...and mean similarity is 0.075
filter(sbsm, Subtree %in% sbti) %>%
  arrange(-Similarity) %>%
  top_n(16, Similarity) %>%
  pull(Similarity) %>%
  mean()
# Even with the top 16, same number as before, mean similarity is just 0.16

# Create Subtree table
oclt = bind_rows(lapply(
  1:length(oclf), function(i){tibble(Subtree = i, Accession = oclf[[i]])}
))

# Check what kind of Subtrees are in there
ocls = oclt %>%
  group_by(Subtree) %>%
  summarise(
    Count = length(Subtree),
    Calvin = sum(Accession %in% posg)/length(Subtree)
  )

# Check the distribution of the Subtrees on the phylogenetic tree
gp = ggtree(batr, layout="fan", )
gp$data = left_join(
  gp$data,
  oclt %>% mutate(Subtree = as.character(Subtree)) %>% rename(label = Accession)
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

ggsave("results/ace_bacterial_subtrees.pdf", gp, w=40, h=40, units="cm")

# Perform correlation and testing on Subtrees larger than 20 genomes
ftpc = bind_rows(lapply(filter(ocls, Count > 20)$Subtree, function(k){
  # Extract subtree
  cltr = drop.tip(
    batr,
    setdiff(batr$tip.label, filter(oclt, Subtree == k)$Accession)
  )

  # Perform Calvin ACE on subtree
  cpos = ifelse(cltr$tip.label %in% posg, "Positive", "Negative")
  names(cpos) = cltr$tip.label
  cACE = ace(cpos, cltr, model="ER", type="discrete")

  # Analyze
  foreach(f=unique(feat$Feature)) %dopar% {
    # Create feature vector
    baft = feat %>%
      filter(Feature == f) %>%
      filter(Accession %in% cltr$tip.label) %>%
      arrange(match(Accession, cltr$tip.label)) %>%
      pull(Count)

    names(baft) = cltr$tip.label

    # Perform ACE on Feature using Brownian Motion model (default)
    baftACE = fastAnc(cltr, baft)

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
        pStudent = t.test(Enzyme~Calvin, baes)$p.value,
        R = cor(baes$Likelihood, baes$Enzyme),
        pCor = cor.test(baes$Likelihood, baes$Enzyme)$p.value,
        Subtree = k,
        Feature = f
      ) %>%
      select(
        Feature, Subtree, Calvin,
        Median, MAD, pWilcox,
        Mean, SD, pStudent,
        R, pCor
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
  mutate(padj = p.adjust(pWilcox, method="BH")) %>%
  filter(padj < 0.05) %>%
  group_by(Subtree) %>%
  top_n(3, abs(R))

# Plot trees for top Features
garbage = foreach(i=1:nrow(topf)) %dopar% {

  iftp = topf[i,]
  f = iftp$Feature
  k = iftp$Subtree
  R = iftp$R
  p = iftp$padj

  # Extract subtree
  cltr = drop.tip(
    batr,
    setdiff(batr$tip.label, filter(oclt, Subtree == k)$Accession)
  )

  # Perform Calvin ACE on subtree
  cpos = ifelse(cltr$tip.label %in% posg, "Positive", "Negative")
  names(cpos) = cltr$tip.label
  cACE = ace(cpos, cltr, model="ER", type="discrete")

  # Create feature vector
  baft = feat %>%
    filter(Feature == f) %>%
    filter(Accession %in% cltr$tip.label) %>%
    arrange(match(Accession, cltr$tip.label)) %>%
    pull(Count)

  names(baft) = cltr$tip.label

  # Perform ACE on Feature using Brownian Motion model (default)
  baftACE = fastAnc(cltr, baft)

  # Prepare extra variables
  bant = tibble(
    Calvin = cACE$lik.anc %>%
      as_tibble() %>%
      pull(Positive),
    Feature = baftACE
  ) %>%
    mutate(node = 1:nrow(.) + length(cltr$tip.label)) %>%
    bind_rows(
      tibble(
        node = 1:length(cltr$tip.label),
        Feature = baft,
        Calvin = ifelse(cpos == "Positive", 1, 0)
      )
    )

  # Plot tree
  gp = ggtree(cltr, aes(colour=Feature), layout="fan")
  gp$data = gp$data %>% inner_join(bant) # Add extra variables
  gp = gp + geom_nodepoint(mapping=aes(fill=Calvin), shape=21)
  gp = gp + geom_tippoint(mapping=aes(fill=Calvin), shape=24)
  gp = gp + scale_colour_viridis_c(option="B")
  gp = gp + scale_fill_gradient2(
    high="#5ab4ac", mid="#f5f5f5", low="#d8b365", midpoint=0.5
  )

  sclf = max(c(sqrt(length(cltr$tip.label) / length(artr$tip.label)), 1))

  ggsave(
    paste(
      "results/ace_bacteria_subtrees_",
      k, "_", f, "_R", round(R, 3),
      "_Wilcoxlog10padj", round(log10(p), 1), ".pdf",
      sep=""
    ),
    gp, w=sclf*15, h=sclf*12, units="cm", limitsize=F
  )
}
