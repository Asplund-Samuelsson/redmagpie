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

# Calculate distances between tips
bdis = cophenetic.phylo(batr)

bdtb = bdis %>%
  as_tibble() %>%
  mutate(Accession = colnames(.)) %>%
  gather(Relative, Distance, -Accession)

# Use hclust
bhcl = hclust(as.dist(bdis))

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

kdis = bdis

# Find optimal number of clusters as close as possible to Arch. tree properties
find_nopt = function(kdis) {

  aksm = foreach(n=1:nrow(kdis)) %dopar% {
    # Determine the subtrees
    clst = cutree(hclust(as.dist(kdis)), k=n)
    clst = tibble(Accession = names(clst), Cluster = clst)
    # Determine non-trivial clusters (size > 1)
    ntrv = clst %>%
      group_by(Cluster) %>%
      summarise(Count = length(Cluster)) %>%
      filter(Count > 1) %>%
      pull(Cluster)
    # For each cluster...
    clsi = lapply(ntrv, function(k){
      # Extract the subtree
      cltr = drop.tip(
        batr, setdiff(batr$tip.label, filter(clst, Cluster == k)$Accession)
      )
      suppressMessages(
        tibble(
          Calvin = ifelse(cltr$tip.label %in% posg, 1, 0)
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
              exp(-abs(log(max(nodeHeights(cltr))/ahgt))),
              exp(-abs(log((sd(cltr$edge.length)/mean(cltr$edge.length)/aecv))))
            )
          )
        ) %>%
        mutate(Cluster = k, Size = length(cltr$tip.label), Clusters = n)
      )
    }) %>% bind_rows()
  } %>% bind_rows()

  # Find the optimal cluster number
  aksa = aksm %>%
    group_by(Clusters, Cluster, Size) %>%
    summarise(Similarity = prod(Factor)) %>%
    arrange(-Similarity)

  akop = aksa %>%
    group_by(Clusters) %>%
    summarise(Similarity = mean(Similarity), Size = mean(Size)) %>%
    arrange(-Similarity)

  nopt = akop$Clusters[1]

  if (nopt == 1){
    return(rownames(kdis))
  } else {
    # Recursion
    clst = cutree(hclust(as.dist(kdis)), k=nopt)
    clst = tibble(Accession = names(clst), Cluster = clst)
    lapply(1:nopt, function(n){
      accs = filter(clst, Cluster == n)$Accession
      find_nopt(kdis[accs, accs])
    })
  }
}

# Get a hierarchical list of lists for the optimal clusters
opcl = find_nopt(bdis)

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
oclf = flatten2(opcl)

# Create cluster table
oclt = bind_rows(lapply(
  1:length(oclf), function(i){tibble(Cluster = i, Accession = oclf[[i]])}
))

# Check what kind of clusters are in there
ocls = oclt %>%
  group_by(Cluster) %>%
  summarise(
    Count = length(Cluster),
    Calvin = sum(Accession %in% posg)/length(Cluster)
  )

# Check the distribution of the clusters on the phylogenetic tree
gp = ggtree(batr, layout="fan", )
gp$data = left_join(
  gp$data,
  oclt %>% mutate(Cluster = as.character(Cluster)) %>% rename(label = Accession)
)
gp = gp + geom_tippoint(mapping=aes(fill=Cluster), shape=24)
gp = gp + scale_fill_manual(
    values = c(
      "#bf812d", "#f6e8c3", "#9970ab", "#e7d4e8",
      "#35978f", "#c7eae5", "#1b7837", "#a6dba0",
      "#8c510a", "#dfc27d", "#762a83", "#c2a5cf",
      "#01665e", "#80cdc1", "#5aae61", "#d9f0d3"
    )
  )

ggsave("results/ace_bacterial_clusters.pdf", gp, w=40, h=40, units="cm")

# Perform correlation and testing on clusters larger than 20 genomes
ftpc = bind_rows(lapply(filter(ocls, Count > 20)$Cluster, function(k){
  # Extract subtree
  cltr = drop.tip(
    batr,
    setdiff(batr$tip.label, filter(oclt, Cluster == k)$Accession)
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
        Cluster = k,
        Feature = f
      ) %>%
      select(
        Feature, Cluster, Calvin,
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
  gzfile("intermediate/feature_history_correlation.bacterial_clusters.tab.gz")
)

# Select only the p-value data
ftpv = ftpc %>%
  select(-Calvin, -Median, -MAD, -Mean, -SD) %>%
  distinct()

# Use correlation R to select 3 most important Features in each cluster
topf = ftpv %>%
  filter(is.finite(pWilcox)) %>%
  mutate(padj = p.adjust(pWilcox, method="BH")) %>%
  filter(padj < 0.05) %>%
  group_by(Cluster) %>%
  top_n(3, abs(R))

# Plot trees for top Features
garbage = foreach(i=1:nrow(topf)) %dopar% {

  iftp = topf[i,]
  f = iftp$Feature
  k = iftp$Cluster
  R = iftp$R
  p = iftp$padj

  # Extract subtree
  cltr = drop.tip(
    batr,
    setdiff(batr$tip.label, filter(oclt, Cluster == k)$Accession)
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
      "results/ace_bacteria_clusters_",
      k, "_", f, "_R", round(R, 3),
      "_Wilcoxlog10padj", round(log10(p), 1), ".pdf",
      sep=""
    ),
    gp, w=sclf*15, h=sclf*12, units="cm", limitsize=F
  )
}
