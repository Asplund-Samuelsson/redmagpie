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

# Perform analysis for each feature
library(foreach)
library(doMC)
registerDoMC(32)

ftpc = foreach(f=unique(feat$Feature)) %dopar% {
  # Create feature vectors
  arft = feat %>%
    filter(Feature == f) %>%
    filter(Accession %in% artr$tip.label) %>%
    arrange(match(Accession, artr$tip.label)) %>%
    pull(Count)

  names(arft) = artr$tip.label

  baft = feat %>%
    filter(Feature == f) %>%
    filter(Accession %in% batr$tip.label) %>%
    arrange(match(Accession, batr$tip.label)) %>%
    pull(Count)

  names(baft) = batr$tip.label

  # Perform ACE on FBP using Brownian Motion model (default)
  arftACE = fastAnc(artr, arft)
  baftACE = fastAnc(batr, baft)

  # Calculate correlations and p-values
  aaes = tibble(
    Likelihood = aACE$lik.anc[,2],
    Calvin = ifelse(aACE$lik.anc[,2] > 0.5, 1, 0),
    Enzyme = arftACE
  )
  baes = tibble(
    Likelihood = bACE$lik.anc[,2],
    Calvin = ifelse(bACE$lik.anc[,2] > 0.5, 1, 0),
    Enzyme = baftACE
  )

  bind_rows(
    aaes %>%
    group_by(Calvin) %>%
    summarise(
      Median = median(Enzyme), MAD = mad(Enzyme),
      Mean = mean(Enzyme), SD = sd(Enzyme)
    ) %>%
    mutate(
      pWilcox = wilcox.test(Enzyme~Calvin, aaes)$p.value,
      pStudent = t.test(Enzyme~Calvin, aaes)$p.value,
      R = cor(aaes$Likelihood, aaes$Enzyme),
      pCor = cor.test(aaes$Likelihood, aaes$Enzyme)$p.value,
      Domain = "Archaea",
      Feature = f
    ),
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
      Domain = "Bacteria",
      Feature = f
    )
  ) %>% select(
    Feature, Domain, Calvin,
    Median, MAD, pWilcox,
    Mean, SD, pStudent,
    R, pCor
  )
} %>% bind_rows()

# Save table; took 3 hours to create
write_tsv(ftpc, gzfile("intermediate/feature_history_correlation.tab.gz"))

# Select only the p-value data
ftpv = ftpc %>%
  select(-Calvin, -Median, -MAD, -Mean, -SD) %>%
  distinct()

# Plot something
gp = ggplot(ftpv, aes(x=log10(p.adjust(pStudent, m="BH")), y=R, colour=Domain))
gp = gp + geom_point(alpha=0.2, size=0.5)
gp = gp + theme_bw()
gp = gp + scale_colour_manual(values=rev(c("#b35806","#542788")))

ggsave(
  "results/ace_feature_correlation_vs_pStudent.png", gp,
  w=16, h=8, units="cm", dpi=200
)

gp = ggplot(ftpv, aes(x=log10(p.adjust(pWilcox, m="BH")), y=R, colour=Domain))
gp = gp + geom_point(alpha=0.2, size=0.5)
gp = gp + theme_bw()
gp = gp + scale_colour_manual(values=rev(c("#b35806","#542788")))

ggsave(
  "results/ace_feature_correlation_vs_pWilcox.png", gp,
  w=16, h=8, units="cm", dpi=200
)

# Use correlation R to select 10 most important Features in either Domain
topf = ftpv %>% group_by(Domain) %>% top_n(10, abs(R))


# Plot trees for top Features
registerDoMC(20)

garbage = foreach(i=1:nrow(topf)) %dopar% {

  iftp = topf[i,]
  f = iftp$Feature
  D = iftp$Domain
  R = iftp$R

  if (D == "Archaea"){

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

    ggsave(
      paste("results/ace_archaea_", f, "_R", R, ".pdf", sep=""), gp,
      w=15, h=12, units="cm"
    )

  } else {

    # Create feature vectors
    baft = feat %>%
      filter(Feature == f) %>%
      filter(Accession %in% batr$tip.label) %>%
      arrange(match(Accession, batr$tip.label)) %>%
      pull(Count)

    names(baft) = batr$tip.label

    # Perform ACE on FBP using Brownian Motion model (default)
    baftACE = fastAnc(batr, baft)

    # Prepare extra variables
    bant = tibble(
      Calvin = bACE$lik.anc %>%
        as_tibble() %>%
        pull(Positive),
      Feature = baftACE
    ) %>%
      mutate(node = 1:nrow(.) + length(batr$tip.label)) %>%
      bind_rows(
        tibble(
          node = 1:length(batr$tip.label),
          Feature = baft,
          Calvin = ifelse(bpos == "Positive", 1, 0)
        )
      )

    # Plot tree
    gp = ggtree(batr, aes(colour=Feature), layout="fan")
    gp$data = gp$data %>% inner_join(bant) # Add extra variables
    gp = gp + geom_nodepoint(mapping=aes(fill=Calvin), shape=21)
    gp = gp + geom_tippoint(mapping=aes(fill=Calvin), shape=24)
    gp = gp + scale_colour_viridis_c(option="B")
    gp = gp + scale_fill_gradient2(
      high="#5ab4ac", mid="#f5f5f5", low="#d8b365", midpoint=0.5
    )

    ggsave(
      paste("results/ace_bacteria_", f, "_R", R, ".pdf", sep=""), gp,
      w=75, h=60, units="cm", limitsize=F
    )

  }
}

# Calculate average importances
avim = bind_rows(ecim, pfim) %>%
  group_by(Feature) %>%
  summarise(Importance = mean(Importance)) %>%
  arrange(-Importance) %>%
  mutate(Type = ifelse(str_starts(Feature, "EC"), "EC", "Pfam"))

# Combine with ancestral character estimation correlation values
imco = inner_join(avim, select(ftpv, Feature, Domain, R)) %>% filter(!is.na(R))

gp = ggplot(imco, aes(x=Importance, y=abs(R), colour=Domain))
gp = gp + geom_point(alpha=0.2, size=0.5)
gp = gp + theme_bw()
gp = gp + scale_colour_manual(values=rev(c("#b35806","#542788")), guide=F)
gp = gp + facet_grid(Domain~Type)
gp = gp + scale_x_log10()
gp = gp + theme(
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  strip.background = element_blank()
)

ggsave(
  "results/ace_feature_correlation_vs_random_forest_importance.png", gp,
  w=15, h=10, units="cm", dpi=200
)

# Check overlap within top n
N = 100
ftpv %>%
  mutate(Type = ifelse(str_starts(Feature, "EC"), "EC", "Pfam")) %>%
  group_by(Domain, Type) %>%
  mutate(nF = length(Feature) - sum(is.na(R))) %>%
  top_n(N, R) %>%
  mutate(
    RandomForest = Feature %in% (
      avim %>% group_by(Type) %>% top_n(N, Importance)
    )$Feature
  ) %>%
  summarise(
    Overlap = sum(RandomForest) / length(RandomForest)
  )

# Sample randomly
registerDoMC(16)

phnf = ftpv %>%
  mutate(Type = ifelse(str_starts(Feature, "EC"), "EC", "Pfam")) %>%
  group_by(Domain, Type) %>%
  summarise(nF = length(Feature) - sum(is.na(R)))
rfnf = tibble(
  Type = c("EC", "Pfam"),
  nF = c(length(unique(dpec$Feature)), length(unique(pfam$Feature)))
)

allf = tibble(Feature = c(unique(dpec$Feature), unique(pfam$Feature))) %>%
  mutate(
    Type = ifelse(str_starts(Feature, "EC"), "EC", "Pfam"),
    Domain = "Both",
    Method = "RF"
  ) %>%
  bind_rows(
    ftpv %>%
      filter(!is.na(R)) %>%
      select(Feature, Domain) %>%
      mutate(
        Type = ifelse(str_starts(Feature, "EC"), "EC", "Pfam"),
        Method = "ACE"
      )
  )

# Sample N features and calculate overlap
smpl = foreach(i=1:1000) %dopar% {
  lapply(seq(10, 550, 10), function(N){
    smpl = allf %>%
      group_by(Type, Domain, Method) %>%
      sample_n(N, replace=F)
    ovlp = smpl %>%
      filter(Method == "ACE") %>%
      mutate(RandomForest = Feature %in% filter(smpl, Method == "RF")$Feature) %>%
      group_by(Type, Domain) %>%
      summarise(Overlap = sum(RandomForest) / length(RandomForest)) %>%
      mutate(i=i, n=N)
    return(ovlp)
  }) %>% bind_rows()
} %>% bind_rows()

ovlp = lapply(
  seq(10, 550, 10),
  function(N){
    ftpv %>%
      mutate(Type = ifelse(str_starts(Feature, "EC"), "EC", "Pfam")) %>%
      group_by(Domain, Type) %>%
      top_n(N, R) %>%
      mutate(
        RandomForest = Feature %in% (
          avim %>% group_by(Type) %>% top_n(N, Importance)
        )$Feature
      ) %>%
      summarise(
        Overlap = sum(RandomForest) / length(RandomForest)
      ) %>%
      mutate(n = N)
  }) %>% bind_rows() %>%
  inner_join(
    smpl %>%
      group_by(Domain, Type, n) %>%
      summarise(Expected = mean(Overlap), Upper=mean(Overlap) + sd(Overlap))
  ) %>%
  rename(Observed = Overlap) %>%
  gather(Overlap, Fraction, -Domain, -Type, -n, -Upper) %>%
  mutate(
    Upper = ifelse(Overlap == "Expected", Upper, Fraction),
    Lower = Fraction
  )

gp = ggplot(
  ovlp,
  aes(x=n, y=Fraction, group=Overlap, fill=Overlap, ymin=Lower, ymax=Upper)
)
gp = gp + geom_col(position="dodge")
gp = gp + geom_errorbar(position="dodge")
gp = gp + facet_grid(Domain~Type)
gp = gp + theme_bw()
gp = gp + theme(
  axis.text = element_text(colour="black"),
  axis.ticks = element_line(colour="black"),
  strip.background = element_blank()
)
gp = gp + scale_fill_manual(values=rev(c("#b35806","#542788")))

ggsave(
  "results/overlap_in_top_N_ace_correlation_vs_rf_importance.pdf", gp,
  w=40, h=18, units="cm"
)
