options(width=150)
library(tidyverse)
library(phytools)
library(MidpointRooter)

# Define infiles
artr_file = "intermediate/archaea.tree"
batr_file = "intermediate/bacteria_midpoint_rooted_tree.Rdata"
posg_file = "data/positive_genomes.txt"
meta_file = "data/gtdb_metadata.tab.gz"

# Load data
artr = read.tree(artr_file)
load(batr_file) # Load "batr" bacterial tree object
posg = scan(posg_file, character())
meta = read_tsv(meta_file) %>%
  select(accession, gtdb_taxonomy, checkm_completeness) %>%
  rename(
    Accession = accession,
    Completeness = checkm_completeness,
    Taxonomy = gtdb_taxonomy
  ) %>%
  mutate(Completeness = Completeness / 100) %>%
  mutate(Domain = ifelse(grepl("d__Archaea", Taxonomy), "Archaea", "Bacteria"))

# Root archaeal tree
artr = midpoint.root2(artr)

# Reduce tree to 100% complete genomes
cmgn = filter(meta, Completeness == 1)$Accession
artr = drop.tip(artr, setdiff(artr$tip.label, cmgn))
batr = drop.tip(batr, setdiff(batr$tip.label, cmgn))

# Add Cyanobacteria to positive genomes
posg = c(posg, filter(meta, grepl("p__Cyanobacteria", Taxonomy))$Accession)

# Function to combine data
comb = function(tr, ac) {
  as_tibble(tr$edge) %>%
  # Extract Parent and Tip nodes from tree, then add Accession from tip labels
  rename(Parent = V1, Tip = V2) %>%
  filter(Tip %in% 1:length(tr$tip.label)) %>%
  mutate(
    Parent = Parent - length(tr$tip.label),
    Accession = tr$tip.label[Tip]
  ) %>%
  # Add the likelihood of Positive or Negative Parent
  inner_join(
    ac$lik.anc %>% as_tibble() %>% mutate(Parent = 1:nrow(ac$lik.anc))
  ) %>%
  # Determine Parent Positivity
  mutate(
    ParentPos  = ifelse(Positive > Negative, 1, 0)
  )
}

# Create vectors with pathway status
apos = ifelse(artr$tip.label %in% posg, "Positive", "Negative")
names(apos) = artr$tip.label
bpos = ifelse(batr$tip.label %in% posg, "Positive", "Negative")
names(bpos) = batr$tip.label

library(foreach)
library(doMC)
registerDoMC(32)

fnfp = foreach(i=1:1000) %dopar% {

  # Sample a completeness value for the positive genomes
  apcm = tibble(
    Accession = names(apos[apos == "Positive"]),
    Completeness = sample(
      filter(meta, Domain == "Archaea")$Completeness, sum(apos == "Positive")
    )
  )

  bpcm = tibble(
    Accession = names(bpos[bpos == "Positive"]),
    Completeness = sample(
      filter(meta, Domain == "Bacteria")$Completeness, sum(bpos == "Positive")
    )
  )

  # Assuming that it is most likely that Rubisco and Prk are associated,
  # and would be missing simultaneously due to missing DNA,
  # sample the CBB state based on Completeness as probability
  apcm = apcm %>%
    group_by(Accession) %>%
    mutate(
      CBB = sample(
        x = c("Positive", "Negative"),
        size = length(Accession),
        prob = c(Completeness, 1 - Completeness)
      )
    )

  bpcm = bpcm %>%
    group_by(Accession) %>%
    mutate(
      CBB = sample(
        x = c("Positive", "Negative"),
        size = length(Accession),
        prob = c(Completeness, 1 - Completeness)
      )
    )

  # Create mock positive/negative vectors
  apom = ifelse(
    artr$tip.label %in% filter(apcm, CBB == "Positive")$Accession,
    "Positive", "Negative"
  )
  names(apom) = artr$tip.label

  bpom = ifelse(
    batr$tip.label %in% filter(bpcm, CBB == "Positive")$Accession,
    "Positive", "Negative"
  )
  names(bpom) = batr$tip.label

  # Perform Ancestral Character Estimation
  aACE = ace(apom, artr, model="ER", type="discrete")
  bACE = ace(bpom, batr, model="ER", type="discrete")

  # Combine data
  apar = comb(artr, aACE)
  bpar = comb(batr, bACE)

  # Check the number of false positives introduced
  arfp = nrow(filter(apar, ParentPos == 1 & !(Accession %in% posg)))
  bafp = nrow(filter(bpar, ParentPos == 1 & !(Accession %in% posg)))

  # ...and the number of false negatives "rescued"
  arfn = nrow(filter(inner_join(apcm, apar), CBB=="Negative" & ParentPos == 1))
  bafn = nrow(filter(inner_join(bpcm, bpar), CBB=="Negative" & ParentPos == 1))

  # Calculate a Benefit which is the number of rescues minus the false positives
  tibble(
    Domain = c("Archaea", "Bacteria"),
    FP = c(arfp, bafp),
    FN_Rescues = c(arfn, bafn)
  ) %>%
    mutate(Benefit = FN_Rescues - FP)
}

fnfp = bind_rows(fnfp)

# Calculate the average benefit
bene = fnfp %>%
  group_by(Domain) %>%
  summarise(CV = sd(Benefit) / mean(Benefit), Benefit = mean(Benefit))
