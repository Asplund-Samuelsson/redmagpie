options(width=150)
library(tidyverse)
library(phytools)
library(MidpointRooter)

# Define infiles
artr_file = "intermediate/archaea.tree"
batr_file = "intermediate/bacteria.tree"
posg_file = "data/positive_genomes.txt"

# Load data
artr = read.tree(artr_file)
batr = read.tree(batr_file)
posg = scan(posg_file, character())

# Create vectors with pathway status
apos = ifelse(artr$tip.label %in% posg, "Positive", "Negative")
names(apos) = artr$tip.label
bpos = ifelse(batr$tip.label %in% posg, "Positive", "Negative")
names(bpos) = batr$tip.label

# Root trees
artr = midpoint.root2(artr)
batr = midpoint.root2(batr)

# Save rooted bacterial tree since it is a very slow process
save(batr, file="intermediate/bacteria_midpoint_rooted_tree.Rdata")

# load("intermediate/bacteria_midpoint_rooted_tree.Rdata")

# Check that there are no branches with length zero
if (
  sum(artr$edge.length == 0) | sum(batr$edge.length == 0)
) {
  message("Warning: Zero length branches detected.")
}

# Perform Ancestral Character Estimation
aACE = ace(apos, artr, model="ER", type="discrete")
bACE = ace(bpos, batr, model="ER", type="discrete")

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
    # Determine Current and Parent Positivity
    mutate(
      CurrentPos = ifelse(Accession %in% posg, 1, 0),
      ParentPos  = ifelse(Positive > Negative, 1, 0)
    ) %>%
    # Classify evolution History
    mutate(
      History = case_when(
        CurrentPos == ParentPos ~ "Unchanged",
        CurrentPos > ParentPos ~ "Gain",
        CurrentPos < ParentPos ~ "Loss"
      )
    )
}

# Combine data
apar = comb(artr, aACE)
bpar = comb(batr, bACE)

# Save new list of positive genomes by interpreting recent "losses" as positive
write(
  filter(bind_rows(apar, bpar), CurrentPos == 1 | ParentPos == 1)$Accession,
  "intermediate/positive_genomes.ace.txt"
)

# Also save the raw analysis output data
write_tsv(bind_rows(apar, bpar), "intermediate/ace.tab")
