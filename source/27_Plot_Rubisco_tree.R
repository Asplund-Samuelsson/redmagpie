options(width=150)
library(tidyverse)
library(ggtree)
library(ggnewscale)

# Define infiles
rbtr_file = "intermediate/rubisco_with_Tabita_examples.tree"
rlps_file = "data/2007_Tabita_Rubisco_examples.tab"
rubi_file = "data/rubisco.txt"
tran_file = "data/2007_Tabita_Rubisco_accessions.tab"
allp_file = "data/positive_genomes.txt"
orgs_file = "intermediate/accession_organism_colours.tab"
meta_file = "data/gtdb_metadata.tab.gz"

# Load data
rbtr = read.tree(rbtr_file)
rlps = read_tsv(rlps_file)
rubi = scan(rubi_file, character())
tran = read_tsv(tran_file, col_names=c("Accession", "label"))
allp = scan(allp_file, character())
orgs = read_tsv(orgs_file)
meta = read_tsv(meta_file)

# Assemble annotation data on Rubisco forms
rann = bind_rows(
  inner_join(filter(tran, !is.na(label)), rlps) %>% select(label, Form),
  tibble(label = rubi) %>%
    mutate(label = unlist(lapply(str_split(label, ":"), "[[", 1))) %>%
    mutate(Genome = unlist(lapply(str_split(label, "\\|"), "[[", 1))) %>%
    mutate(Form = ifelse(Genome %in% allp, "CBB-positive", "Other")) %>%
    select(-Genome)
) %>%
  mutate(Form = ifelse(grepl("III", Form), "III", Form)) %>%
  as.data.frame()
row.names(rann) = rann$label

# Find most recent common ancestor of all RLPs
library(phytools)

mrca = findMRCA(
  rbtr,
  filter(
    tran,
    Accession %in% filter(rlps, grepl("^IV", Form))$Accession,
    !is.na(label)
  )$label
)

# Reroot at RLP
rbtr = reroot(rbtr, mrca)

mrca = findMRCA(
  rbtr,
  filter(
    tran,
    Accession %in% filter(rlps, grepl("^IV", Form))$Accession,
    !is.na(label)
  )$label
)

# Set up color scale
rubc = c(
  "#e0e0e0",
  "#00441b", "#1b7837", "#5aae61", "#a6dba0", "#d9f0d3",
  "#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8",
  "#000000"
)

# Get phylogenetic group
meta = meta %>% select(accession, gtdb_taxonomy)
meta = meta %>%
  mutate(
    Domain = unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 1)),
    Phylum = unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 2)),
    Class = unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 3)),
    Order = unlist(lapply(str_split(gtdb_taxonomy, ";"), "[[", 4))
  ) %>%
  select(-gtdb_taxonomy) %>%
  mutate(
    Domain = str_replace(Domain, "d__", ""),
    Phylum = str_replace(Phylum, "p__", ""),
    Class = str_replace(Class, "c__", ""),
    Order = str_replace(Order, "o__", "")
  )

meta = mutate(meta, Group = ifelse(Phylum == "Proteobacteria", Class, Phylum))
meta = mutate(meta, Group = unlist(lapply(str_split(Group, "_"), "[[", 1)))

meta = rename(meta, Accession = accession)

# Add colours to all accessions
meta = meta %>%
  left_join(distinct(select(orgs, -Accession, -CBB, -Organism, -Colour))) %>%
  mutate(
    Organism0 = ifelse(is.na(Organism0), paste("Other", Group), Organism0)
  ) %>%
  left_join(distinct(select(orgs, Organism0, Organism))) %>%
  mutate(Organism = replace_na(Organism, "Other")) %>%
  inner_join(distinct(select(orgs, Organism, Colour)))

# Add colour for Cyanobacteria
meta = meta %>%
  mutate(
    Organism = ifelse(Group == "Cyanobacteria", "Cyanobacteria", Organism),
    Colour = ifelse(Group == "Cyanobacteria", "#35978f", Colour)
  )

# Plot tree
gp = ggtree(rbtr, layout="fan", size=0.2)

gp = gheatmap(
  gp, select(rann, -label) %>% filter(Form != "Other"),
  offset = 0.05, width = 0.05, colnames = F, color=NA
)
gp = gp + scale_fill_manual(values=rubc, name="Rubisco form")
gp = gp + geom_treescale(x=1.5, y=0, fontsize=3, offset=20, linesize=0.2)
gp = gp + geom_hilight(node=mrca, fill="#fae7fb")
gp = gp + new_scale_fill()

# Filter taxonomic colours to organisms in Rubisco tree
txfl = gp$data %>%
filter(isTip) %>%
select(label) %>%
mutate(Accession = unlist(lapply(str_split(label, "\\|"), "[[", 1))) %>%
left_join(select(meta, Accession, Colour))

# Create dataframe with group association for heatmap
txfl = data.frame(
  row.names = txfl$label,
  Organism = txfl$Colour,
  stringsAsFactors=F
)

gp = gheatmap(
  gp, txfl,
  offset = 0.5, width = 0.05, colnames = F, color=NA
)
gp = gp + scale_fill_identity()

ggsave(
  "results/rubisco_with_Tabita_examples.pdf", gp,
  width=24, height=24, units="cm"
)
