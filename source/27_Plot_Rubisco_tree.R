options(width=150)
library(tidyverse)
library(ggtree)

# Define infiles
rbtr_file = "intermediate/rubisco_with_Tabita_examples.tree"
rlps_file = "data/2007_Tabita_Rubisco_examples.tab"
rubi_file = "data/rubisco.txt"
tran_file = "data/2007_Tabita_Rubisco_accessions.tab"
allp_file = "data/positive_genomes.txt"

# Load data
rbtr = read.tree(rbtr_file)
rlps = read_tsv(rlps_file)
rubi = scan(rubi_file, character())
tran = read_tsv(tran_file, col_names=c("Accession", "label"))
allp = scan(allp_file, character())

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

# Set up color scale
rubc = c(
  "#e0e0e0",
  "#00441b", "#1b7837", "#5aae61", "#a6dba0", "#d9f0d3",
  "#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8",
  "#000000"
)

# Plot tree
gp = ggtree(rbtr, layout="fan", size=0.2)
gp = gheatmap(
  gp, select(rann, -label) %>% filter(Form != "Other"),
  offset = 0.05, width = 0.1, colnames = F, color=NA
)
gp = gp + scale_fill_manual(values=rubc, name="Rubisco form")
gp = gp + geom_treescale(x=1.5, y=0, fontsize=3, offset=20, linesize=0.2)
gp = gp + geom_hilight(node=mrca, fill="#fae7fb")

ggsave(
  "results/rubisco_with_Tabita_examples.pdf", gp,
  width=24, height=24, units="cm"
)
