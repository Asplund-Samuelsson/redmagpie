options(width=150)
library(tidyverse)

# Define infiles
supE_file = "results/Supplementary_Enrichment.tab"
supA_file = "results/Supplementary_ACE.tab"
supR_file = "results/Supplementary_Random_forest.tab"

# Load data
supE = read_tsv(supE_file)
supA = read_tsv(supA_file)
supR = read_tsv(supR_file)

# Count DUFs
filter(supE, Feature_Type == "Pfam") %>%
	select(Feature, Name, Rank) %>%
	distinct() %>%
	top_n(100, -Rank) %>%
	filter(startsWith(Name, "DUF") | startsWith(Name, "UPF")) %>%
	pull(Feature) %>%
	unique %>%
	length

filter(supA, Feature_Type == "Pfam") %>%
	select(Feature, Name, Rank) %>%
	distinct() %>%
	top_n(100, -Rank) %>%
	filter(startsWith(Name, "DUF") | startsWith(Name, "UPF")) %>%
	pull(Feature) %>%
	unique %>%
	length

filter(supR, Feature_Type == "Pfam") %>%
	select(Feature, Name, Rank) %>%
	filter(startsWith(Name, "DUF") | startsWith(Name, "UPF")) %>%
  arrange(Rank)
