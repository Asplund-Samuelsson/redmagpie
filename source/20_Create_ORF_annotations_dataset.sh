(
  echo -e "Accession\tORF\tFeature";
  zcat intermediate/deepec.tab.gz;
  zcat intermediate/pfam.tab.gz
) | pigz > results/Supplementary_ORF_annotations.tab.gz
