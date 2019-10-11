#!/usr/bin/env bash

# Perform random forest analysis using EC count features
DATA="EC_count_features"
echo "Using ${DATA}"

source/11_Rank_features_with_random_forest.py \
  --features      intermediate/${DATA}.X.csv.gz \
  --classes       intermediate/${DATA}.y.txt \
  --feature_names intermediate/${DATA}.feature_names.txt \
  --accession_ids intermediate/${DATA}.accession_ids.txt \
  --taxonomy      intermediate/accession_taxonomy.tab \
  --importances   intermediate/${DATA}.importance.tab \
  --predictions   intermediate/${DATA}.prediction.tab \
  >               intermediate/${DATA}.random_forest.log \
  2>              intermediate/${DATA}.random_forest.error

# Compress output
pigz intermediate/${DATA}.importance.tab intermediate/${DATA}.prediction.tab

# Perform random forest analysis using alignment features
DATA="alignment_features"
echo "Using ${DATA}"

source/11_Rank_features_with_random_forest.py \
  --features      intermediate/${DATA}.X.csv.gz \
  --classes       intermediate/${DATA}.y.txt \
  --feature_names intermediate/${DATA}.feature_names.txt \
  --accession_ids intermediate/${DATA}.accession_ids.txt \
  --taxonomy      intermediate/accession_taxonomy.tab \
  --importances   intermediate/${DATA}.importance.tab \
  --predictions   intermediate/${DATA}.prediction.tab \
  >               intermediate/${DATA}.random_forest.log \
  2>              intermediate/${DATA}.random_forest.error

# Compress output
pigz intermediate/${DATA}.importance.tab intermediate/${DATA}.prediction.tab
