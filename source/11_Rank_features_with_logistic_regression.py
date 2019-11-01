#!/usr/bin/env python3

# Import libraries
from sklearn.linear_model import LogisticRegression
import numpy as np
import pandas as pd
from decimal import Decimal
from random import sample
import argparse
from collections import Counter
import gzip

# # Parse command line arguments
# parser = argparse.ArgumentParser()
# parser.add_argument('--feature_names', type=str, help='Load feature names.')
# parser.add_argument('--accession_ids', type=str, help='Load accession IDs.')
# parser.add_argument('--features', type=str, help='Load features.')
# parser.add_argument('--classes', type=str, help='Load classes.')
# parser.add_argument('--taxonomy', type=str, help='Load taxonomy.')
# parser.add_argument('--importances', type=str, help='Save importances.')
# parser.add_argument('--predictions', type=str, help='Save predictions.')
# args = parser.parse_args()
#
# # Define infiles
# feature_name_file = args.feature_names
# accession_id_file = args.accession_ids
# X_file            = args.features
# y_file            = args.classes
# taxonomy_file     = args.taxonomy

feature_name_file = "intermediate/EC_count_features.feature_names.txt"
accession_id_file = "intermediate/EC_count_features.accession_ids.txt"
X_file            = "intermediate/EC_count_features.X.csv.gz"
y_file            = "intermediate/EC_count_features.y.txt"
taxonomy_file     = "intermediate/accession_taxonomy.tab"

# # Define outfiles
# importance_file   = args.importances
# prediction_file   = args.predictions

importance_file   = "intermediate/EC_count_features.lr_coefficient.tab"
prediction_file   = "intermediate/EC_count_features.lr_prediction.tab"

# Define functions
def round_to(x, n=1):
    """Round number to n significant digits"""
    # https://stackoverflow.com/a/3413529/6018441
    if x == 0:
        return 0
    else:
        return np.round(x, n - int(np.floor(np.log10(abs(x)))) - 1)

# def find_bias(accns, features):
#     # Determine what taxonomic groups have more than 50 examples
#     tax = Counter(list(taxonomy[taxonomy['Accession'].isin(accns)]['Group']))
#     groups = [x for x in tax if tax[x] > 50]
#     # Filter out features that identify large taxonomic groups (>50 members)
#     all_biased_features = set()
#     for group in groups:
#         group_biased_features = set()
#         group_acc = list(taxonomy[taxonomy['Group'] == group]['Accession'])
#         group_labels = np.array([int(x in group_acc) for x in accns])
#         for n in range(0, 100):
#             # Train random forest
#             rf = RandomForestClassifier(
#                 n_jobs = 16, n_estimators = 1000, oob_score = True,
#                 class_weight = "balanced"
#             )
#             junk = rf.fit(features, group_labels)
#             # Determine what features make up 50% of the importance
#             f = [
#                 x for _, x in sorted(zip(rf.feature_importances_, \
#                 feature_list), reverse = True)
#             ]
#             i = sorted(rf.feature_importances_, reverse = True)
#             forest_biased_features = set(np.array(f)[np.cumsum(i) <= 0.5])
#             # Add biased features to set for the group
#             group_biased_features.update(forest_biased_features)
#             print(group + " " + str(n))
#             print(" Biased features: " + str(len(forest_biased_features)))
#             print("             OOB: " + str(round_to(rf.oob_score_, 3)))
#             print("")
#         # Add identified biased features to the final set
#         all_biased_features.update(group_biased_features)
#         n_group = str(len(group_biased_features))
#         print("Biased features for " + group + ": " + n_group)
#         print("")
#     return all_biased_features

# Read data
X_data = []
for line in gzip.open(X_file, 'rt'):
    line = [float(x) for x in line.strip().split(",")]
    X_data.append(line)


X = np.array(X_data)
y = np.array([int(x.strip()) for x in open(y_file).readlines()])
feature_list = [x.strip() for x in open(feature_name_file).readlines()]
accession_list = [x.strip() for x in open(accession_id_file).readlines()]
taxonomy = pd.read_csv(taxonomy_file, header = 0, sep = "\t")

# # Find phylogenetically biased features for positive genomes
# bias_1 = find_bias(np.array(accession_list)[y == 1], X[y == 1,:])
# print("Found " + str(len(bias_1)) + " biased features for positive genomes.")
# print("")
#
# # Find phylogenetically biased features for negative genomes
# bias_0 = find_bias(np.array(accession_list)[y == 0], X[y == 0,:])
# print("Found " + str(len(bias_0)) + " biased features for negative genomes.")
# print("")
#
# # Report the total number of phylogenetically biased features identified
# phylogenetically_biased_features = bias_1 | bias_0
# print("Total biased features: " + str(len(phylogenetically_biased_features)))
# print("")

B = "intermediate/EC_count_features.biased.txt"
phylogenetically_biased_features = set([x.strip() for x in open(B).readlines()])

# Filter the input features to those that are not phylogenetically biased
f = np.array([x not in phylogenetically_biased_features for x in feature_list])
X_f = X[:,f]
feature_list_f = list(np.array(feature_list)[f])

# Open feature importance file and write header
imp_file = open(importance_file, 'w')
imp_file.write("Classifier\tFeature\tCoefficient\n")
prd_file = open(prediction_file, 'w')
prd_file.write("Classifier\tAccession\tClass\tPrediction\n")

# Train and test 100 classifiers, storing feature coefficients
for n in range(0, 100):
    print("Classification " + str(n))
    # Get randomly sampled indices for test and train datasets
    test_indices = sorted(sample(range(0,len(y)), int(len(y)*0.25)))
    train_indices = sorted(list(set(range(0,len(y))).difference(test_indices)))
    # Extract the accession IDs for the test dataset
    test_accessions = list(np.array(accession_list)[test_indices])
    # Extract test labels and features
    test_labels = y[test_indices]
    test_features = X_f[test_indices,:]
    # Extract training labels and features
    train_labels = y[train_indices]
    train_features = X_f[train_indices,:]
    # Train random forest
    rf = LogisticRegression(C=0.01)
    junk = rf.fit(train_features, train_labels)
    # Predict test labels
    predictions = rf.predict(test_features)
    # Extract feature importances
    i = rf.coef_[0]
    # Write to outfiles
    imp = zip([n] * len(i), feature_list_f, i)
    imp_out = "\n".join(["\t".join([str(a) for a in x]) for x in imp]) + "\n"
    junk = imp_file.write(imp_out)
    prd = zip([n] * len(test_labels), test_accessions, test_labels, predictions)
    prd_out = "\n".join(["\t".join([str(a) for a in x]) for x in prd]) + "\n"
    junk = prd_file.write(prd_out)
    # Report accuracy
    accu = sum([int(x) for x in predictions == test_labels]) / len(test_labels)
    print(" Accuracy: " + str(round_to(accu, 3)))
    # print("      OOB: " + str(round_to(rf.oob_score_, 3)))
    print("")

# Close files
imp_file.close()
prd_file.close()
