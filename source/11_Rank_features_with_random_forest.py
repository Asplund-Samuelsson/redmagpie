#!/usr/bin/env python3

# Import libraries
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
import numpy as np
import pandas as pd
from decimal import Decimal
from random import sample
import argparse
from collections import Counter
import gzip

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--feature_names', type=str, help='Load feature names.')
parser.add_argument('--accession_ids', type=str, help='Load accession IDs.')
parser.add_argument('--features', type=str, help='Load features.')
parser.add_argument('--classes', type=str, help='Load classes.')
parser.add_argument('--taxonomy', type=str, help='Load taxonomy.')
parser.add_argument('--importances', type=str, help='Save importances.')
parser.add_argument('--predictions', type=str, help='Save predictions.')
parser.add_argument(
    '--var_cut', action='store_true', default=False,
    help='Select features with variance before phylogenetic filter.'
)
parser.add_argument(
    '--lgr_cut', action='store_true', default=False,
    help='Select features with logistic regression after phylogenetic filter.'
)
args = parser.parse_args()

# Define infiles
feature_name_file = args.feature_names
accession_id_file = args.accession_ids
X_file            = args.features
y_file            = args.classes
taxonomy_file     = args.taxonomy

# Define outfiles
importance_file   = args.importances
prediction_file   = args.predictions

# Store flags
var_cut = args.var_cut
lgr_cut = args.lgr_cut

# Define functions
def round_to(x, n=1):
    """Round number to n significant digits"""
    # https://stackoverflow.com/a/3413529/6018441
    if x == 0:
        return 0
    else:
        return np.round(x, n - int(np.floor(np.log10(abs(x)))) - 1)

def find_bias(accns, features):
    # Determine what taxonomic groups have more than 50 examples
    tax = Counter(list(taxonomy[taxonomy['Accession'].isin(accns)]['Group']))
    groups = [x for x in tax if tax[x] > 50]
    # Filter out features that identify large taxonomic groups (>50 members)
    all_biased_features = set()
    for group in groups:
        group_biased_features = set()
        group_acc = list(taxonomy[taxonomy['Group'] == group]['Accession'])
        group_labels = np.array([int(x in group_acc) for x in accns])
        for n in range(0, 100):
            # Train random forest
            rf = RandomForestClassifier(
                n_jobs = 16, n_estimators = 1000, oob_score = True,
                class_weight = "balanced"
            )
            junk = rf.fit(features, group_labels)
            # Determine what features make up 50% of the importance
            f = [
                x for _, x in sorted(zip(rf.feature_importances_, \
                feature_list), reverse = True)
            ]
            i = sorted(rf.feature_importances_, reverse = True)
            forest_biased_features = set(np.array(f)[np.cumsum(i) <= 0.5])
            # Add biased features to set for the group
            group_biased_features.update(forest_biased_features)
            print(group + " " + str(n))
            print(" Biased features: " + str(len(forest_biased_features)))
            print("             OOB: " + str(round_to(rf.oob_score_, 3)))
            print("")
        # Add identified biased features to the final set
        all_biased_features.update(group_biased_features)
        n_group = str(len(group_biased_features))
        print("Biased features for " + group + ": " + n_group)
        print("")
    return all_biased_features

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

if var_cut:
    # Filter features based on variance
    varf = np.var(X, axis=0).argsort()[-3000:][::-1]
    X = X[:,varf]
    feature_list = list(np.array(feature_list)[varf])

# Find phylogenetically biased features for positive genomes
bias_1 = find_bias(np.array(accession_list)[y == 1], X[y == 1,:])
print("Found " + str(len(bias_1)) + " biased features for positive genomes.")
print("")

# Find phylogenetically biased features for negative genomes
bias_0 = find_bias(np.array(accession_list)[y == 0], X[y == 0,:])
print("Found " + str(len(bias_0)) + " biased features for negative genomes.")
print("")

# Report the total number of phylogenetically biased features identified
phylogenetically_biased_features = bias_1 | bias_0
print("Total biased features: " + str(len(phylogenetically_biased_features)))
print("")

# Filter the input features to those that are not phylogenetically biased
f = np.array([x not in phylogenetically_biased_features for x in feature_list])
X_f = X[:,f]
feature_list_f = list(np.array(feature_list)[f])

if lgr_cut:
    # Select most promising features using logistic regression
    rfe = RFE(
        estimator=LogisticRegression(),
        n_features_to_select=int(np.ceil(0.25*X_f.shape[1])),
        step=10, verbose=0
    )
    rfe.fit(X_f, y)
    rfe_support = rfe.get_support()
    X_f = X_f[:,rfe_support]
    feature_list_f = list(np.array(feature_list_f)[rfe_support])

# Open feature importance file and write header
imp_file = open(importance_file, 'w')
imp_file.write("Forest\tFeature\tImportance\n")
prd_file = open(prediction_file, 'w')
prd_file.write("Forest\tAccession\tClass\tPrediction\n")

# Train and test 100 forests, storing feature importances
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
    rf = RandomForestClassifier(n_jobs=16, n_estimators=500, oob_score=True)
    junk = rf.fit(train_features, train_labels)
    # Predict test labels
    predictions = rf.predict(test_features)
    # Extract feature importances
    i = rf.feature_importances_
    # Write to outfiles
    imp = zip([n] * len(i), feature_list_f, i)
    imp_out = "\n".join(["\t".join([str(a) for a in x]) for x in imp]) + "\n"
    junk = imp_file.write(imp_out)
    prd = zip([n] * len(test_labels), test_accessions, test_labels, predictions)
    prd_out = "\n".join(["\t".join([str(a) for a in x]) for x in prd]) + "\n"
    junk = prd_file.write(prd_out)
    # Report accuracy
    accu = sum([int(x) for x in predictions == test_labels]) / len(test_labels)
    if n == 0:
        asum = accu
    else:
        asum += accu
    print(" Accuracy: " + str(round_to(accu, 3)))
    print("      OOB: " + str(round_to(rf.oob_score_, 3)))
    print(" Av. acc.: " + str(round_to(asum/(n+1), 3)))
    print("")

# Close files
imp_file.close()
prd_file.close()
