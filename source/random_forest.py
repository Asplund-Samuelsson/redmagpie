# Import libraries
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import pandas as pd
from decimal import Decimal
from random import sample

# Define infiles
feature_name_file = "intermediate/rpp_examples_feature_names.no_cyano.txt"
accession_id_file = "intermediate/rpp_examples_accession_ids.no_cyano.txt"
X_file            = "intermediate/rpp_examples_X.no_cyano.csv"
y_file            = "intermediate/rpp_examples_y.no_cyano.txt"

# Read data
X  = np.array(pd.read_csv(X_file, header = None))
y  = np.array([int(x.strip()) for x in open(y_file).readlines()])
feature_list = [x.strip() for x in open(feature_name_file).readlines()]
accession_list = [x.strip() for x in open(accession_id_file).readlines()]

# Define functions
def random_forest(
        train_features, train_labels, rs = None, F = "auto", L = 1, N = 1000
    ):
    """"Train random forest"""
    rf = RandomForestClassifier(
            n_estimators = N, n_jobs = 16, random_state = rs,
            max_features = F, min_samples_leaf = L
        )
    rf.fit(train_features, train_labels)
    return rf

def test_random_forest(rf, test_features, test_labels):
    """Test accuracy of random forest"""
    predictions = rf.predict(test_features)
    return sum(predictions == test_labels) / len(test_labels)

def print_top_features(rf, feature_list, n = 10):
    """Report top features of random forest ordered by feature importance"""
    f = [
            x for _, x in sorted(zip(rf.feature_importances_, \
            feature_list), reverse=True)
        ][:n]
    i = sorted(rf.feature_importances_, reverse = True)[:n]
    for x in zip(f, i):
        print(x[0].lstrip("EC") + "\t" + str(round(x[1], 4)))

def round_to(x, n=1):
    """Round number to n significant digits"""
    # https://stackoverflow.com/a/3413529/6018441
    if x == 0:
        return 0
    else:
        return np.round(x, n - int(np.floor(np.log10(abs(x)))) - 1)

def round_to_same(a, b):
    """Round a to same decimal places as b"""
    # https://stackoverflow.com/q/6189956/6018441
    e = abs(Decimal(str(b)).as_tuple().exponent)
    return np.round(a, e)

def test_repeatedly(X, y, N = 10, F = "auto", L = 1, en = 1000):
    """Test random forest repeatedly and report accuracy with SD"""
    accuracies = []
    for n in range(0, N):
        # Split into training and testing datasets
        train_features, test_features, train_labels, test_labels = \
        train_test_split(X, y, test_size = 0.25)
        # Train random forest
        rf = random_forest(train_features, train_labels, F = F, L = L, N = en)
        # Calculate accuracy
        accuracies.append(test_random_forest(rf, test_features, test_labels))
    # Calculate standard deviation and round to one significant digit
    sd = round_to(np.std(accuracies), 1)
    # Round accuracy to the same number of decimal places as SD
    accuracy = round_to_same(np.mean(accuracies), sd)
    # Report accuracy and standard deviation
    return (accuracy, sd)

# Test hyperparameters

# Number of features
for f in range(0,89):
    print(str(0.02 + f / 100) + "\t", end="")
    print("±".join([str(x) for x in test_repeatedly(X, y, F = 0.02 + f / 100)]))

# Leaf size
for l in range(1, 100, 5):
    print(str(l) + "\t", end="")
    print("±".join([str(x) for x in test_repeatedly(X, y, L = l)]))

# Number of trees
for n in range(3, 13):
    N = int(10**(n/3))
    print(str(N) + "\t", end="")
    print("±".join([str(x) for x in test_repeatedly(X, y, en = N)]))


# Report feature importance repeatedly
outfile = open("results/feature_importance.no_cyano.tab", "w")

for n in range(0, 10):
    rf = random_forest(X, y)
    f = feature_list
    i = rf.feature_importances_
    outfile.write("\n".join([x[0]+"\t"+str(x[1]) for x in zip(f, i)]) + "\n")

outfile.close()

# Calculate accuracy per accession
accession_classification = []

for n in range(0, 1000):
    print(n)
    # Get randomly sampled indices for test and train datasets
    test_indices = sorted(sample(range(0,len(y)), int(len(y)*0.25)))
    train_indices = sorted(list(set(range(0,len(y))).difference(test_indices)))
    # Extract the accession IDs for the test dataset
    test_accessions = list(np.array(accession_list)[test_indices])
    # Extract test labels and features
    test_labels = y[test_indices]
    test_features = X[test_indices,:]
    # Extract training labels and features
    train_labels = y[train_indices]
    train_features = X[train_indices,:]
    # Train random forest
    rf = random_forest(train_features, train_labels)
    # Predict test labels
    predictions = rf.predict(test_features)
    # Encode correct (1) and incorrect classifications (0)
    classifications = [int(x) for x in predictions == test_labels]
    # Save accession and classification correctness
    accession_classification.extend(zip(test_accessions, classifications))

# Summarise per-accession accuracy
accession_accuracy = {}

for test in accession_classification:
    try:
        accession_accuracy[test[0]][0] += test[1]
        accession_accuracy[test[0]][1] += 1
    except KeyError:
        accession_accuracy[test[0]] = [test[1], 1]

# Calculate accuracy per accession and write to file
outfile = open("results/accession_accuracy_n.no_cyano.tab", "w")

for accession in accession_accuracy:
    i = accession_accuracy[accession]
    accuracy = str(i[0] / i[1])
    count = str(i[1])
    outfile.write(accession + "\t" + accuracy + "\t" + count + "\n")

outfile.close()
