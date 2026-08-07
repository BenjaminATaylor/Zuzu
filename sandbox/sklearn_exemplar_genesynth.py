from sklearn import svm
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import OneHotEncoder, LabelEncoder

scale = True 
classifier = "SVC"
usetuner = True

synthcounts = pd.read_csv('input/synthtest/synth_counts.csv', index_col = 0)
synthmetadata = pd.read_csv('input/synthtest/synth_metadata.csv', index_col = 0)
synthannot = pd.read_csv("input/synthtest/synth_annotations.csv", index_col = 0)

# format for sklearn
X = synthcounts.transpose()
le = LabelEncoder()
le.fit(synthmetadata.phenotype.array)
y = le.transform(synthmetadata.phenotype.array)
    
if scale:
    # Standardize the features
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
else:
    X = np.array(X)

## Here we plot to check that the correct genes are showing differences in means, as expected, and visualise the effect sizes
# Calculating the means for each group of phenotypes
group1 , group2 = np.mean(X[:30, :], axis=0) , np.mean(X[30:, :], axis=0)
# Calculating the differences between the means
absdiffs = abs(group1-group2)
# Plotting
fig, ax = plt.subplots()
ax.bar(np.arange(len(absdiffs)) + 1, absdiffs)
ax.set_ylabel('Difference in Means')
plt.show()

# What's clear here is that the effects sizes produced by compCodeR, after scaling, are quite strong

# How do they compare to effect sizes observed in a real dataset?
realcounts = pd.read_table('./input/test3/GSE91061_BMS038109Sample.hg19KnownGene.raw.tsv',sep='\t',index_col=0) 
realmetadata = pd.read_csv('./input/test3/samplesheet_phenotype.csv')

# format for sklearn
X = realcounts.transpose()
le = LabelEncoder()
le.fit(realmetadata.phenotype.array)
y = le.transform(realmetadata.phenotype.array)

if scale:
    # Standardize the features
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
else:
    X = np.array(X)

size1 = y.tolist().count(1)

# Calculating the means for each group of phenotypes
group1 , group2 = np.mean(X[:size1, :], axis=0) , np.mean(X[size1:, :], axis=0)
# Calculating the differences between the means
absdiffs = abs(group1-group2)
# Plotting sorted values
fig, ax = plt.subplots()
ax.bar(np.arange(len(absdiffs)) + 1, np.sort(absdiffs)[::-1])
ax.set_ylabel('Difference in Means')
plt.show()
# Okay, interesting. For a real dataset scaled means rarely differ by more than 0.4-0.5, so maybe this is where we should be setting the value for the gensynth data







min_features_to_select = 1  # Minimum number of features to consider

if classifier=="SVC":
    grid = {'C': [0.01, 0.1, 1, 10, 100]}
    estimator = SVC(kernel='linear')
    thisclf = GridSearchCV(estimator, param_grid=grid, cv=5)
    thisclf.fit(X, y)
    clf = thisclf.best_estimator_
else:
    clf = LogisticRegression()
    classifier=""

cv = StratifiedKFold(5)
    
rfecv = RFECV(
    estimator=clf,
    step=1,
    cv=cv,
    scoring="accuracy",
    min_features_to_select=min_features_to_select,
    n_jobs=2,
)
rfecv.fit(X, y)

print(f"Optimal number of features: {rfecv.n_features_}")

n_scores = len(rfecv.cv_results_["mean_test_score"])
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Mean test accuracy")
plt.errorbar(
    range(min_features_to_select, n_scores + min_features_to_select),
    rfecv.cv_results_["mean_test_score"],
    yerr=rfecv.cv_results_["std_test_score"],
)
#plt.title(str(nsamples) + " samples and " + str(")
plt.savefig("tempfigs/genesynth_60samples_1000_200informative_depth1e6_effect0.5.png")
