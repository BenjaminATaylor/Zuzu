import pandas as pd
import numpy as np
import os
import csv

from sklearn import preprocessing
from sklearn.feature_selection import RFECV
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold

os.chdir("/Users/benjamin/Repositories/Zuzu")
print("Hello world")

# Reading the Dataset
# A real test dataset
count_data=pd.read_table('./input/test3/GSE91061_BMS038109Sample.hg19KnownGene.raw.tsv',sep='\t',index_col=0) 
samplesheet=pd.read_csv('./input/test3/samplesheet_phenotype.csv')

# Generate a small training set
count_data=count_data.head(n=10000)
# In real life we'll be using filtered data
# For testing purposes apply a meager filter
count_data = count_data.loc[count_data.mean(axis=1)>5,:]

# Transpose (sklearn expects to find samples as rows and features as columns)
count_data = count_data.transpose()

# Apply a scaler 
scaler = preprocessing.StandardScaler().fit(count_data)
count_data_scaled = scaler.transform(count_data)

# Switch to numeric encodings
# (TODO: It would be nice if these were ordered with reference level = 0)
samplesheet['phenotype'] = preprocessing.LabelEncoder().fit_transform(samplesheet.phenotype)

# Rename to fit sklearn conventions
X,y = count_data_scaled,samplesheet.phenotype

# Okay, we're set up to run support vector classification with recusrive feature elimination

# Setup classification and RFE parameters
clf = SVC(kernel='linear')
cv = StratifiedKFold(5)
step = 100
rfecv = RFECV(
    estimator=clf,
    step=step,
    cv=cv,
    scoring="accuracy",
    min_features_to_select=1,
    n_jobs=-1,
    verbose=2
)

# Run FRFECV (this may take a while depending on number of features (higher = longer),
# cv (higher = longer), step size (higher = shorter))
rfecv.fit(X, y)

# Get outputs
print(f"Optimal number of features: {rfecv.n_features_}")
# Get the names of the selected features
DEG_indices = rfecv.get_support(indices=True)
DEGnames = count_data.axes[1][DEG_indices].tolist()

# Write list of DEG names
pd.DataFrame(data={"DEGs": DEGnames}).to_csv("./test.csv",sep =',',index=False)


csv.writer(DEGnames)

pd.DataFrame.to_csv(DEGnames)

# An output plot for testing, if one is required
#import matplotlib.pyplot as plt
#
#n_scores = len(rfecv.cv_results_["mean_test_score"])
#plt.figure()
#plt.xlabel("Number of features selected")
#plt.ylabel("Mean test accuracy")
#plt.errorbar(
#    range(min_features_to_select, n_scores + min_features_to_select),
#    rfecv.cv_results_["mean_test_score"],
#    yerr=rfecv.cv_results_["std_test_score"],
#    ecolor='black',
#    elinewidth=0.1
#)
#plt.title("LR, 8 groups, CV 50")
#plt.show()