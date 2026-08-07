from operator import concat
import warnings  
import pandas as pd
import numpy as np
import os
import sklearn
from sklearn import preprocessing
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score

from sklearn.pipeline import Pipeline
from sklearn.pipeline import make_pipeline



from sklearn.naive_bayes import GaussianNB



os.chdir("/Users/benjamin/Repositories/Zuzu")
print("Hello world")

# Reading the Dataset
# A real test dataset
#count_data=pd.read_table('./input/test3/GSE91061_BMS038109Sample.hg19KnownGene.raw.tsv',sep='\t',index_col=0) 
#samplesheet=pd.read_csv('./input/test3/samplesheet_phenotype.csv')
# A synthetic test dataset
#count_data=pd.read_table('./input/synthtest/synth_counts.csv',sep=',',index_col=0)
count_data=pd.read_table('./input/synthtest/synth_counts_semiperm.csv',sep=',',index_col=0)
samplesheet=pd.read_csv('./input/synthtest/synth_metadata.csv', index_col=0)

# Check null values
count_data.isnull().values.any()
samplesheet.isnull().values.any()
# No null values in data

# Grab a small subset of the data for testing purposes
# Generate a small training set
#count_data=count_data.head(n=10000)
# In real life we'll be using filtered data
# For testing purposes apply a meager filter
#count_data = count_data.loc[count_data.mean(axis=1)>5,:]

#samplesheet_train = pd.concat([
#    samplesheet.query('phenotype == "Pre-treatment"').head(n=20),
#    samplesheet.query('phenotype == "On-treatment"').head(n=20)
#])

samplesheet_train = samplesheet

count_data_train = count_data.loc[:,samplesheet_train.loc[:,"sample"]].transpose()
# Make sure there are no 0-variance genes in training set
count_data_train = count_data_train.loc[:,count_data_train.std(axis=0) != 0]

# Apply a scaler 
scaler = preprocessing.StandardScaler().fit(count_data_train)
count_data_train_scaled = scaler.transform(count_data_train)

# Switch to numeric encodings
# (TODO: It would be nice if these were ordered with reference level = 0)
for x in [samplesheet_train]:
    x['phenotype'] = preprocessing.LabelEncoder().fit_transform(x.phenotype)
    # This doesn't need to be loop, but leave like thsi in case we decide to split into train/test




















# Generate a small testing set
samplesheet_test = pd.concat([
    samplesheet.query('phenotype == "On-treatment"').tail(n=10),
    samplesheet.query('phenotype == "Pre-treatment"').tail(n=10)
])
count_data_test = count_data.loc[:,samplesheet_test.loc[:,"sample"]].transpose()


# Get CV error for the training
# Note that by putting both the scaling and SVM into a pipeline here, 
# we perform CV on the whole pipeline instead of just the SVM 
clf = make_pipeline(preprocessing.StandardScaler(), SVC(kernel='rbf'))
scores = cross_val_score(clf, count_data_train, samplesheet_train.phenotype, cv=5)
print("%0.2f accuracy with a standard deviation of %0.2f" % (scores.mean(), scores.std()))
clf.fit(count_data_train, samplesheet_train.phenotype)









# SVC
clf = make_pipeline(StandardScaler(), SVC(kernel="rbf"))
clf.fit(count_data_train, samplesheet_train.phenotype)
clf.score(count_data_train, samplesheet_train.phenotype)
clf.score(count_data_test, samplesheet_test.phenotype)

clf.fit(count_data_train, samplesheet_train.phenotype).score(count_data_test, samplesheet_test.phenotype)
clf.set_params(svc__C=10).fit(count_data_train, samplesheet_train.phenotype).score(count_data_test, samplesheet_test.phenotype)


clf.get_params(count_data_test, samplesheet_test.phenotype)

classifier = SVC(C=1,kernel="rbf")

# CV error for the training data 
cross_val_score(classifier,count_data_train_scaled,samplesheet_train.phenotype, cv=5).mean()


clf = classifier.fit(count_data_train,samplesheet_train.phenotype)
y_pred = clf.predict(count_data_test)
accuracy_score(y_pred, samplesheet_test.phenotype)
