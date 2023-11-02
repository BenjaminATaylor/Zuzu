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
count_data=pd.read_table('./input/test3/GSE91061_BMS038109Sample.hg19KnownGene.raw.tsv',sep='\t',index_col=0)
samplesheet=pd.read_csv('./input/test3/samplesheet_phenotype.csv')

# Check null values
count_data.isnull().values.any()
samplesheet.isnull().values.any()
# No null values in data

# Grab a small subset of the data for testing purposes
# Generate a small training set
count_data=count_data.head(n=100)
samplesheet_train = pd.concat([
    samplesheet.query('phenotype == "Pre-treatment"').head(n=10),
    samplesheet.query('phenotype == "On-treatment"').head(n=10)
])
count_data_train = count_data.loc[:,samplesheet_train.loc[:,"sample"]].transpose()
# Generate a small testing set
samplesheet_test = pd.concat([
    samplesheet.query('phenotype == "On-treatment"').tail(n=10),
    samplesheet.query('phenotype == "Pre-treatment"').tail(n=10)
])
count_data_test = count_data.loc[:,samplesheet_test.loc[:,"sample"]].transpose()

# Switch to numeric encodings
# (TODO: It would be nice if these were ordered with reference level = 0)
for x in [samplesheet_train,samplesheet_test]:
    x['phenotype'] = preprocessing.LabelEncoder().fit_transform(x.phenotype)

# Scale data before SVM
count_data_test_scaled = preprocessing.StandardScaler().fit_transform(count_data_test)
count_data_train_scaled = preprocessing.StandardScaler().fit_transform(count_data_train)

# Check that scaling has worked (data should have 0 mean and unit variance)
count_data_test_scaled.mean(axis=0) ; count_data_test_scaled.std(axis=0)
count_data_train_scaled.mean(axis=0) ; count_data_train_scaled.std(axis=0)
# Means are very close to 0, and STDs are either 0 or 1 (0 when there was no variance anyway)

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
