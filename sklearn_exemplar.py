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

nsamples = 60
nfeatures = 15
ninformative = 5
sep = 2
#nclasses = 2
scale = False # This is actually redundant at the moment since make_classification produces normalized data anyway
classifier = "none"
usetuner = True
tunereport = ""

classarray = [2,3,4,8]

for nclasses in classarray:

    X, y = make_classification(
        n_samples=nsamples,
        n_features=nfeatures,
        n_informative=ninformative,
        n_redundant=0,
        n_repeated=0,
        n_classes=nclasses,
        n_clusters_per_class=1,
        shuffle = False,
        class_sep=sep,
        random_state=0,
    )
    
    if scale:
        # Standardize the features
        scaler = StandardScaler()
        X = scaler.fit_transform(X)

    min_features_to_select = 1  # Minimum number of features to consider
    
    if usetuner:
        tunereport = "withtuner_"
        grid = {'C': [0.1, 1, 10, 100]}
        estimator = SVC(kernel='linear')
        thisclf = GridSearchCV(estimator, param_grid=grid, cv=5)
        thisclf.fit(X, y)
        clf = thisclf.best_estimator_
        classifier=""
    elif classifier=="SVC_":
        clf = SVC(kernel='linear')
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

    figname = "Python_" + str(nsamples) + "samples_" + str(nfeatures) + "_" + str(ninformative) + "informative_" + str(sep) + "sep_" + str(nclasses) + "class"
    if scale:
       figname = figname + "_scaled"
    
    n_scores = len(rfecv.cv_results_["mean_test_score"])
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Mean test accuracy")
    plt.errorbar(
        range(min_features_to_select, n_scores + min_features_to_select),
        rfecv.cv_results_["mean_test_score"],
        yerr=rfecv.cv_results_["std_test_score"],
    )
    plt.title(str(nsamples) + " samples and " + str(nfeatures) + " features of which " + str(ninformative) + " informative, with " + str(nclasses) + " classes")
    plt.savefig('tempfigs/' + classifier + figname + tunereport + ".png")
