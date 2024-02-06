from sklearn import svm
from sklearn.datasets import load_iris
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification

nfeatures = 15
ninformative = 3
sep = 0.8
#nclasses = 2
scale = False # This is actually redundant at the moment since make_classification produces normalized data anyway

classarray = [2,3,4,8]

for nclasses in classarray:

    X, y = make_classification(
        n_samples=500,
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

    ## save simulated data for R testing
    #dfy,dfx = pd.DataFrame(y),pd.DataFrame(X)
    #dfy.to_csv('test_y.csv') 
    #dfx.to_csv('test_x.csv') 
    #X.all
    #
    ## Standardize the features
    #scaler = StandardScaler()
    #X = scaler.fit_transform(X)
    #
    ## also save the scaled version for testing purposes
    #dfy,dfx = pd.DataFrame(y),pd.DataFrame(X)
    #dfx.to_csv('test_x_transform.csv') 


    # Train an SVC model with a linear kernel on the entire dataset
    svc_model = svm.SVC(kernel='linear')
    svc_model.fit(X, y)

    # Initialize the minimum number of features
    min_num_features = 1

    # Lists to store results for plotting
    num_removed_features_list = []
    cv_error_list = []

    while X.shape[1] > min_num_features:
        # Get the coefficients and normalize them
        coef = svc_model.coef_
        coef = np.abs(coef)
        
        # Identify the least informative feature
        least_informative_feature_index = np.argmin(np.sum(coef, axis=0))
        
        # Remove the least-informative feature
        X = np.delete(X, least_informative_feature_index, axis=1)
        
        # Train an SVC model with a linear kernel on the modified dataset
        svc_model = svm.SVC(kernel='linear', coef0=0.0, C=1.0)
        svc_model.fit(X, y)
        
        # Cross-validation score after feature removal
        cv_score = np.mean(cross_val_score(svc_model, X, y, cv=5))
        
        # Store results for plotting
        num_removed_features_list.append(X.shape[1])
        cv_error_list.append(cv_score)

    def thisfunction(i):
        return nfeatures - i

    included_features_list =  list(map(thisfunction, num_removed_features_list))

    figname = "Python_" + str(nfeatures) + "_" + str(ninformative) + "informative_" + str(sep) + "sep_" + str(nclasses) + "class"

    # Plotting the results
    plt.figure(figsize=(10, 6))
    #plt.plot(included_features_list, cv_error_list, marker='o')
    plt.plot(num_removed_features_list, cv_error_list, marker='o')
    plt.title(str(nfeatures) + " features of which " + str(ninformative) + " informative, with " + str(nclasses) + " classes")
    plt.xlabel('Number of Included Features')
    plt.ylabel('Cross-Validation Score')
    plt.grid(True)
    plt.savefig('tempfigs/manual_' + figname + ".png")

