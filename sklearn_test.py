from sklearn import svm
from sklearn.datasets import load_iris
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification

totalfeatures = 500

X, y = make_classification(
    n_samples=50,
    n_features=totalfeatures,
    n_informative=50,
    n_redundant=0,
    n_repeated=0,
    n_classes=2,
    n_clusters_per_class=1,
    class_sep=0.5,
    random_state=0,
)

# save simulated data for R testing
dfy,dfx = pd.DataFrame(y),pd.DataFrame(X)
dfy.to_csv('test_y.csv') 
dfx.to_csv('test_x.csv') 

# Standardize the features
scaler = StandardScaler()
X = scaler.fit_transform(X)

# Train an SVC model with a linear kernel on the entire dataset
svc_model = svm.SVC(kernel='linear', coef0=0.0, C=1.0)
svc_model.fit(X, y)

# Initialize the minimum number of features
min_num_features = 5

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
    return totalfeatures - i

included_features_list =  list(map(thisfunction, num_removed_features_list))

# Plotting the results
plt.figure(figsize=(10, 6))
#plt.plot(included_features_list, cv_error_list, marker='o')
plt.plot(num_removed_features_list, cv_error_list, marker='o')
plt.title('Number of Removed Features vs Cross-Validation Score')
plt.xlabel('Number of Included Features')
plt.ylabel('Cross-Validation Score')
plt.grid(True)
plt.show()
plt.close()
