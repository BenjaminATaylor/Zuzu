process SVC_FULLSYNTH{

    input:
    tuple val(depth), 
          path(truedegs), //"trueDEGs.RData"
          path(samplesheet), //"synthsheet.csv"
          path(countsframe), //"synthcounts.csv"
          val(refnum)

    output:
    tuple val(depth), 
          path(truedegs),
          val(refnum),
          path("svc_table.csv")

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import os
    import csv

    from sklearn import preprocessing
    from sklearn.feature_selection import RFECV
    from sklearn.svm import SVC
    from sklearn.model_selection import StratifiedKFold

    # Input data
    # Reading the Dataset
    count_data=pd.read_table('$countsframe',sep=',',index_col=0) 
    samplesheet=pd.read_csv('$samplesheet')

    DEBUG=True

    # Keep dataset small while testing
    if(DEBUG):
        count_data=count_data.head(n=2000)

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

    # Write DEG output to table
    allgenes=pd.DataFrame(count_data.axes[1].to_list())
    genescores = allgenes.isin(DEGnames).astype('int')
    outframe = pd.concat([allgenes,genescores], axis=1)
    outframe.columns = ['gene', 'DEGstatus']
    outframe.to_csv("./svc_table.csv",sep =',',index=False)
    """
}