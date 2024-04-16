process SVC_PERMUTE{

    publishDir "$params.outdir/svctables", pattern: '*.png'
    label 'big_job'

    input: 
    tuple path(samplesheet), path(countsframe) 
    val x

    output:
    path "svc_table.csv", emit: outfile
    path "nDEGs.txt", emit: nDEGs
    path "svc_permute_fs_*.png"
    path "rfecv.pickle"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import os
    import csv
    import matplotlib.pyplot as plt
    import pickle

    from sklearn import preprocessing
    from sklearn.feature_selection import RFECV
    from sklearn.svm import SVC
    from sklearn.model_selection import StratifiedKFold

    # Input data
    # Reading the Dataset
    count_data=pd.read_table('$countsframe',sep=',',index_col=0) 
    samplesheet=pd.read_csv('$samplesheet')

    # Keep dataset small while testing
    pare = True if ('$params.pareml' == "true") else False
    print("Paring: " + str(pare))
    if(pare):
        count_data=count_data.head(n=100)
        
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
    
    # Use a larger step size if debugging, to speed the pipeline
    step = 50 if pare else $params.stepsize
    
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
    # Save for later
    with open('rfecv.pickle', 'wb') as handle:
        pickle.dump(rfecv, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

    # Get outputs
    print(f"Optimal number of features: {rfecv.n_features_}")
    # Get the names of the selected features
    DEG_indices = rfecv.get_support(indices=True)
    DEGnames = count_data.axes[1][DEG_indices].tolist()

    # Write number of DEGs
    with open('nDEGs.txt','w') as out:
        out.write(f"{rfecv.n_features_}\\n")
        out.close()

    # Write output frame with all genes and their DEG status (1=DEG)
    allgenes=pd.DataFrame(count_data.axes[1].to_list())
    genescores = allgenes.isin(DEGnames).astype('int')
    outframe = pd.concat([allgenes,genescores], axis=1)
    outframe.columns = ['gene', 'DEGstatus']
    outframe.to_csv("./svc_table.csv",sep =',',index=False)
    
    n_scores = len(rfecv.cv_results_["mean_test_score"])
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Mean test accuracy")
    plt.errorbar(
        range(1, n_scores + 1),
        rfecv.cv_results_["mean_test_score"],
        yerr=rfecv.cv_results_["std_test_score"],
    )
    #plt.title(str(nsamples) + " samples and " + str(nfeatures) + " features of which " + str(ninformative) + " informative, with " + str(nclasses) + " classes")
    plt.savefig("svc_permute_fs_" + str($x) + ".png", format="png")
    plt.close("all")
    
    print("test")
    """
}
