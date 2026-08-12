process SVC_FULLSYNTH{

    publishDir "$params.outdir/svctables", pattern: '*.png'
    label 'big_job'

    input:
    tuple val(depth), 
          path(truedegs), //"trueDEGs.RData"
          path(samplesheet), //"synthsheet.csv"
          path(countsframe), //"synthcounts.csv"
          val(refnum),
          val(altnum),
          val(x)

    output:
    tuple val(depth), 
          path(truedegs),
          val(refnum),
          val(altnum),
          val(x),
          path("svc_table.csv"),
          path("svc_fullsynth_fs_*.png"),
          path("rfecv.pickle")
          

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
    all_gene_names = count_data.index.to_list()

    # Check that the requested group sizes and sample order reached SVC unchanged
    required_sample_columns = {'sample', 'phenotype'}
    missing_sample_columns = required_sample_columns.difference(samplesheet.columns)
    if missing_sample_columns:
        raise ValueError(
            "Synthetic sample sheet is missing columns: " +
            ", ".join(sorted(missing_sample_columns))
        )
    if samplesheet[['sample', 'phenotype']].isna().any().any():
        raise ValueError("Synthetic sample sheet contains missing sample or phenotype values")
    if samplesheet['sample'].duplicated().any():
        raise ValueError("Synthetic sample sheet contains duplicate sample identifiers")
    if count_data.columns.duplicated().any():
        raise ValueError("Synthetic count matrix contains duplicate sample identifiers")

    sample_order = samplesheet['sample'].astype(str).tolist()
    count_order = count_data.columns.astype(str).tolist()
    if count_order != sample_order:
        raise ValueError("Synthetic count columns do not match sample sheet order")

    reflevel = '$params.reflevel'
    expected_refnum = int($refnum)
    expected_altnum = int($altnum)
    if expected_refnum < 1 or expected_altnum < 1:
        raise ValueError("Synthetic group sizes must both be positive")
    observed_refnum = int((samplesheet['phenotype'] == reflevel).sum())
    observed_altnum = len(samplesheet) - observed_refnum
    if samplesheet['phenotype'].nunique() != 2:
        raise ValueError("Synthetic sample sheet must contain exactly two phenotypes")
    if len(samplesheet) != expected_refnum + expected_altnum:
        raise ValueError("Synthetic sample count does not match requested group sizes")
    if count_data.shape[1] != expected_refnum + expected_altnum:
        raise ValueError("Synthetic count matrix does not match requested group sizes")
    if observed_refnum != expected_refnum or observed_altnum != expected_altnum:
        raise ValueError(
            "Synthetic phenotype counts do not match requested group sizes: "
            f"observed ref={observed_refnum}, alt={observed_altnum}; "
            f"expected ref={expected_refnum}, alt={expected_altnum}"
        )

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
        n_jobs=$params.njobs,
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

    # Write DEG output to table
    # Preserve the complete synthetic gene universe even in pare mode. Genes
    # omitted from the debug fit remain valid uncalled features in the output.
    allgenes=pd.DataFrame(all_gene_names)
    genescores = allgenes.isin(DEGnames).astype('int')
    outframe = pd.concat([allgenes,genescores], axis=1)
    outframe.columns = ['gene', 'DEGstatus']
    outframe.to_csv("./svc_table.csv",sep =',',index=False)

    # Plot RFECV performance across feature-selection steps
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
    plt.savefig("svc_fullsynth_fs_depth" + str($depth) +
                "_refnum_" + str($refnum) +
                "_altnum_" + str($altnum) +
                "_replicate_" + str($x) + ".png", format="png")
    plt.close("all")
    
    print("test")
    """
}
