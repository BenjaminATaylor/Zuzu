process SVM_DATA_QUASI{
  
  //publishDir "$params.outdir/svmtables", pattern: '*.png'
  
  input: 
  tuple path(samplesheet), path(countsframe)
  tuple val(x), val(refnum), val(altnum)

  output:
  tuple path("countsframe_quasi.csv"), path("trueDEGs.RData"), path(samplesheet), val(refnum)
  //path "svm_quasi_fs.png"

  script:
  """
  #!/usr/bin/env Rscript
  library(tidyverse)
  library(e1071)
  library(pracma) # for movavg
    
  # Read in raw data 
  data.y = read.csv(file = "$samplesheet")
  data.x = read.csv(file = "$countsframe", row.names = 1, check.names = F)
  
  # Pare for testing
  data.x = head(data.x, 7000)
  if("$params.pareml" == "true"){data.x = head(data.x,100) }

  set.seed($x)
  ref.num = $refnum
  alt.num = $altnum

  # Take a random subset of the samples based on our chosen subsample size
  phenos = relevel(factor(unique(data.y\$phenotype)), ref = "$params.reflevel")
  chosen.samples = c(sample(subset(data.y, phenotype == levels(phenos)[1])\$sample,ref.num),
                     sample(subset(data.y, phenotype == levels(phenos)[2])\$sample,alt.num))
  data.y.sub = subset(data.y, sample %in% chosen.samples)
  data.x.sub = data.x[,data.y.sub\$sample]
  stopifnot(colnames(data.x.sub) == data.y.sub\$sample)

  
  # Reencode the phenotype data
  data.y.sub\$phenotype = as.numeric(as.factor(data.y.sub\$phenotype))-1
  
  #testing
  readcounts = data.x.sub
  traindata = data.y.sub
  referencelevel = 0
  kerneltype = "linear"
  crossfold = 5
  
  #### Define SVM function
  svm.train = function(readcounts, traindata, referencelevel = 0, kerneltype = "linear", crossfold = 5){
    
    #
    svm.counts = readcounts
  
    # Perform a grid search to optimise SVM parameters
    svm.counts.tuneResult = tune("svm", 
                                 train.x = t(readcounts),
                                 train.y = as.numeric(traindata\$phenotype == referencelevel),
                                 probability = TRUE, 
                                 scale = TRUE,
                                 kernel = kerneltype,
                                 tunecontrol = tune.control(sampling = "cross",
                                                            cross = crossfold),
                                 ranges = list(gamma = 10^(-6:-6),
                                               cost = 2^(4:4))
    )
    
    # Final classifier
    svm.counts.classifier = svm.counts.tuneResult\$best.model
  
    #output prediction for test data and cross-validation error for training data
    svm.result = list("validation_error" = signif(svm.counts.tuneResult\$best.performance,4),
                      "traincounts" = svm.counts)
    
    #return results
    return(svm.result)
  }
  
  # apply svm to entire set of genes
  svm.full = svm.train(data.x.sub,
                       data.y.sub, 
                       crossfold = 5)
  print(paste0("Root mean cross-validation error rate for full model: ",svm.full\$validation_error))
  
  #### Perform feature selection
  # create copy of training data that we can subject to repeated trimming while preserving original frame
  svm.counts.train.iterate = svm.full\$traincounts
  #record original number of features
  nfeatures = nrow(svm.counts.train.iterate)
  #target number of features 
  nfeatures_target = 5
  traindata = data.y.sub
  #instantiate data frame to hold data on the error of each model
  iterations = data.frame(feature = character(),
                          error_before_removal = numeric())
  # choose step size for feature eleimination
  stepsize = 5
  
  #iteratively remove features until target number is reached
  while(nfeatures >= (nfeatures_target+stepsize)){
    
    error = c()
    
    #run repeatedly to account for stochasticity in cross-validation
    for(i in 1:5){
      
      # Perform a grid search to optimise SVM parameters
      svm.counts.tuneResult = tune("svm", 
                                   train.x = t(svm.counts.train.iterate), 
                                   train.y =  as.numeric(traindata\$phenotype == 0),
                                   probability = TRUE, 
                                   scale = TRUE,
                                   kernel = "radial", 
                                   tunecontrol = tune.control(sampling = "cross", 
                                                              cross = 3),
                                   ranges = list(gamma = 10^(-6:-6), cost = 2^(4:4)))
      #record error
      error = c(error, svm.counts.tuneResult\$best.performance)
    }
    #sample classifier
    svm.counts.classifier = svm.counts.tuneResult\$best.model
    #return mean error value
    error = signif(mean(error),4)
    #extract feature weights
    weights = (t(svm.counts.classifier\$coefs) %*% svm.counts.classifier\$SV)
    #calculate feature, or features, with lowest weight (for ties, choose arbitrarily)
    weakfeatures = colnames(weights)[which(order(abs(weights)) %in% c(1:stepsize))]
    #remove lowest-weight feature(s) from data frame
    weakindices = which(row.names(svm.counts.train.iterate) %in% weakfeatures)
    svm.counts.train.iterate = svm.counts.train.iterate[-weakindices,]
    
    data.frame(feature = weakfeatures, error_before_removal = error)
    
    #in a dataframe, store removed feature name and error value before removing that feature
    iterations = rbind(iterations, data.frame(feature = weakfeatures,
                                          error_before_removal = error))
    #tick down
    nfeatures = (nfeatures-stepsize)
    #output every 20 runs to track progress. TODO: fix this to work with any initial number of starting features. Atm it only works if initial number - stepsize sometimes converges with a number divisible by 20.
    if((nfeatures/20)%%1==0){print(paste0("Features remaining: ",nfeatures))}
  }
  
  write.csv(iterations, "iterations.csv")
  
  errors = unique(iterations\$error_before_removal)
  
  iterLength = 1:length(errors)
  # take moving average to smooth out variation, if desired
  use_movavg = T
  if(use_movavg){moving_avg = movavg(errors, 5, "s") }
  
  # plot data to ensure we have the expected 'hockeystick' shape 
  errors = if(use_movavg){ moving_avg } else { errors }
  
  hockeyData = data.frame(num = iterLength, error = errors)
  hockeyData_plot = hockeyData
  hockeyData_plot\$num = abs(iterLength - (max(iterLength)+1))
  
  hockeyplot = ggplot(hockeyData_plot, aes(x = num, y = error)) +
    geom_point() +
    scale_x_reverse() +
    theme_bw()
    
  #get index of lowest error
  # TODO: this is currently a very crude/imprecise way of obtaining this index. This MUST be fixed
  minerrindex = which(errors == min(errors))
  dropgenes = iterations[1:(minerrindex*stepsize),]\$feature
  
  outframe = data.frame(gene = row.names(svm.full\$traincounts), 
                        DEGstatus = !(row.names(svm.full\$traincounts) %in% dropgenes))
  
  write.csv(outframe, row.names = F, file = "svm_table.csv")

    
  ggsave(hockeyplot,
         filename = "svm_quasi_fs.png",
         device = "png", path = "./",
         width =  38, height = 23, units = "cm")
   
   # identify the top X (100?) genes in the SVM, to be defined as 'true' DEGs for the quasipermutation
   # TODO: find a more nuanced way of defining this for FS
   truedegs = colnames(weights)[which(order(abs(weights)) %in% c(1:100))]
   
   # Select half of true degs
   keepnum = round(length(truedegs)/2)
   
   #Select a subset of the 'true' degs to keep constant
   keepdegs = sample(truedegs, size = keepnum, replace = FALSE)
   
   # Permute everything...
   permuterow = function(x){ sample(x, size = ncol(data.x)) }
   data.x.quasi = apply(X = data.x, MARGIN = 1, FUN = permuterow) %>% 
     t() %>% 
     data.frame() %>% 
     `colnames<-`(colnames(data.x))
   # Then restore true counts for the 'true' DEGs
   data.x.quasi[keepdegs,] = data.x[keepdegs,]
   
   # Check that this has worked as expected: 
   # which(rowSums(data.x == data.x.quasi) == ncol(data.x))
   
   # Output new quasi-permuted frame *and* the set of true degs associated with that frame
   write.csv(data.x.quasi,"countsframe_quasi.csv")
   save(keepdegs, file="trueDEGs.RData")        
  """
}
