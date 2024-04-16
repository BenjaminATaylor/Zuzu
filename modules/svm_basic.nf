process SVM_BASIC{

  label 'big_job'
  publishDir "$params.outdir/svmtables", pattern: '*.png'
  //publishDir ".", pattern: 'svc_table.csv'
  
  input: 
  tuple path(samplesheet), path(countsframe) 

  output:
  path "svm_table.csv", emit: table
  path "svm_basic_fs.png"
  path "iterations.csv"

  script:
  """
  #!/usr/bin/env Rscript
  library(tidyverse)
  library(e1071)
  library(pracma) # for movavg
    
  # Read in raw data 
  data.y = read.csv(file = "$samplesheet")
  data.x = read.csv(file = "$countsframe", row.names = 1)
  
  # Pare for testing
  data.x = head(data.x, 7000)
  if("$params.pareml" == "true"){data.x = head(data.x,100) }

  # Reencode the phenotype data
  data.y\$phenotype = as.numeric(as.factor(data.y\$phenotype))-1
  
  #testing
  readcounts = data.x
  traindata = data.y
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
  svm.full = svm.train(data.x,
                       data.y, 
                       crossfold = 5)
  print(paste0("Root mean cross-validation error rate for full model: ",svm.full\$validation_error))
  
  #### Perform feature selection
  # create copy of training data that we can subject to repeated trimming while preserving original frame
  svm.counts.train.iterate = svm.full\$traincounts
  #record original number of features
  nfeatures = nrow(svm.counts.train.iterate)
  #target number of features 
  nfeatures_target = 5
  traindata = data.y
  #instantiate data frame to hold data on the error of each model
  iterations = data.frame(feature = character(),
                          error_before_removal = numeric())
  # choose step size for feature elimination
  stepsize = $params.stepsize
  
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
         filename = "svm_basic_fs.png",
         device = "png", path = "./",
         width =  38, height = 23, units = "cm")
           
  """
}
