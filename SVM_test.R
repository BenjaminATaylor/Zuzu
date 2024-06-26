library(tidyverse)
library(e1071)
library(DESeq2)
library(pracma) # for movavg

setwd("/Users/benjamin/Repositories/Zuzu/")

# Do we want to use pre-scaled data to match Python's scaling?
prenormalized = F
informativefeatures = 200

# Read in data (either pre-scaled or 'raw')
data.y = read.csv(file = "test_y.csv", col.names = c("ID","phenotype"))
if(prenormalized){
  data.x = read.csv(file = "test_x_transform.csv", row.names = 1) %>% t()
} else {
  data.x = read.csv(file = "test_x.csv", row.names = 1) %>% t()
}

# ALternatively, read in genesynth data for testing
# data.y = read.csv("synthsheet.csv", col.names = c("ID","phenotype")) %>% 
#   mutate(phenotype = as.numeric(as.factor(phenotype))-1)
# data.x = read.csv("synthcounts.csv", row.names = 1)


#Alternatively, read tumor data for testing
data.y = read.csv("input/test4/tumorcleandata.csv", col.names = c("ID","phenotype")) %>%
  mutate(phenotype = as.numeric(as.factor(phenotype))-1)
data.x = read.csv("input/test4/tumorcleancounts.csv", row.names = 1)


# NB here we can visually check that the first X features are the meaningful ones, as expected
data.x0 = data.x[,which(data.y$phenotype==0)]
data.x1 = data.x[,which(data.y$phenotype==1)]
foo = rowSums(data.x0) - rowSums(data.x1)
plot(foo) # for the gene expression data, I guess part of the issue must be the overdispersion of certain genes. If that's the case, it'd be worth seeing what this looks like when scaled and taking means
data.x.scale = t(apply(data.x,1,scale)) %>% `colnames<-`(colnames(data.x)) 
data.x0 = data.x.scale[,which(data.y$phenotype==0)]
data.x1 = data.x.scale[,which(data.y$phenotype==1)]
foo = abs(rowMeans(data.x0) - rowMeans(data.x1))
plot(foo)
data.x = data.x.scale ; prenormalized = T
# Oookay so these values are higher than expected given the parameters we used, I think?

# Note which features are the informative ones
truefeatures = row.names(data.x)[1:informativefeatures]

# data for manual testing
# readcounts = data.x
# traindata = data.y
# testdata = NA
# crossfold = 5
# referencelevel = 0
# kerneltype = "radial"

#### Define SVM function
svm.train = function(readcounts, traindata, testdata = NA, referencelevel = 0, kerneltype = "radial", crossfold = 5, vstCheck = T){
  
  svm.counts.test=NA
  
  # normalise data
  svm.counts = readcounts
  ## normalize counts between samples 
  # svm.counts.vst.quantiles = normalizeBetweenArrays(svm.counts.vst, method = "quantile")
  # scale counts and remove zero-variance features
  # svm.counts.vst.quantiles.scale = t(scale(t(svm.counts.vst.quantiles)))
  # svm.counts.vst.quantiles.scale = na.omit(svm.counts.vst.quantiles.scale)
  
  # Divide transcriptomic data into training set (queens and workers from control) and test set (individuals from treatment) 
  svm.counts.train = svm.counts[,which(colnames(svm.counts) %in% traindata$ID)]
  if(length(testdata)>1){
    svm.counts.test = svm.counts[,which((colnames(svm.counts) %in% testdata$ID))]
  }
  
  # Perform a grid search to optimise SVM parameters
  svm.counts.tuneResult = tune("svm", 
                               train.x = t(svm.counts.train),
                               train.y = as.numeric(traindata$phenotype == referencelevel),
                               probability = TRUE, 
                               scale = !(prenormalized),
                               kernel = kerneltype,
                               tunecontrol = tune.control(sampling = "cross",
                                                          cross = crossfold),
                               ranges = list(gamma = 10^(-6:-6),
                                             cost = 2^(4:4))
  )
  
  # Final classifier
  svm.counts.classifier = svm.counts.tuneResult$best.model
  
  svm.counts.prediction = NULL
  if(length(testdata)>1){
    # Make predictions for the test data, if test data were provided.
    svm.counts.prediction = predict(svm.counts.classifier,
                                    t(svm.counts.test),
                                    type = "class", 
                                    probability = TRUE)
  }
  
  #output prediction for test data and cross-validation error for training data
  svm.result = list("prediction" = svm.counts.prediction,
                    "validation_error" = signif(svm.counts.tuneResult$best.performance,4),
                    "traincounts" = svm.counts.train,
                    "testcounts" = svm.counts.test)
  
  #return results
  return(svm.result)
}

#### Perform initial classification
# Divide data into training (control) set and test (queen removal) set


# apply svm to entire set of genes
svm.full = svm.train(data.x,
                     data.y, 
                     crossfold = 5)
print(paste0("Root mean cross-validation error rate for full model: ",svm.full$validation_error))









#### Perform feature selection
# create copy of training data that we can subject to repeated trimming while preserving original frame
svm.counts.train.iterate = svm.full$traincounts
#record original number of features
nfeatures = nrow(svm.counts.train.iterate)
#target number of features 
nfeatures_target = 5
traindata = data.y
#instantiate data frame to hold data on the error of each model
iterations = data.frame(feature = character(),
                        error_before_removal = numeric())
#iteratively remove features until target number is reached
while(nfeatures > nfeatures_target){
  
  error = c()
  
  #run repeatedly to account for stochasticity in cross-validation
  for(i in 1:5){
    
    # Perform a grid search to optimise SVM parameters
    svm.counts.tuneResult = tune("svm", 
                                 train.x = t(svm.counts.train.iterate), 
                                 train.y =  as.numeric(traindata$phenotype == 0),
                                 probability = TRUE, 
                                 scale = FALSE,
                                 kernel = "radial", 
                                 tunecontrol = tune.control(sampling = "cross", 
                                                            cross = 3),
                                 ranges = list(gamma = 10^(-7:-6), cost = 2^(5:6)))
    #record error
    error = c(error, svm.counts.tuneResult$best.performance)
  }
  #sample classifier
  svm.counts.classifier = svm.counts.tuneResult$best.model
  #return mean error value
  error = signif(mean(error),4)
  #extract feature weights
  weights = (t(svm.counts.classifier$coefs) %*% svm.counts.classifier$SV)
  #calculate feature with lowest weight (for ties, choose arbitrarily)
  weakfeature = colnames(weights)[which(abs(weights) == min(abs(weights)))[1]]
  #remove lowest-weight feature from data frame
  svm.counts.train.iterate = subset(svm.counts.train.iterate, !(rownames(svm.counts.train.iterate) %in% c(weakfeature)))
  #in a dataframe, store removed feature name and error value before removing that feature
  iterations = rbind(iterations, tibble(feature = weakfeature,
                                        error_before_removal = error))
  #tick down
  nfeatures = (nfeatures-1)
  #output every 20 runs to track progress
  if((nfeatures/20)%%1==0){print(paste0("Features remaining: ",nfeatures))}
}
iterLength = 1:nrow(iterations)
# take moving average to smooth out variation, if desired
use_movavg = F
if(use_movavg){moving_avg = movavg(iterations$error_before_removal, 3, "s") }

# plot data to ensure we have the expected 'hockeystick' shape 
error = if(use_movavg){ moving_avg } else { iterations$error_before_removal }

hockeyData = data.frame(num = iterLength, error = error, truefeature = (iterations$feature %in% truefeatures))
hockeyData_plot = hockeyData
hockeyData_plot$num = abs(iterLength - (max(iterLength)+1))

ggplot(hockeyData_plot, aes(x = num, y = error, colour = truefeature)) +
  geom_point() +
  scale_x_reverse() +
  scale_colour_manual(values = c("Black","Red")) +
  theme_bw()

plot(hockeyData_plot$num,hockeyData_plot$error,colour = hockeyData_plot$feature,
     xlim = rev(c(0, length(hockeyData_plot$error)+5)))





outframe.columns = ['gene', 'DEGstatus']
outframe.to_csv("./svc_table.csv",sep =',',index=False)


# get minimum of this curve to find the point at which the error window is at its minimum
optimal_removal = which(moving_avg == min(moving_avg));
# list the features to be removed from the original set of genes
features_to_remove = iterations$feature[1:optimal_removal]
# new dataframe with less-useful features removed
counts_clean_subsample = subset(counts_clean_subsample,
                                !(rownames(counts_clean_subsample) %in% features_to_remove))
# re-perform support vector classification using the new, optimally caste-separating set of features
svm.optimal = svm.train(counts_clean_subsample, 
                        referencelevel = "queen",
                        svm.data.train, 
                        svm.data.test,
                        crossfold = 3,
                        vstCheck = F)

print(paste0("Number of genes included in optimised model: ", nrow(counts_clean_subsample)))
print(paste0("Root mean cross-validation error rate for optimised model: ", svm.optimal$validation_error))


