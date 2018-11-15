


library(ROCR)
library(glmnet)
library(nnet)
rhome="."
source(paste(rhome, "/fspls.R",sep="")) ## VARIABLE SELECTION CODE

options("family"="binomial") ## for linear regression

options("lambda" = NULL)    ## if 0, no shrinkage. IF specified, then uses that shrinkage, if NULL, then uses CV
options("lambda1" = 1)      ## for pvalue adjustment
options("debug" = "on")     ## a lot of debugging information printed

family = getOption("family","gaussian")


data = NULL
dataTest = NULL


## gaussian dataset
elapse_start = proc.time()


##### be careful about importing data
## 'data' for running the program should be a list with two fields 'data' and 'y' 
## representing the input data matrix and output, respectively
# datf = "/shares/common/groups/Group-Coin/c.zhou/fspls/data/hugeData.RData"
# load(datf)

## this may not be applicable to your data file
## when reading large data file 'fread' function from library 'data.table' could be a better choice
datf = "data.csv" #paste(group_coin, "Users/c.zhou/fspls/data/data.csv"
t = read.csv(datf, sep=",", header=T)
datasets = splitData(t)
data = datasets$train
dataTest = datasets$test



print("LOADING DATA ELPASED TIME:")
print(proc.time()-elapse_start)

## p-value threshold will determine the number of selected features
## you may want to choose a much smaller p-value
pv_thresh = 0.05

options("method" = "fspls")

elapse_start = proc.time()
model_fspls = trainModel(data,max=dim(data$data)[2],pv_thresh=pv_thresh,test=dataTest,refit=FALSE)
print("OFS ELPASED TIME:")
print(proc.time()-elapse_start)

print("VARIABLES SELECTED:")
vars = model_fspls$variables
spvals = model_fspls$spvals
evals = rbind(vars, spvals)
colnames(evals) = colnames(data$data)[vars]
rownames(evals) = c("Variables","p-Values")
print(evals)

save(model_fspls, file=paste0("model_fspls_",gsub("-|:| ","_",Sys.time()),".RData"))
