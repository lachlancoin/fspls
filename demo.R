library(ROCR)
library(glmnet)
library(nnet)

source("fspls.R") ##VARIABLE SELECTION CODE


## extrat variables selected from a cv.glmnet object
.vars <- function(glmnetObj){
    vars = c()
    coeff = coef(glmnetObj,s='lambda.min')
    if(getOption("family","gaussian")=="multinomial"){
        for(i in 1:length(coeff)) {vars = c(vars,which(coeff[[i]][-1]!=0))}
        vars = unique(vars)
    }else{
        vars = which(coeff[-1]!=0)
    }
    vars
}

## extrat variables and rmse/auc/acc
.eval <- function(obj,print=TRUE){
    if(is.null(obj)) return(NULL)
    family = getOption("family","gaussian")
    res = NULL
    var = NULL
    if(family=="gaussian"){
        eva = "RMSE"
        pos = c(2,6)
    }else if(family=="binomial"){
        eva = "AUC"
        pos = c(3,7)
    }else if(family=="multinomial"){
        eva = "ACC"
        pos = c(4,8)
    }
    tr = paste0("Tr",eva)
    ts = paste0("Ts",eva)
    t = c(tr,ts,"No.VARs")
    if(class(obj)=="list"){
        eval = obj$eval
        var = obj$variables
        res = c(eval[dim(eval)[1],pos],length(var))
        names(res) = t
    }else if(class(obj)=="cv.glmnet"){
        var = .vars(obj)
        res = c(obj$eval,length(var))
        names(res) = t
    }
    if(print){
        print(paste0(deparse(substitute(obj)),":"))
        print(res)
    }
    list(res=res,var=var)
}


options("family"="gaussian") ## for linear regression
#options("family"="binomial") ## for logistic regression
#options("family"="multinomial")  ## for multinomial probit regression
categories = 5

options("lambda" = NULL)    ## if 0, no shrinkage. IF specified, then uses that shrinkage, if NULL, then uses CV
options("lambda1" = 1)   ## for pvalue adjustment

family = getOption("family","gaussian")

pv_thresh = NULL

refit=FALSE

## -- load demo dataset -- ##
if(family=="gaussian"){
    ## gaussian dataset
    datf = "data/z.tbod.gaussian.11.dat.rdat"
    pv_thresh = 0.05
}else if(family=="binomial"){
    ## binomial dataset
    datf = "data/w.tbod.binomial.11.dat.rdat"
    pv_thresh = 0.01
}else if(family=="multinomial"){
    ## 5 categories multinomial dataset
    datf = "data/r.tbod.multinomial.totcat5.11.dat.rdat"
    if(categories==3)
    ## 3 categories multinomial dataset
    datf = "data/r.tbod.multinomial.totcat3.11.dat.rdat"
    pv_thresh = 0.001
    refit = TRUE
}
load(datf)

model_fspls_refit_coeff_t = NULL
model_fspls_refit_coeff_f = NULL
model_fs_refit_coeff_t = NULL
model_cv_lasso = NULL
model_cv_enet = NULL

options("method" = "fspls")
## run new version fspls
model_fspls_refit_coeff_t = trainModel(data,max=30,pv_thresh=pv_thresh,test= dataTest,refit=refit)
## run old version fspls
#model_fspls_refit_coeff_f = trainModel(data,max=30,pv_thresh=pv_thresh,test= dataTest,fit_coeff=FALSE,refit=refit)

options("method" = "fs")
## run fs
model_fs_refit_coeff_t = trainModel(data,max=30,pv_thresh=0.001,test= dataTest,refit=refit)

## run lasso
model_cv_lasso = cv.glmnet(data$data,data$y,alpha=1,family=family)
## run elastic net
model_cv_enet = cv.glmnet(data$data,data$y,alpha=0.5,family=family)

calcACC <- function(ypred,y) sum(ypred==y)/length(y)

#### evaluate lasso and enet
if(family=="gaussian"){ ## evaluated by RMSE
    if(!is.null(model_cv_lasso)){
        ptr_lasso = predict(model_cv_lasso,data$data,type="response")
        pts_lasso = predict(model_cv_lasso,dataTest$data,type="response")
        model_cv_lasso$eval = c(calcRMS(ptr_lasso,data$y),calcRMS(pts_lasso,dataTest$y))
    }
    if(!is.null(model_cv_enet)){
        ptr_enet = predict(model_cv_enet,data$data,type="response")
        pts_enet = predict(model_cv_enet,dataTest$data,type="response")
        model_cv_enet$eval = c(calcRMS(ptr_enet,data$y),calcRMS(pts_enet,dataTest$y))
    }
}else if(family=="binomial"){ ## evaluated by the area under curve(AUC)
    if(!is.null(model_cv_lasso)){
        ptr_lasso = predict(model_cv_lasso,data$data,type="response")
        pts_lasso = predict(model_cv_lasso,dataTest$data,type="response")
        model_cv_lasso$eval = c(calcAUC(ptr_lasso,data$y),calcAUC(pts_lasso,dataTest$y))
    }
    if(!is.null(model_cv_enet)){
        ptr_enet = predict(model_cv_enet,data$data,type="response")
        pts_enet = predict(model_cv_enet,dataTest$data,type="response")
        model_cv_enet$eval = c(calcAUC(ptr_enet,data$y),calcAUC(pts_enet,dataTest$y))
    }
}else if(family=="multinomial"){ ## evaluated by accuracy(ACC)
    if(!is.null(model_cv_lasso)){
        ptr_lasso = predict(model_cv_lasso,data$data,type="class")
        pts_lasso = predict(model_cv_lasso,dataTest$data,type="class")
        model_cv_lasso$eval = c(calcACC(ptr_lasso,data$y),calcACC(pts_lasso,dataTest$y))
    }
    if(!is.null(model_cv_enet)){
        ptr_enet = predict(model_cv_enet,data$data,type="class")
        pts_enet = predict(model_cv_enet,dataTest$data,type="class")
        model_cv_enet$eval = c(calcACC(ptr_enet,data$y),calcACC(pts_enet,dataTest$y))
    }

}

## extract variables and rmse/auc/acc from models
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

res_fspls_t = .eval(model_fspls_refit_coeff_t)
#res_fspls_f = .eval(model_fspls_refit_coeff_f)
res_fs_t = .eval(model_fs_refit_coeff_t)
res_cv_lasso = .eval(model_cv_lasso)
res_cv_enet = .eval(model_cv_enet)

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
