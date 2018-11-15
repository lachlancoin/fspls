
##can be used to refit a new offset for a set of vars
#return coefficients
refit_model<-function(data,variables,beta,const){
    family = family=getOption("family","binomial")
    if(length(variables)>1){
        x = as.matrix(data$data[,variables])
        if(family!="multinomial"){
	        lm<-glm(data$y~x,family=family)
            coeff = coef(lm)
            const = as.matrix(coeff[1])
            beta = as.matrix(coeff[-1])
        }else{
            lm<-cv.glmnet(x,data$y,family=family)
            coeff = coef(lm,s='lambda.min')
            for(i in 2:length(coeff)){
                const[i-1] = coeff[[i]][1]-coeff[[1]][1]
                beta[,i-1] = coeff[[i]][-1]-coeff[[1]][-1]
            }
        }

    }
    list(beta=beta,const=const)
}

##Aim is to project out all variations corresponding to vars variables (this is just indices)
## from the matrix Dall.
projOut<-function(Dall,vars){
  P = Dall[,vars,drop=F]
  alias = apply(P,1,.cntNA)==0
  P = P[alias,,drop=F]
  ncol = dim(Dall)[[2]]
  nrow = dim(Dall)[[1]]
  svd = svd(P)
  U = svd$u
  V = t(svd$v)
  if(length(vars)==1) {
      D = as.matrix(svd$d)
  }else {
      D = diag(svd$d)
  }
  Vinv = solve(V)
  Dinv = solve(D)
  colU = dim(U)[[2]]
  proj = matrix(0,nrow = dim(U)[[2]],ncol = ncol)
  for(j in 1:ncol){ 	
	 for(k in 1:colU){
           proj[k,j] = U[,k] %*% Dall[alias,j]
	}	
  }
  Vinv %*% Dinv %*% proj
}

##this just calculates the constant term
getOffset<-function(data){
 family=getOption("family","binomial")
 yTr = data$y
  if(family!="multinomial"){
   model <- glm(yTr~rep(1,length(yTr)),family=family) 
   coeffNew = summary(model)$coeff[1,1]
  }else{
   invisible(capture.output({model <- multinom(yTr~rep(1,length(yTr)))}))
   coeffNew <- summary(model)$coeff[,1]
  }
  coeffNew
}

.centralise <-
function(vec,na.replace = NA){
  vec1 = (vec-mean(vec,na.rm=TRUE))
  vec1[is.na(vec1)] = na.replace
  vec1
}

.cntNA<-function(vec) length(which(is.na(vec)))

 fitModel<-function(yTr,x,offset, indices=attr(yTr,"indices"),family=getOption("family","binomial"),  useBF = getOption("useBF",FALSE)){
    ll=0
    n = dim(x)[2]
    if(family!="multinomial"){
	    if(n==0) lm<-glm(yTr~offset(offset),family=family)
	    else lm<-glm(yTr~x+offset(offset),family=family)
    }else{
        offset = cbind(rep(0,dim(offset)[1]),offset)
        if(n==0) invisible(capture.output({lm<-multinom(yTr~offset(offset))}))
        else invisible(capture.output({lm<-multinom(yTr~x+offset(offset))}))
    }
    if(useBF){
        coe = summary(lm)$coeff
        ll =ll+abf(coe[2,])
    }else{
        ll=ll+logLik(lm)
    }
    ll
}

 fitModelShrink<-function(yTr,x,offset,family=getOption("family","binomial"),  useBF = getOption("useBF",FALSE), lambda =0){
    #if((family=="binomial" || family=="multinomial") && getOption("method","fspls")=="fspls") alpha1 = 0
    #else alpha1 = 1
    alpha = getOption("alpha",0)
    ones = rep(1,length(yTr))
    x2 = cbind(x,ones)
    y2 = yTr
    if(is.null(lambda)){
        aa = cv.glmnet(x2,y2,family=family, alpha = alpha, offset=offset)
        aamin = which(aa$cvm==min(aa$cvm))
        indsa = which(aa$cvm>=aa$cvlo[aamin] & aa$cvm <=aa$cvup[aamin])
        lambda =  max(aa$lambda[aamin])
    }
    if(family!="multinomial"){
        ridge = glmnet(x2,y2,family=family,lambda = lambda, alpha = alpha, offset=offset, intercept = F)
        rbeta <- coef(ridge,s="lambda.1se")
        coeffNew = as.matrix(rbeta[-length(rbeta)])
    }else{
        ridge = glmnet(x2,y2,family=family,lambda = lambda, alpha = alpha, offset=offset, intercept = F)
        rbeta <- coef(ridge,s="lambda.1se")
        rbeta = matrix(unlist(lapply(rbeta,as.matrix)),nrow=dim(x)[2]+2)
        rbeta = rbeta-rbeta[,1]
        coeffNew = as.matrix(rbeta[-(dim(x)[2]+2),-1])
    }
    coeffNew
}

#This finds the next best variable to select.
##data is an object with data$y and data$data. Rows are different individuals, and columns diff variables 
#variables is a list of variables already selected
#beta are the weights for those selected variables
#returns updated beta
find_best_glm<-function(data, variables, beta, Wall, const_term, var_thresh = 0.001, missThresh = 0.5, fit_coeff = TRUE, refit = FALSE, best_i1 = NULL, excl = c()){
   family=getOption("family","binomial") 
   yTr = data$y
   leny = length(yTr)
   LI = max(1,length(attr(yTr,"levs"))-1)
   pls = (getOption("method","fspls")=="fspls")
   Dall = as.matrix(data$data)
   #offset = .rep(const_term,leny)
   if(fit_coeff || family=="gaussian") offset = .rep(const_term,leny)
   else offset = .rep(rep(0,LI),leny)
   W = matrix(0,nrow=length(variables),ncol = dim(Dall)[[2]])
   todoInds = 1:dim(Dall)[[2]]
   if(length(variables>=1)){
     todoInds = todoInds[-variables]
     if(pls){
        W = projOut(Dall,variables)
     }else if(!fit_coeff){
        offset = Dall[,variables,drop=F] %*% beta
     }
   }

   R = Dall[,variables,drop=F] %*% W
   spll=-Inf
   best_i=0
   if(pls) ll0 = fitModel(yTr,Dall[,c(),drop=F],offset)
   else  ll0 = fitModel(yTr,Dall[,variables,drop=F],offset)
   start_time = proc.time()
    excl1 = match(excl,todoInds);
   if(length(excl1)>0) todoInds = todoInds[-excl1]
  
   if(!is.null(best_i1)) todoInds = todoInds[which(todoInds==best_i1)];
   for (i1 in todoInds){
	    if(!pls) i = c(variables,i1)
            else i = i1
	    x=Dall[,i,drop=F] - R[,i,drop=F]
        if(var(x[,1],na.rm=T)>var_thresh & .cntNA(x[,1])/leny < missThresh){
          ll = fitModel(yTr,x,offset,indices)
		

		      if( ll>spll ){
			      spll=ll
			      best_i=i
		      }	
		
	    }
  }
  
  print(paste("+find_best_var,",proc.time()[3]-start_time[3]))
  df=1
  if(family=="multinomial") df=length(unique(yTr))-1
  pval=pchisq(2*(spll-ll0),df,lower.tail=F)
  print(paste0(">>pval ",pval))
  print(paste0(">>pval_df1 ",pchisq(2*(spll-ll0),1,lower.tail=F)))
  
  pval = .applySidakCorrection(pval,getOption("lambda1",1))
  if(length(variables)==0){
    Wall = matrix(1)
  }else{
    Wall3 = cbind(Wall,-W[,best_i[length(best_i)]])
    Wall = rbind(Wall3,c(rep(0,ncol(Wall3)-1),1))
  }
  
  variables = c(variables,best_i[length(best_i)])
  coeff =matrix(0,nrow=length(variables),ncol=LI)
  lambda = getOption("lambda",NULL)

  if(!refit){
    x = Dall[,best_i,drop=F] - R[,best_i,drop=F]
   
    start_time = proc.time()
    if(fit_coeff){
        coeffNew = fitModelShrink(yTr,x,offset,lambda=lambda)
    }else{
        if(family!="multinomial"){
            model<-glm(yTr~x,family=family,offset=offset)
        
            coeffNew = summary(model)$coeff[,1,drop=F]
            delta = coeffNew[dim(x)[2]+1]
            beta = beta-delta*W[,best_i[length(best_i)]]
        }else{
            invisible(capture.output({model<-multinom(yTr~x,offset=cbind(rep(0,dim(offset)[1]),offset))}))
            coeffNew = t(summary(model)$coeff)
            if(dim(beta)[1]>0){
                delta = .rep(coeffNew[dim(x)[2]+1,],dim(beta)[1])*t(.rep(W[,best_i[length(best_i)]],dim(beta)[2]))
                beta = beta-delta
            }
          }
    }
    print(paste("+find_new_coeff,",proc.time()[3]-start_time[3]))

    const_term = const_term+coeffNew[1,]

    coeffNew = coeffNew[-1,,drop=F]

    if(!pls) coeff=coeffNew
    else coeff=rbind(beta,coeffNew[1,])
  }
  list(pval=pval, variable=variables, coeff=coeff, var = var(Dall[,best_i[length(best_i)]]), Wall = Wall, const_term=const_term, lambda = lambda)
}


.rep<-function(v,n) t(matrix(rep(v,n),ncol=n))

.applySidakCorrection<-function(pvs,effPhe){
  pv = min(pvs,na.rm=T)
if(getOption("mPhen.log10p",FALSE)){
   if(pv<1e-15) res = pv+log10(effPhe)
   else res = 1-(1-(10^pv))^effPhe
}else{
  if(pv<1e-15) res = pv*effPhe
  else res = 1-(1-pv)^effPhe
}
res
}
.nyholdtSidak<-function(pheno){
M = dim(pheno)[[2]]
if(M==1){
  return(1)
}else{
cor = cor(pheno,use="pairwise.complete.obs")
excl_ind = (apply(apply(cor,c(1,2),is.na),1,max)>0)
cor1 = cor[!excl_ind,!excl_ind]
if(dim(cor1)[1]==0) return(1)
v =var(eigen(cor1)$val)
return((M-1)*(1-v/M)+1)
}
}


shrink<-function(coe,lambda){
  #metagen(c(0,coe[1]),c((1/(log(lambda)/log(2)))*coe[1],coe[2]))$TE.fixed
  p = .applySidakCorrection(coe[4],lambda)
  effSize = qnorm(p/2,sd =coe[2],lower.tail=FALSE)
  effSize*sign(coe[1])
}

tikhonov<-function(A,b,gamma){
AA = solve((t(A) %*%A)+ (t(gamma) %*% gamma))
AA %*% (t(A) %*% b)
}

abf <- function( fit, w=.21*.21 ){
#  fit <- coeff[2,]
  v <- fit[2]
  z <- fit[3]
  abf <- as.numeric(sqrt(v/(w+v)) * exp(0.5*z*z*w/(v+w)))
  return(abf)
}
.prob<-function(vec)vec/sum(vec)
.samp<-function(pr, levs){
  ind = which(pr==max(pr,na.rm=T))
  return(levs[ind[1]])
  #sample(levs,1,pr=pr)
  ##levs[which(pr==max(pr))[1]] 
}


### this calculates the predicted linear values
pred<-function(coeff, vars,data, Wall = diag(rep(1,length(vars))), const = 0, 
	means = if(length(vars)==0) c(0) else apply(data$data[,vars,drop=F],2,mean,na.rm=T)){
   testD = data$data1
   if(is.null(testD)) testD=data$data
   len = dim(testD)[1]
   xM = const
   if(length(vars)>0) {
  	   xM = .rep(xM - means %*% (Wall %*% coeff),len)
       #xM = rep(xM,len)
       xM = xM + (as.matrix(testD[,vars])%*% Wall %*%  coeff)
   }else{
 	   xM = matrix(NA,ncol = length(const),nrow = len)
	   for(k1 in 1:len){
		 xM[k1,] = const
            }
      }
   attr(xM,"levs") = attr(data$y,"levs")
   xM
}

##this takes linear estimates of risk and applies threshold model
liability<-function(xM){
 levs = attr(xM,"levs")
 alias = apply(xM,1,.cntNA)==0
 M = xM[alias,,drop=F]
 len =dim(M)[1]
 M1 = cbind(rep(1,len),exp(M))
 dimnames(M1)[[2]] = levs
 su = apply(M1,1,sum)
 pr = t(apply(M1,1,.prob))
 res = rep(NA,length(alias))
 res[alias] = apply(pr,1,.samp,levs)
 attr(xM,"levs") = levs
 res
}


countCorrect<-function(pr, thresh,y){
  yfact =  as.factor(y)
  levs = levels(yfact)
  
   res = array(0,dim=c(3,length(levs)), dimnames = list(c("TPR","FPR","missing"),levs))
   levs1 = c(levs,"unclassified")
   res1 = array(NA,dim = c(length(levs),length(levs1)),dimnames = list(levs,levs1))

  for(i in 1:length(levs)){
   
     indices = which(yfact==levs[i] & pr[[2]][,1]>=thresh)  
     res[3,i] = 1-length(indices)/length(  which(yfact==levs[i])) 
     v = pr[[1]][indices,1]    ## all true
     positives = which(v==levs[i])  ## true positives

     for(j in 1:length(levs)){
       res1[i,j] = length(which(v==levs[j]))
     }
     res1[i,1+length(levs)] = length(which(yfact==levs[i] & pr[[2]][,1]<thresh) ) 

   
   
     res[1,i] = length(positives)/length(v)

     indices = which(yfact!=levs[i] & pr[[2]][,1]>=thresh)  ## all false
     v = pr[[1]][indices,1]    ## all false
     positives = which(v==levs[i])  ## false positives
     res[2,i] = length(positives)/length(v)
    
  }
   list(summary = apply(res,c(1,2),round,3), confusion_table=res1)
}

##This function takes a data matrix with names and turns into the list object we use for training
makeTrainObject<-function(data, incl=NULL, merge=NULL){
  y =(as.matrix(data.frame(strsplit(dimnames(data)[[1]],"\\.")))[1,])
  if(!is.null(incl)){
     indices = c()
     for(i in incl){
        indices = c(indices,which(y==i))
     }
     data = data[indices,,drop=F]
     y = y[indices]
  
  } 
  if(!is.null(merge)){
    nme = names(merge)
    for(k in 1:length(merge)){
      indices = c()
      for(i in merge[[k]]){
        indices = c(indices,which(y==i))
     }
	y[indices] = nme[k]
    }
  }
  list(data=data,y = y)
}

##this function converts multinomial data to binary, and returns 
## new y with the indices of different binomial comparison to a pivot.
## if continuous, then just returns the original  y
convertToBinary<-function(y, pivot = NULL, levs = levels(as.factor(y))){
   
    if(getOption("family","binomial")=="binomial" || getOption("family","binomial")=="multinomial"){
     yfact = as.factor(y)
     ncol =length(levs)-1
     if(!is.null(pivot)){
       kk = which(levs==pivot)
       levs1 = c(levs[kk],levs[-kk])
       levs = levs1
     }
     LI= length(levs)-1
     indices = list()
     for(k in 1:LI){
       indices[[k]] = which(yfact==levs[1]| yfact==levs[k+1]) 
     }
     y[yfact==levs[1]]=0
     for(k in 2:length(levs)){
       y[yfact==levs[k]]=1
     }
     y = as.numeric(y)
   }else{
     indices = list(1:length(y))
     names(indices) = NULL
     levs = NULL
     ncol = 1
  }
  
  attr(y,"indices") = indices
  attr(y,"levs") = levs
  y
}

##this just pre-processes data
preprocess<-function(trainOriginal, centralise = TRUE, pivot = 1){
  if(is.null(trainOriginal)) return(NULL)
  dat = trainOriginal$data
  if(centralise){
    means = apply(dat,2,mean,na.rm=T) 
    data = apply(dat,2,.centralise)
  }else{
    means = rep(0,dim(dat)[2])
    data = dat
  }
  #train=list(data=data,y=convertToBinary(trainOriginal$y), data1 = dat, y1 = trainOriginal$y)
  train=list(data=data,y=trainOriginal$y, data1 = dat)
  levs=levels(as.factor(trainOriginal$y))
  if(getOption("family","binomial")=="multinomial"){
    train$y <- relevel(as.factor(train$y), ref = which(levs==pivot))
  }  
  attr(train,"means") = means
  if(getOption("family","binomial")!="gaussian"){
    attr(train$y,"levs") = levels(as.factor(train$y))
  }
  train
}

evalModels<-function(coeff,vars,data, Wall, const =0,
                  means = if(length(vars)==0) c(0) else apply(data$data[,vars,drop=F],2,mean,na.rm=T),ind=0, fit_coeff=TRUE){
   if(is.null(data)) return(NULL)
   if(is.null(Wall)) Wall= diag(rep(1,length(vars)))
   if(fit_coeff || getOption("family","binomial")=="gaussian"){
     ypred = pred(coeff,vars,data,Wall,const,means)
   }else{
     testD = data$data1
     if(is.null(testD)) testD=data$data
     ypred = as.matrix(testD[,vars])%*%coeff-as.numeric(means%*%coeff)
     attr(ypred,"levs") = attr(data$y,"levs")
   }
   #ypred = pred(coeff,vars,data,Wall,const,means)
   levs = attr(data$y,"levs")

   auc = NA
   acc = NA
   rms = NA
   if(getOption("family","binomial")=="binomial"){
       ybin = liability(ypred)
	#print('calc AUC');
	#print(data$y);
       auc = calcAUC(ypred,data$y)
       rms = calcRMS(ybin,data$y)
       acc = calcACC(ypred,data$y)
   }else if(getOption("family","binomial")=="multinomial"){
       ybin = liability(ypred)
       ypred0 = as.numeric(ybin)
       acc = sum(ypred0==data$y)/length(ypred0)
   }else{
       rms = calcRMS(ypred,data$y)
   }
   eval = c(ind,rms,auc,acc)
   names(eval) = c("ind","rms","auc","acc")
   eval
}
                    
                    
#train is an object with train$data and train$y  where y is outcome and data is data matrix (cols are variables)
#var_thresh ignores variables with var < var_thresh
#project - do we project out variables we already selected for search
#refit - do we keep the coeff trained for each variable separately, or refit as we go
# you can already specificy variables and beta selected

trainModel<-function(trainOriginal, pv_thresh = 0.01, max = min(25, dim(trainOriginal$data)[[2]]), var_thresh = 0.01, project = TRUE, 
	 info = list(),pivot = 0, fit_coeff = TRUE, testOriginal = NULL, refit = FALSE, log = NULL, append = FALSE, sel_vars = NULL, excl = c()){
  proc_start_time = proc.time()
  pv_thresh1 = 1.0;
  toremove = which(sel_vars %in% excl)
  if(length(toremove)>0) sel_vars = sel_vars[-toremove];
     
  if(!is.null(log)) sink(log,append)
  spval=0.0
  if(max>dim(trainOriginal$data)[[2]]) max = dim(trainOriginal$data)[[2]]
  if(getOption("family","binomial")=="multinomial" && is.null(pivot)) stop("Invalid pivot provide.")
  if(fit_coeff || getOption("family","binomial")=="gaussian"){
    train = preprocess(trainOriginal,pivot = pivot)
  }else{
    train = preprocess(trainOriginal, centralise = FALSE, pivot=pivot)
  }
  if(!is.null(testOriginal)){
	  test = preprocess(testOriginal, centralise = FALSE, pivot=pivot)
 	 attr(test$y,"levs") = attr(train$y,"levs")
  }
  means = attr(train,"means")
  levs = attr(train$y,"levs")
  ncol = max(1,length(levs)-1)
  beta = matrix(ncol=ncol,nrow=0)
  betaGlobal = list()
  Wall =NULL   ## this object is the transformation
  constant_term1 = getOffset(train)
  print(constant_term1)
  variables = c()
  contains = c()
  eval =  evalModels(beta, variables,train,Wall, const = constant_term1, means = means[variables], ind = 0, fit_coeff = fit_coeff)
  if(!is.null(test)){
 	 evalT = evalModels(beta, variables,test,Wall, const = constant_term1, means = means[variables], ind = 0, fit_coeff = fit_coeff)
	print(c(eval,evalT))
  }else{
	print(eval)
	}
  
  constant_term = constant_term1
  global_index=length(variables)+1
  spvals = c()
  lambdas = c()
  while(global_index<=max && spval <=pv_thresh1 && length(contains)<=1){
	  print(paste("index",global_index))
	best_i = NULL;
	if(global_index<=length(sel_vars)) best_i = sel_vars[global_index]
	if(!is.null(best_i)) pv_thresh1 = 1.0 else pv_thresh1 = pv_thresh      

      start_time <- proc.time()
      model_fit=find_best_glm(train, variables, beta, Wall, const_term=constant_term, var_thresh =var_thresh, fit_coeff = fit_coeff, refit=refit, best_i = best_i, excl = excl)
	
      print(paste("&find_best_glm,",proc.time()[3]-start_time[3]))
      lambdas = c(lambdas,model_fit$lambda)
      varia = model_fit$variable
	  lastVar = varia[length(varia)]
      contains = which(varia==lastVar)
	  spval=model_fit$pval[1] # pval of selected model
      if(length(contains) <=1 && spval <= pv_thresh1){
        variables = varia
        spvals = c(spvals,spval)
	    beta = model_fit$coeff
	    Wall = model_fit$Wall
        print(Wall)
        constant_term = model_fit$const_term
        start_time = proc.time()
        if(refit){
            model_refit = refit_model(train,variables,beta,constant_term)
            beta = model_refit$beta
            constant_term = model_refit$const
            Wall = diag(rep(1,length(variables)))
        }
	betaGlobal[[global_index]] = beta
        print(paste("&find_refit_coeff,",proc.time()[3]-start_time[3]))
        cat("\n")
        print(paste("var",paste(variables,collapse=" ")))
        print(paste("beta",paste(beta,collapse=" ")))
        print(paste("const",constant_term))
        print(paste("pvals",paste(spvals,collapse=" ")))
        eval1 = evalModels(beta, variables,train,Wall, const = constant_term, means = means[variables], ind = global_index, fit_coeff = fit_coeff)
	    eval1T1 =evalModels(beta, variables,test,Wall, const = constant_term, means = means[variables], ind = global_index, fit_coeff = fit_coeff)
	    eval = rbind(eval,eval1)
	    evalT = rbind(evalT,eval1T1)
	 	print(cbind(eval,evalT))
	 	global_index=global_index+1
        cat("\n")
      }
    if(length(contains)>1) print("STOPPED DUE TO REPEAT")
  }
  cat("\n@")
  print(paste("&proc_running_time,",proc.time()[3]-proc_start_time[3]))
  if(!is.null(log)) sink()
  mv = means[variables] * beta
  list(variables=variables,beta=beta,betaGlobal = betaGlobal, levels=levs,constantTerm=constant_term,
  means=means[variables],spvals=spvals,Wall=Wall,eval=cbind(eval,evalT),lambdas=lambdas)
}

centralise<-function(data){
 dimd = dim(data)
 print(dimd)
 mean = apply(data,2,mean)
 sd = apply(data,2,sd)
  for(i in 1:dimd[[2]]){
  data[,i] = (data[,i] - rep(mean[i],dimd[[1]]))/sd[i]
 }
 data
}


#test is as for trainModel
#results is results from train model

testModels<-function(path,cutoff,label, model, setOriginal){
	set<-setOriginal
	yfact = as.factor(setOriginal$y)
	levs =  levels(yfact)
	LI= length(levs)-1
	indices = list()
	for(k in 1:LI){
		indices[[k]] = which(yfact==levs[1]| yfact==levs[k+1]) 
	}
	set$y[yfact==levs[1]]=0
	for(k in 2:length(levs)){
		set$y[yfact==levs[k]]=1
	}
	
	print(getwd())
	beta = model$beta             #betas
  	variables = model$variables   #indeces
  	info = model$info 
	
    print(paste("varid","aucTest LTBI", "spval LTBI", "aucTest OD", "spval OD", "variance", "pval model"))
	res=array(0, dim=c(dim(beta)[1],6))	
	for(m in 1:dim(beta)[1]){
		par(new=TRUE)
		aucTest=calcAUC(beta[1:m,], variables[1:m], set, indices)
		
		res[m,]=c(variables[m],aucTest$auc[1],aucTest$model.pval[1],aucTest$auc[2],aucTest$model.pval[2],info[[m]][1])
		
	print(paste(variables[m],aucTest$auc[1],aucTest$model.pval[1],aucTest$auc[2],aucTest$model.pval[2],info[[m]] ,sep="  "))

	}

	k=length(beta)

	#if(length(levels(factor(set$y)))==2){
	#	print(length(levels(factor(set$y))))
	#	aucTest=calcAUC(beta[1:k], variables[1:k], set)
	#	}

	#if(length(levels(factor(test$y)))>2)aucTest=calcAUCs(beta[1:k], variables[1:k], test)
	colnames(res)<-c("varid","aucTest LTBI", "spval LTBI", "aucTest OD", "spval OD", "model" )
	write.table(res, file=paste(path,cutoff,"_",label,"model_test.txt",sep=""), sep="\t")
	res
}


##calcs RMS from yTs and predicted y, from linear model
calcRMS<-function( predy,yTs){
     rms = NA
     if(getOption("family","binomial")!="gaussian"){
        y = predy
        alias = !is.na(y)
        rms = (sqrt(length(which(y[alias]!=yTs[alias]))/length(y[alias])))
     }else{
        y = predy
	    rms[1] = sqrt(sum((yTs - y)^2)/length(yTs))
     }
     rms
}

##this just returns the auc
#vars are indices referring to data$data
calcAUC<-function(ypred,y){
    yTs = as.numeric(y)
    x = as.numeric(ypred)
    auc = NA
    if(var(x, na.rm=T)>=0 ){
        model<-glm(yTs~x,family=getOption("family","binomial")) 
	    predictions=predict(model,x=x,type="response")
        pred <- prediction(predictions, yTs)
	    perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
	    auc <- performance(pred, measure = "auc")@y.values[[1]]
    }
   auc
}

#this just returns the acc
#vars are indices referring to data$data
calcACC<-function(ypred,y){
    p=rep(NA,length(ypred))
    yTs = as.numeric(y)
    x = as.numeric(ypred)
    acc = NA
    if(var(x, na.rm=T)>0 ){
        model<-glm(yTs~x,family=getOption("family","binomial"))
        predictions=predict(model,x=x,type="response")
        pred <- prediction(predictions, yTs)
        pref <- performance(pred,measure = "acc")
        acc <- pref@y.values[[1]][max(which(pref@x.values[[1]]>=.5))]
    }
    acc
}

if(FALSE){

##this is the code that will run it
###EXAMPLE 0 - all 

train = makeTrainObject(data,incl=NULL,merge=NULL)
model = trainModel(train,max=30,pv_thresh=0.01,pivot="LumA")
assignments  = getInvLogitPr(train,model)
countCorrect(assignments,0.0,train$y)

###EXAMPLE 1
train = makeTrainObject(data,incl=c("LumB","LumA","Basal","Her2"),merge=list(Luminal = c("LumA","LumB")))
model = trainModel(train,max=30,pv_thresh=0.01,pivot="Luminal")
assignments  = getInvLogitPr(train,model)
countCorrect(assignments,0.0,train$y)

###EXAMPLE 2
train = makeTrainObject(data,incl=c("LumB","LumA"),merge=NULL)
model = trainModel(train,max=30,pv_thresh=0.01,pivot="LumA")
assignments  = getInvLogitPr(train,model)
countCorrect(assignments,0.0,train$y)


###EXAMPLE 3
train = makeTrainObject(data,incl=NULL,merge=list(cancer=c("LumA","LumB","Basal","Her2")))
model = trainModel(train,max=30,pv_thresh=0.01,pivot="Normal")
assignments  = getInvLogitPr(train,model)
countCorrect(assignments,0.0,train$y)
}
