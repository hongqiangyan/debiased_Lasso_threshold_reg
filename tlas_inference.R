# Simulation function for the threshdolded Lasso

mc_ttlas <- function( lpar1,lpar2=NULL,s_0,vobs ,vtau,
                      taugrid=taugrid,iter,
                      ncores,sig_eps,corrXQ,
                      intercept,thd_intercept=TRUE,ic)
{	
  
  tt <- proc.time()
  
  # Names of the stats
  sn <- c('size_th','size_lm','power_th','power_lm','CRnz','CRz','CRnzt','CRzt','CRnzlm','CRzlm')
  # Creating the storage array
  astat <- array(NA, dim=c(length(vobs),length(vtau),length(lpar1),length(sn))
                 ,dimnames = list('nobs'=vobs,
                                  'tau'=vtau,'par'=NULL,'stats'=sn) )
  
  # Looping over the parameter specifications
  pcount <- 0 # specification counter
  for(par1 in lpar1){
    pcount <- pcount + 1
    if(is.null(lpar2))par2 <- par1[-c(1:2)]
    if(!is.null(lpar2))par2 <- lpar2[[pcount]]
    # Looping over the threshold locations:
    for(tau in vtau){
      # Looping over number of observations:
      for(nobs in vobs){
        # Estimate the model
        rt <- proc.time()
        mctlas <- mclapply(1:iter,ittlas, par1,par2,nobs,tau,taugrid=taugrid,
                           sig_eps=sig_eps,corrXQ=corrXQ
                           ,mc.cores=ncores,intercept,thd_intercept,ic)
        # gather the results, compute the sstatistics, store in array
        vpar <- c(par1,par2)
        lstat <- mcstat(mctlas,vpar,tau,sig_eps,s_0,ic)
        cat(paste0('\nn=',nobs,' par:',pcount,' tau:',tau,' corrX_2&Q:',corrXQ,' debiasing:',round(proc.time()-rt,2)[3],' seconds. \n'))
        
        astat[as.character(nobs),as.character(tau),pcount,'size_th'] <- mean(lstat$testh_0)
        astat[as.character(nobs),as.character(tau),pcount,'size_lm'] <- mean(lstat$testh_0lm)
        astat[as.character(nobs),as.character(tau),pcount,'power_th'] <- mean(lstat$testh_a)
        astat[as.character(nobs),as.character(tau),pcount,'power_lm'] <- mean(lstat$testh_alm)
        astat[as.character(nobs),as.character(tau),pcount,'CRnz'] <- mean(lstat$CRnz)
        astat[as.character(nobs),as.character(tau),pcount,'CRnzlm'] <- mean(lstat$CRnzlm)
        astat[as.character(nobs),as.character(tau),pcount,'CRz'] <- mean(lstat$CRz)
        astat[as.character(nobs),as.character(tau),pcount,'CRzlm'] <- mean(lstat$CRzlm)
        astat[as.character(nobs),as.character(tau),pcount,'CRnzt'] <- mean(lstat$CRnzt)
        astat[as.character(nobs),as.character(tau),pcount,'CRzt'] <- mean(lstat$CRzt) 
        
        
      }
    }
  }
  
  cat(paste0('\n\nTotal run time:',round((proc.time()-tt)/60,2)[3],' minutes. \n'))
  
  return(astat)
  
}

# A function taking a MC run, returning an array of stats
mcstat <- function(mctlas,vpar,tau,sig_eps,s_0,ic)
{
  
  #init
  iter <- length(mctlas)
  lmbd <- c()
 CRnzlm<-CRzlm<-CRnz<- CRz<-CRnzt<-CRzt<-testh_0lm<-testh_alm<-testh_0<-testh_a<-matrix(0,1,iter)
  
  #loop over iterations
  for(i in 1:iter)
  {
    mc <- mctlas[[i]]
    # Estimation parameters th
    gmin <- which.min(mc$Vobj)
    lmbd <- c(lmbd,mc$lambda[gmin])
    alphath<- mc$alphath[,gmin]
    npar  <- ncol(mc$X2)
    # debiased lasso th
    # inverse martix
    xltauhat=mc$XG[which(mc$thdvar<mc$taugrid[gmin]),]
    xgtauhat=mc$XG[-which(mc$thdvar<mc$taugrid[gmin]),]
    #scaled the inverse to make it consistent with sample cov
    A=est_ndwcov(xltauhat,ic)[[2]]/nrow(xltauhat)
    B=est_ndwcov(xgtauhat,ic)[[2]]/nrow(xgtauhat)
    thetath=rbind(cbind(B,-B),cbind(-B,A+B))*length(mc$y)	
    dalphath=alphath[-c(1:2)]+thetath%*%t(cbind(mc$XG,mc$X2))%*%(mc$y -alphath[1]-alphath[2]-cbind(mc$XG,mc$X2)%*%alphath[-c(1:2)])/length(mc$y)		
    
    # debiased lasso linear
    # inverse martix
    lmalphals<-mc$lmalphals
    thetalm=est_ndwcov(mc$XG,ic)[[2]]
    dalphalm=lmalphals[-1]+thetalm%*%t(mc$XG)%*%(mc$y-lmalphals[1]- mc$XG%*%lmalphals[-1])/length(mc$y)	
    ######hypothesis
    varth<-thetath%*%t(cbind(mc$XG,mc$X2))%*%(mc$y -alphath[1]-alphath[2]-cbind(mc$XG,mc$X2)%*%alphath[-c(1:2)])%*%
      t(mc$y -alphath[1]-alphath[2]-cbind(mc$XG,mc$X2)%*%alphath[-c(1:2)])%*%cbind(mc$XG,mc$X2)%*%t(thetath)/length(mc$y)	
    # h<- c(1/sqrt(2),rep(0,s_0-1),1/sqrt(2),rep(0,ncol(mc$X2)-s_0-1),rep(0,ncol(mc$X2)))
    # ## H_0 IS TRUE SIZE
    # testh_0[i]<-as.integer(length(mc$y)	*(t(h)%*%(dalphath-vpar[-c(1:2)]))^2/t(h)%*%varth%*% h>qchisq(1-2*sig_eps, df=2))
    # ### H_A IS TRUE POWER
    # htparh0<-c(rep(1,s_0),0.5,rep(0,ncol(mc$X2)-s_0-1),rep(0,ncol(mc$X2)))
    # testh_a[i]<-as.integer(length(mc$y)*(t(h)%*%(dalphath-htparh0))^2/(t(h)%*%varth%*% h)>qchisq(1-2*sig_eps, df=2))
    ### H_0 IS TRUE SIZE
    h<- c(1/2,rep(0,s_0-1),1/2,rep(0,npar	-s_0-1),1/2,rep(0,s_0-1),1/2,rep(0,npar-s_0-1))
    testh_0[i]<-as.integer(length(mc$y)	*(t(h)%*%(dalphath-vpar[-c(1:2)]))^2/(t(h)%*%varth%*% h)>qchisq(1-2*sig_eps, df=4))
    ### H_A IS TRUE POWER
    htparh0<-c(rep(1,s_0),rep(0,npar-s_0),rep(1,s_0),1,rep(0,npar-s_0-1))
    testh_a[i]<-as.integer(length(mc$y)*(t(h)%*%(dalphath-htparh0))^2/(t(h)%*%varth%*% h)>qchisq(1-2*sig_eps, df=4))
    
    #######Coverage rate
   CRnz[i]<- as.integer(abs(dalphath[1]-vpar[3]) <=sqrt(varth[1,1])*qnorm(1-sig_eps)/sqrt(length(mc$y)))
   CRz[i]<-as.integer(abs(dalphath[1+s_0]-vpar[3+s_0]) <=sqrt(varth[1+s_0,1+s_0])*qnorm(1-sig_eps)/sqrt(length(mc$y)))
    CRnzt[i]<-as.integer(abs(dalphath[1+npar]-vpar[3+npar]) <=sqrt(varth[1+npar,1+npar])*qnorm(1-sig_eps)/sqrt(length(mc$y)))
    CRzt[i]<-as.integer(abs(dalphath[1+npar+s_0]-vpar[3+npar+s_0]) <=sqrt(varth[1+npar+s_0,1+npar+s_0])*qnorm(1-sig_eps)/sqrt(length(mc$y)))
  
    ########linear
    varlm<-thetalm%*%t(mc$XG)%*%(mc$y -lmalphals[1]-mc$XG%*%lmalphals[-1])%*%
      t(mc$y -lmalphals[1]-mc$XG%*%lmalphals[-1])%*%mc$XG%*%t(thetalm)/length(mc$y)
    hlm<- c(1/sqrt(2),rep(0,s_0-1),1/sqrt(2),rep(0,npar-s_0-1))
    
    ### H_0 IS TRUE SIZE
    testh_0lm[i]<-as.integer(length(mc$y)*(t(hlm)%*%(dalphalm-vpar[3:(npar+2)]))^2/(t(hlm)%*%varlm%*% hlm)>qchisq(1-2*sig_eps, df=2))
    ### H_A IS TRUE POWER
    htparh0lm<-c(rep(1,s_0),1,rep(0,npar-s_0-1))
    testh_alm[i]<-as.integer(length(mc$y)*(t(hlm)%*%(dalphalm-htparh0lm))^2/(t(hlm)%*%varlm%*% hlm)>qchisq(1-2*sig_eps, df=2))
    #######Coverage rate
    CRnzlm[i]<- as.integer(abs(dalphalm[1]-vpar[3]) <=sqrt(varlm[1,1])*qnorm(1-sig_eps)/sqrt(length(mc$y)))
    CRzlm[i]<-as.integer(abs(dalphalm[1+s_0]-vpar[3+s_0]) <=sqrt(varlm[1+s_0,1+s_0])*qnorm(1-sig_eps)/sqrt(length(mc$y)))
    }
  
  lstat <- list('testh_0lm'=testh_0lm,'testh_alm'=testh_alm,'testh_0'=testh_0,'testh_a'=testh_a,'CRnz'=CRnz,'CRz'=CRz,'CRnzt'=CRnzt,'CRzt'=CRzt,'CRnzlm'=CRnzlm,'CRzlm'=CRzlm)
  
  return(lstat)	
}

# Constructs a variable such that corr(x,y) ~= r
# returns a data frame of two variables which correlate with a population correlation of rho
# If desired, one of both variables can be fixed to an existing variable by specifying x
getBiCop <- function(n, rho, mar.fun=runif, x = NULL, ...) {
  if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")
  df<-matrix(rho, nrow = n, ncol = 2)
  X1 <- mar.fun(n)
  I<-rbinom(n,1, rho)
  X2 <-mar.fun(n)
  X3<-mar.fun(n)
  df[,1] <- I*X1+(1-I)*X2
  df[,2]<- I*X1+(1-I)*X3
  return(df)
}


# generates the data
# calls the ttlas function
# called by the mc_ttlas function
ittlas <- function(i, par1,par2,nobs,tau, taugrid=taugrid,
                   sig_eps=sig_eps,corrXQ=corrXQ,intercept,
                   thd_intercept,ic)
{
  nvars <- length(par1) # nbr vars (incl. cste)
  
  #covariates
  XG <- matrix(rt(nobs*(nvars-2),10),ncol=nvars-2)# generating the regessors
  
  #threshold variable
  if(corrXQ==0) {
    thdvar <- matrix(runif(nobs),ncol=1) # uncorr threshold variable
  }
  else {
    if(corrXQ >1) {
      thdvar <-matrix(runif(nobs),ncol=1)
      XG[,2]<- thdvar# threshold variable is X1
    }
    else {
      df<-getBiCop(n = nobs, rho = corrXQ, mar.fun=runif)
      thdvar <- df[,2]
      XG[,2]<-df[,1]
    }
  }
  
  # Threshold covariates
  
  if(thd_intercept) X1<- cbind(matrix(thdvar<tau,ncol=1,nrow=nobs),XG)		
  if(!thd_intercept) X1<- cbind(0,XG)	
  if(intercept) X1<- cbind(1,X1)		
  if(!intercept) X1<- cbind(0,X1)		
  X2 <- XG
  X2[thdvar>=tau,] <- 0 # threshold regressors
  y <- X1%*%par1 + X2%*%par2 + rt(nobs,10) #generate the data
  ############Threshold
  estth <- ttlas(y , XG=XG, x=NULL , thdvar = thdvar, taugrid=taugrid, intercept=intercept,
                 thd_intercept=thd_intercept, ic=ic, standardize = FALSE)
  estlm<-DSLasso (y ,XG=XG,intercept = intercept,ic=ic, standardize = FALSE)
  
  
  estth$X <- cbind(X1,X2)
  estth$XG <- XG
  estth$X2 <- X2
  estth$y <- y
  
  est<-c(estth,estlm)
  return(est)
}


###########################


# core function to estimate the threshold model by Lasso. 
ttlas <- function(y, XG, x=NULL , thdvar ,
                  taugrid , nfolds=10,
                  intercept, thd_intercept, ic='CV', standardize = FALSE)
{
  #init 
  n <- length(y)	
  # Storage
  lmbdcv <- c()
  alphath<- alp <- NULL
  Vobj <- c()
  
  # Looping over the taus
  for(cnt in 1:length(taugrid)){
    tau <- taugrid[cnt]
    
    # Construction of x tau
    if(thd_intercept) X1<- cbind(matrix(thdvar<tau,ncol=1,nrow=n),XG)	#X1<- cbind(1,XG)
    if(!thd_intercept) X1<- cbind(0,XG)	
    X2 <- XG
    X2[thdvar>=tau,] <- 0
    
    
    # Concatenation
    xall <- cbind(X1,X2)
    if(!is.null(x)) xall <- cbind(xall,x)	
    # Constructing the weights
    D <- apply(xall,2,function(x)sqrt(mean((x-mean(x))^2)))
    if(!is.null(x)) D[length(D)-c(ncol(x):1)+1] <- 0
    
    # cross validation for lambda selection
    if(ic=='CV')cgn	<- cv.glmnet(x = xall, y=y, penalty.factor = D, nfolds=nfolds, standardize = standardize, intercept = intercept)
    if(ic== 'BIC'){cgn	<- glmnet(x = xall, y=y, penalty.factor= D, intercept=intercept, standardize = standardize) 
    yhat    <- predict(cgn, newx=xall, type = "response")
    sigmah	<- colSums((matrix(y,nrow(yhat),ncol(yhat))-yhat)^2)/n
    bic     <- log(sigmah) + cgn$df *log(n)/n  
    cgn$lambda.min 	<- cgn$lambda[which.min(bic)]
    }	
    if(ic=='GIC'){cgn	<- glmnet(x = xall, y=y, penalty.factor= D, intercept=intercept, standardize = standardize) 
    yhat    <- predict(cgn, newx=xall, type = "response")
    sigmah	<- colSums((matrix(y,nrow(yhat),ncol(yhat))-yhat)^2)/n
    gic     <- log(sigmah) + cgn$df *log(log(n))/n * log(ncol(xall)-1) # Fan & Tang JRSS-B 2004
    cgn$lambda.min 	<- cgn$lambda[which.min(gic)] 
    }	

    # Storing
    lmbdcv[cnt] <- cgn$lambda.min
    if(ic=='CV')alp <- matrix(coef(cgn,s='lambda.min'))
    if(ic=='GIC')alp <- c(cgn$a0[which.min(gic)],cgn$beta[,which.min(gic)])
    if(ic=='BIC')alp <- c(cgn$a0[which.min(bic)],cgn$beta[,which.min(bic)])
    #if(intercept & BIC) alp<-c(cgn$a0[which.min(bic)],alp)
    alphath <- cbind(alphath,alp)	
    
    # Computing the value of the objective function:
    if(ic=='CV')yhat	 <- predict(cgn,newx = xall,s='lambda.min')
    RSS 	 <- sum((matrix(y,nrow(yhat),ncol(yhat))-yhat)^2)/(2*n)
    Vobj[cnt]<- RSS + cgn$lambda.min*sum(abs(alphath[-1,cnt]))
    
  }
  
  return(list('Vobj'=Vobj,'alphath'= alphath,'taugrid'=taugrid,'thdvar'=thdvar,'lambda'=lmbdcv))
}
# core function to estimate the threshold model by linear  desparstify Lasso. 

DSLasso <- function(y,XG, nfolds=10 ,intercept = FALSE ,ic='CV', standardize = FALSE)
{
  #init 
  n <- length(y)	
  # Storage
  alphals <-  Vobj<-  lmbdcv <-NULL
  
  
  # Construction of x tau  # Concatenation
  
  # cross validation for lambda selection
  if(ic=='CV')cgn	<- cv.glmnet(x = XG, y=y, nfolds=nfolds, standardize = standardize, intercept = intercept)
  if(ic=='BIC'){cgn	<- glmnet(x = XG, y=y,  intercept=intercept, standardize = standardize) 
  yhat    <- predict(cgn, newx=XG, type = "response")
  sigmah	<- colSums((matrix(y,nrow(yhat),ncol(yhat))-yhat)^2)/length(y)
  bic     <- log(sigmah) + cgn$df *log(n)/n 
  cgn$lambda.min 	<- cgn$lambda[which.min(bic)] 
  }	
  if(ic=='GIC'){cgn	<- glmnet(x = XG, y=y,  intercept=intercept, standardize = standardize) 
  yhat    <- predict(cgn, newx=XG, type = "response")
  sigmah	<- colSums((matrix(y,nrow(yhat),ncol(yhat))-yhat)^2)/n
  gic     <- log(sigmah) + cgn$df *log(log(n))/n * log(ncol(XG)) # Fan & Tang JRSS-B 2004
  cgn$lambda.min 	<- cgn$lambda[which.min(gic)]
  }	

  # Storing
  lmbdcv<- cgn$lambda.min
  if(ic=='CV') alphals  <- matrix(coef(cgn,s='lambda.min'))
  if(ic=='BIC') alphals <- c(cgn$a0[which.min(bic)],cgn$beta[,which.min(bic)])
  if(ic=='GIC') alphals <- c(cgn$a0[which.min(gic)],cgn$beta[,which.min(gic)])
  #if(intercept & BIC)  alphals <-c(cgn$a0[which.min(bic)], alphals)
  
  
  
  
  # Computing the value of the objective function:
  if(ic=='CV')yhat	 <- predict(cgn,newx = XG,s='lambda.min')
  RSS 	 <- sum((matrix(y,nrow(yhat),ncol(yhat))-yhat)^2)/(2*n)
  Vobj<- RSS + cgn$lambda.min*sum(abs(alphals[-1]))
  
  return(list('lmVobj'=Vobj,'lmalphals'=alphals,'lmlambda'=lmbdcv))
}
# Nodewise estimation of the covariance matrix
est_ndwcov <- function(Y,ic){
  # initialization
  p <- ncol(Y)
  n <- nrow(Y)
  C <- matrix(0,p,p)
  diag(C) <- 1
  tau <- NULL
  if(ic=='CV') ic<-'GIC'
  # Loop over the assets
  for(j in 1:p){
    # Estimate the Lasso
    jlas <- glmnet(x=Y[,-j],y=Y[,j],family = 'gaussian')
    # Get fit
    jfit <- predict(jlas, newx=Y[,-j], type="response")    
    # residuals
    jres <- matrix(Y[,j],n,length(jlas$lambda)) - jfit
    # std err
    jsig <- colSums(jres^2)/n
    # Computing information criterion
    if(ic=='WIC') jbic  <- log(jsig) + jlas$df * log(n)/n * log(log(p)) # BIC (Wang,2010)
    if(ic=='BIC') jbic  <- log(jsig) + jlas$df * log(n)/n  #BIC
    if(ic=='GIC') jbic  <- log(jsig) + jlas$df * log(p) * log(log(n))/n # Fan & Tang JRSS-B 2004
    if(ic=='AIC') jbic  <- log(jsig) + 2 * jlas$df # AIC 
    # Index of selected model 
    jind  <- which.min(jbic)
    # Get the parameters
    jpar <- jlas$beta[,jind]
    # Computing tau squared. Two formulas in text
    jtau <- sum(jres[,jind]^2)/n + jlas$lambda[jind]*sum(abs(jpar)) # using (12)
    
    # using the msgps package
    # Not used because it's very slow!
    #jlas <- msgps(Y[,-j],as.vector(Y[,j]))
    #jpar <- jlas$dfbic_result$coef[-1]
    #jtun <- jlas$dfbic_result$tuning
    #jres <- Y[,j]-predict(jlas,Y[,-j],tuning = jlas$dfbic_result$tuning)
    #jtau <- sum(jres^2) + jtun*sum(abs(jpar))
    
    # Storing the parameters
    C[j,-j] <- -jpar
    tau <- c(tau,jtau)
  }
  
  # Construct T-squared inverse
  T2inv <- diag(1/tau)
  
  # Construct Theta-hat
  Theta <- T2inv %*% C
  
  # sparsity
  sp <- sum(Theta==0)/(p^2)
  
  return(list(NULL,Theta,sp))
}
