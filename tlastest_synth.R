rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.
library(Matrix)
library(glmnet)
library(parallel)
library(reshape2)
library(ggplot2)
library(MASS)
library(lattice)
library(igraph)
library(flare)
setwd("/Users/hongqiangyan/Desktop/th lasso")
source('tlas_inference.R')


set.seed(17) #duh

# Number of cores for mclapply
ncores <-4
vtau <- c(0.5)
vcorrXQ <- c(0,0.9)
vobs <- c(500)
# Search grids
taugrid <- seq(0.15,0.85,0.05)
# Monte Carlo iterations
iter <-100
s_0<-10
numvars<-250
lpar1 <-list(c(0,0,rep(1,s_0),rep(0,numvars-s_0)))#,c(0,0,rep(1.5,s_0),rep(0,numvars-s_0)))
lpar2<-NULL
#lpar2 <- list(c(rep(1,s_0),rep(0,numvars-s_0)),c(rep(1.5,s_0),rep(0,numvars-s_0)))
#vpar=c(lpar1,lpar2)
h<- c(1/2,rep(0,s_0-1),1/2,rep(0,numvars-s_0-1),1/2,rep(0,s_0-1),1/2,rep(0,numvars-s_0-1))
intercept<-FALSE
thd_intercept <-FALSE
sig_eps=0.025
ic<- 'GIC'#######'CV' 'GIC''BIC' IF CV nfold=10&using GIC for nodewise reg
mctmp<- NULL
for(corrXQ in vcorrXQ){
  mca <- mc_ttlas(lpar1=lpar1,lpar2=lpar2,s_0=s_0,vobs=vobs,vtau=vtau,
                  taugrid=taugrid,iter=iter,
                  ncores = ncores,sig_eps=sig_eps,corrXQ=corrXQ,
                  intercept=intercept,thd_intercept=thd_intercept,ic=ic)
  
  if(is.null(mctmp)) mctmp <- array(NA,dim=c(dim(mca),length(vcorrXQ)),
                                    dimnames = c(dimnames(mca),list('corrXQ' = vcorrXQ)))
  
  mctmp[,,,,as.character(corrXQ)] <- mca
  
}

mmc <- melt(mctmp)
tmc <- acast(mmc,corrXQ+tau+nobs~stats)

file_name <- paste0(Sys.time(),' npar:',numvars,ic,s_0,'save.RData')
save(tmc,file = file_name)