# Attempt multivariate DLM by Pole method, for Paleo data
# SRC 2020-01-20 and Modified for Squeal Experiments
# (c) Stephen R. Carpenter

rm(list = ls())
graphics.off()

# ONLINE DYNAMIC LINEAR MODEL (DLM) ESTIMATION

DLM <- function(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat) {
  
  # Online algorithm for Dynamic linear regression 
  # Copyright 2016 by Stephen R. Carpenter
  
  # Description and definitions:
  
  # Observation equation is
  # Y_t = F_t'*theta_t + eta_t where
  # Y_t is the prediction
  # F_t is a vector of predictors at the beginning of the time step
  # theta_t is the parameter vector
  # eta_t is an individual observation error
  
  # System equation is:
  # theta_t = theta_t-1 + omega_t
  # where theta is defined above and omega_t is an individual process error
  
  # Inputs to the function are:
  # delta, the discount factor
  # n.gamma, the initial number of observations (usually 1)
  # d.gamma, the initial shape parameter for prediction errors
  #  (prior estimate of prediction variance = d.gamma / n.gamma)
  # mvec, the initial guess of regression coefficients
  # Cpar, the initial guess of the covariance matrix of regression coefficients
  # Yvec, the vector of the observed response variate
  # Fmat, the matrix of predictors
  
  # Outputs are:
  # predix, the one-step-ahead predictions of the response variate
  # varpredix, the prediction variance at start of time step before error is measured
  # pars, the updated parameter estimates using the most recent prediction error
  # parvar, the variances of the parameters
  # Svec, the update (after error is measured within a time step) of varpredix
  
  # Updating follows the equations on p. 176-179 of Carpenter 2003,
  # Regime Shifts in Lake Ecosystems: Pattern and Variation
  
  # Determine constants
  npar <- length(mvec)
  Nobs <- length(Yvec)
  S0 <- d.gamma/n.gamma
  
  # Set up vectors to hold results
  predix <- rep(0,Nobs)
  varpredix <- rep(0,Nobs)
  Svec = rep(0,Nobs)
  pars <- matrix(0,nrow=Nobs,ncol=npar)
  parvar = matrix(0,nrow=Nobs,ncol=npar)
  
  for(i in 1:Nobs)  {  #Start DLM loop
    # Generate predictions
    Fvec <- Fmat[i,] # vector of predictors
    predix[i] <- sum(Fvec*mvec)
    # Compute error and update estimates
    error <- Yvec[i]-predix[i]
    Rmat <- Cpar/delta
    varpredix[i] <- (t(Fvec) %*% Rmat %*% Fvec) + S0
    n.gamma <- (delta*n.gamma)+1
    d.gamma <- (delta*d.gamma)+(S0*error*error/varpredix[i])
    S1 <- d.gamma/n.gamma
    Svec[i] = S1  # save updated variance
    Avec <- (Rmat %*% Fvec)/varpredix[i]
    mvec <- mvec + (Avec*error)
    pars[i,] <- mvec
    Cpar <- (S1/S0)*(Rmat - (Avec %*% t(Avec))*varpredix[i])
    # Disallow negative variances on the diagonal
    for(idiag in 1:npar) {
      Cpar[idiag,idiag] <- max(0,Cpar[idiag,idiag])
    }
    parvar[i,] = diag(Cpar)
    S0 <- S1 # roll over S
  } # End DLM loop
  
  DLM.out <- list(predix,varpredix,pars,parvar,Svec)
  return(DLM.out)
} # END DLM FUNCTION

# Main program -----------------------------------------------------------------------------------

# Load the data
# Save block:
# variates: DoY Year darkBGA dark.lBGA  dDOsat  dpH 
#   Chl_Manual Zmix_daily      DOC
#save(Rdaily,Tdaily,file='Daily_PeterTuesday2015_allvars.Rdata')
load(file='Daily_PeterTuesday2015_allvars.Rdata')
dat0 = Tdaily
dat0$lchl = log10(dat0$Chl_Manual)
title = c('B. Tuesday Lake')
# Make a matrix with the time series as columns
#Xmat.0 = dat0[,4:6]  # keep log10 phycocyanin, Chl, delta DOsat, delta pH
Xmat.0 = subset(dat0,select=c(dark.lBGA,dDOsat,dpH,lchl))
Xmat.1 = as.matrix(Xmat.0)
dimX = dim(Xmat.1)
nX = dimX[1] # number of observations
nV = dimX[2] # number of variates
# make a time variable
tvec = dat0[,1]

windows(height=10,width=5)
par(mfrow=c(4,1),mar=c(2.5, 4.5, 3, 1) + 0.1, cex.axis=1.5,cex.lab=1.8,
    cex.main=1.6)
plot(tvec,Xmat.1[,1],type='l',lwd=2,col='forestgreen',xlab='',
     ylab='log10 Phycocyanin',main=title)
par(mar=c(2.5, 4.5, 1, 1) + 0.1)
plot(tvec,Xmat.1[,4],type='l',lwd=2,col='red',xlab='Day of Year',
     ylab='log10 Chlorophyll')
plot(tvec,Xmat.1[,2],type='l',lwd=2,col='blue',xlab=' ',
     ylab='Daily Range, D.O. Sat.')
par(mar=c(4, 4.5, 1, 1) + 0.1)
plot(tvec,Xmat.1[,3],type='l',lwd=2,col='red',xlab='Day of Year',
     ylab='Daily Range, pH')

# optional transformations
######################
#Xmat.1 <- log(1+Xmat.1)  # log transform (optional)
######################
unit = rep(1,nX)
cmean = colMeans(Xmat.1)
Xmat.2 = Xmat.1 - (unit%*%t(cmean)) # center (optional)
######################
csd = apply(Xmat.1,2,sd,na.rm=T)
invcsd = 1/csd
Xmat.3 = Xmat.1*(unit%*%t(invcsd)) # make z scores
#####################

### What matrix will be analyzed?
#Xmat = Xmat.3
Xmat = as.matrix(cbind(unit,Xmat.3)) # include intercept

# set number of parameters (remember to add 1 if there is an intercept)
npar = nV+1

# Quick parameter estimate using LS 
print('Multivariate least squares regression to check the data',quote=F)
Xmat1 = Xmat[2:nX,]
Xmat0 = Xmat[1:(nX-1),]
Xinv = solve(t(Xmat0)%*%Xmat0)
bmat = Xinv%*%t(Xmat0)%*%Xmat1
print('',quote=F)
print('b matrix',quote=F)
print(bmat)
print('eigenvalues of b matrix',quote=F)
beig = eigen(bmat,only.values=T)
print(beig$values)
print('prediction error variances',quote=F)
yhat = Xmat0%*%bmat
err = Xmat1 - yhat
verr = apply(err,2,var)
print(verr)
print('parameter covariance matrices',quote=F)
print('each of the regressions has the same t(x)%*%x but unique error variance',quote=F)
parcov = array(data=0,dim=c(npar,npar,npar))
for(i in 1:npar) {
  parcov[i,,] = Xinv*verr[i]
  #print(parcov[i,,])
}
print('printing suppressed if there are lots of covariance matrices',quote=F)

# Set up for multivariate DLM --------------------------------------------------------------

# Additional inputs
delta = rep(0.95,npar) # vector of 1 value per parameter, e.g. c(0.9,0.9,0.9)
n.gamma = rep(1,npar)
d.gamma = verr

# Output arrays to hold results
par.dlm = array(data=0,dim=c((nX-1),npar,npar))

# Run the npar DLMs and save results
yhatmat = matrix(0,nr=(nX-1),nc=npar)
for(idlm in 1:npar) {
  MAR.est = DLM(delta[idlm],n.gamma[idlm],d.gamma[idlm],bmat[,idlm],parcov[idlm,,],Xmat1[,idlm],Xmat0)
  par.dlm[,,idlm] = MAR.est[[3]]
  # Plot predix and obs
  yhat = MAR.est[[1]]
  yhatmat[,idlm] = yhat
  if(idlm > 1) {
    print(c('Corr of prediction and observation for column ',idlm),quote=F)
    print(cor(yhat,Xmat1[,idlm]),quote=F)
  }
  windows()
  par(mfrow=c(1,1),mar=c(4, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
  plot(tvec[2:nX],yhat,type='l',lwd=3,col='deepskyblue',xlab='DOY',ylab='Y and Yhat',
       main=bquote('Predix & Obs for column'~ .(idlm)) )
  points(tvec[2:nX],Xmat1[,idlm],type='p',pch=19,col='darkblue')
}

# Predictions & observations of log10 BGA
windows()
par(mfrow=c(1,1),mar=c(4, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(tvec[2:nX],yhatmat[,1],type='l',lwd=3,col='deepskyblue',xlab='DOY',
     ylab='log10 Phycocyanin',main='Tuesday Lake 2015, r = 0.96')
points(tvec[2:nX],Xmat1[,1],type='p',pch=19,col='darkblue')
legend('topleft',legend=c('One-step prediction','Observed'),lwd=c(3,NA),
       cex=1.4,pch=c(NA,19),col=c('deepskyblue','darkblue'),bty='n') 

# Predictions & observations of log10 Chl
windows()
par(mfrow=c(1,1),mar=c(4, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(tvec[2:nX],yhatmat[,4],type='l',lwd=3,col='deepskyblue',xlab='DOY',
     ylab='log10 Chlorophyll',main='Tuesday Lake 2015, r = 0.95')
points(tvec[2:nX],Xmat1[,4],type='p',pch=19,col='darkblue')
legend('topleft',legend=c('One-step prediction','Observed'),lwd=c(3,NA),
       cex=1.4,pch=c(NA,19),col=c('deepskyblue','darkblue'),bty='n') 

windows()
par(mfrow=c(1,1),mar=c(4, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(yhatmat[,2],yhatmat[,3],type='p',pch=19,cex=0.8,col='darkblue',
     xlab='delta DO, r = 0.34',ylab='delta pH, r = 0.61',
     main='One-step predictions, Tuesday Lake 2015')
arrows(x0=yhatmat[1:(nX-2),2],x1=yhatmat[2:(nX-1),2],
       y0=yhatmat[1:(nX-2),3],y1=yhatmat[2:(nX-1),3],lwd=1,col='blue')

# Construct the b matrices at each time and compute eigenvalues
eigvals = rep(0,(nX-1)) # vector to hold maximum eigenvalues
weights1 = matrix(0,nr=(nX-1),nc=(npar-1)) # subtract 1 column due to removing intercept
for(it in 1:(nX-1)) {
  b.it = par.dlm[it,,]
  b.noint = b.it[2:npar,2:npar] # remove intercept
  lam = eigen(b.noint,only.values=F)
  lam.max = max(Mod(lam$values))
  imax = which.max(Mod(lam$values))
  eigvals[it] = sign(Re(lam$values[imax]))*Re(lam.max)
  weights1[it,] = lam$vectors[,1]
}

windows()
par(mfrow=c(1,1),mar=c(5, 5, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(tvec[2:nX],eigvals,type='b',lwd=2,pch=19,col='red',
     xlim=(range(tvec[2:nX])),xlab='Day of Year 2015',ylab='Eigenvalues',
     main=title)
abline(h=1,lty=2,lwd=2)
abline(v=152,lty=2,lwd=2,col='blue')
abline(v=240,lty=2,lwd=2,col='blue')
