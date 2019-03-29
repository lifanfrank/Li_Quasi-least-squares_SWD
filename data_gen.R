#############################################################
# Generate correlated Gaussian outcomes
# for cohort stepped wedge CRTs

# n: Number of clusters
# m: Cohort size
# t: Number of periods
# delta: Effect size
# s2: Dispersion parameter or total variance (assume to be 1
#     and suppressed in the input)
# beta: Vector of period effects
# alpha: Vector of correlations 
#        tau=within period correlation
#        rho=decay parameter
#############################################################

contGEN<-function(n,m,t,delta,beta,alpha){
  require(mvtnorm)
  
  ########################################################
  # Create proportional decay correlation matrix.
  ########################################################
  
  pdcorr<-function(alpha){
    tau<-alpha[1]
    rho<-alpha[2]
    Fmat<-rho^abs(outer(0:(t-1),0:(t-1),"-"))
    Gmat<-(1-tau)*diag(1,m)+tau*matrix(1,m,m)
    return(kronecker(Gmat,Fmat))
  }
  
  ########################################################
  # returns variance matrix of Gaussian variables with dispersion
  # parameter s2 and corr matrix r[,].
  ########################################################
  
  cor2var<-function(r,s2){
    return(s2*r)
  }
  
  # Create treatment sequences
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  g<-n/(t-1) # number of clusters per step
  
  # Simulate correlated Gaussian outcomes
  s2<-1
  y<-NULL
  r<-pdcorr(alpha)
  v<-cor2var(r,s2)   # v <- cov matrix
  for(i in 1:(t-1)){
    u_c<-c(beta+delta*trtSeq[i,])
    u<-rep(u_c,m)
    y<-cbind(y,t(rmvnorm(g,u,v))) # simulate data matrix
  }
  
  # Return simulated data matrix
  return(y)
}