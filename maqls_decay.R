#####################################################################################
# QLS analysis of cohort stepped wedge cluster randomized trials 
# with continuous outcomes and a decaying correlation structure

# This program implements the regular QLS and the "matrix-adjusted" QLS 
# for estimating the marginal mean parameters along with the proportional
# decay correlation structure

# Fan Li; March 2019

# The following variance estimators for the marginal mean model are output from the program:
# (1) MB: Model-based estimator
# (2) BC0: Uncorrected sandwich estimator that extends Liang and Zeger, 1986
# (3) BC1: Bias-corrected sandwich estimator that extends Kauermann and Carroll, 2001
# (4) BC2: Bias-corrected sandwich estimator that extends Mancl and DeRouen, 2001
# (5) BC3: Bias-corrected sandwich estimator that extends Fay and Graubard, 2001

# Note that all inputs are required.
# ID's should be integers from 1 to K. Data should be sorted by ID before 
# calling contMAEE.  Z matrix should be all pairs j < k for each cluster i, ie
# (1,2) (1,3).... (1,n)... (2,3).... (n-1, n) for each cluster.

# INPUT
# y: The continuous outcome variable
# X: Marginal mean covariates (design matrix including intercept)
# id: Cluster identifier
# n: Vector of cluster sample sizes (across periods)
# m: Vector of cohort sizes
# t: Number of periods
# maxiter: Maximum number of iterations for Fisher Scoring
# epsilon: Tolerance for convergence
# alpadj: Bias adjustment to correlation estimating equations: 
#         1 - no adjustment (Shults and Morrow, 2002)
#         3 - matrix-adjustment described in the article  (Li, 2019+)    
#####################################################################################

contQLS=function(y, X, id, n, m, t, maxiter, epsilon, alpadj){
  require(MASS)
  require(rootSolve)
  
  #####################################################################################
  # MODULE: BEGINEND
  # Creates two vectors that have the start and end points for each cluster
  
  # INPUT
  # n: Vector of cluster sample sizes
  
  # OUTPUT
  # first: Vector with starting row for cluster i
  # last: Vector with ending row for cluster i
  #####################################################################################
  
  BEGINEND=function(n){
    last=cumsum(n)
    first=last-n+1
    return(cbind(first,last))
  }
  
  #####################################################################################
  # Module: IS_POS_DEF
  # A = symmetric matrix
  # returns 1 if A is positive definite
  # 0 otherwise
  #####################################################################################
  
  is_pos_def=function(A){
    return(min(eigen(A)$values)>1e-13)
  }
  
  #####################################################################################
  # MODULE: CREATEA
  # Creates residual for the marginal mean estimating equation, (Y - mu)
  
  # INPUT
  # mu: Vector of n_i marginal means
  # y: Outcome vector for ith cluster
  
  # OUTPUT
  # residuals for the marginal mean estimating equation
  #####################################################################################
  
  CREATEA=function(mu,y){
    return(y-mu)
  }
  
  #####################################################################################
  # MODULE: CREATEB
  # Creates covariance matrix for the marginal mean estimating equation, var(Y)
  
  # INPUT
  # s2: Dispersion
  # alpha: (alpha_0, alpha_1), parameterizing the propritional decay corr structure
  # m: Cohort size (scalar) for cluster i
  # t: Number of periods
  
  # OUTPUT
  # covariance matrix for the marginal mean estimating equation
  #####################################################################################
  
  CREATEB=function(s2,alpha,m,t){
    alpha0=alpha[1]
    alpha1=alpha[2]
    Fmat=alpha1^abs(outer(0:(t-1),0:(t-1),"-"))
    Gmat=(1-alpha0)*diag(1,m)+alpha0*matrix(1,m,m)
    R=kronecker(Gmat,Fmat)
    Ahalf=diag(rep(sqrt(s2),m*t))
    B=Ahalf%*%R%*%Ahalf
    return(B)
  }
  
  # TBC
  #####################################################################################
  # MODULE: INVBMAT
  # Invert covariance matrix for the marginal mean estimating equation, var(Y)
  
  # INPUT
  # s2: Dispersion
  # alpha: (alpha_0, alpha_1), parameterizing the propritional decay corr structure
  # m: Cohort size (scalar) for cluster i
  # t: Number of periods
  
  # OUTPUT
  # Inverse of the covariance matrix for the marginal mean estimating equation
  #####################################################################################
  
  INVBMAT=function(s2,alpha,m,t){
    alpha0=alpha[1]
    alpha1=alpha[2]
    INVFmat=solve(alpha1^abs(outer(0:(t-1),0:(t-1),"-")))
    INVGmat=(1/(1-alpha0))*diag(1,m)-(alpha0/((1-alpha0)*(1+(m-1)*alpha0)))*matrix(1,m,m)
    INVR=kronecker(INVGmat,INVFmat)
    INVAhalf=diag(rep(1/sqrt(s2),m*t))
    INVB=INVAhalf%*%INVR%*%INVAhalf
    return(INVB)
  }
  
  #####################################################################################
  # MODULE: SCORE
  # Generates the score matrix for each cluster and approximate information
  # to be used to estimate parameters and generate standard errors
  
  # INPUT
  # beta: Vector of marginal mean parameters
  # alpha: Vector of marginal correlation parameters
  # y: Vector of continuous outcomes
  # X: Marginal mean covariates
  # n: Vector of cluster sample sizes (across periods)
  # m: Vector of cohort sizes
  # t: Number of periods
  # p: Number of marginal mean parameters
  # s2: Dispersion
  
  # OUTPUT
  # U: Score vector
  # UUtran: Sum of U_i*U_i` across all clusters
  # Ustar: Approximate information matrix
  #####################################################################################
  
  SCORE=function(beta, alpha, y, X, n, m, t, p, s2){
    U=rep(0,p)
    UUtran=Ustar=matrix(0,p,p)
    locx=BEGINEND(n)
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      
      U_c=rep(0,p)
      Ustar_c=matrix(0,p,p)
      mu_c=c(X_c%*%beta)
      
      C=X_c
      A=CREATEA(mu_c,y_c)
      B=CREATEB(s2,alpha,m[i],t)
      INVB=INVBMAT(s2,alpha,m[i],t)
      
      U_c=t(C)%*%INVB%*%A
      UUtran_c=tcrossprod(U_c)
      Ustar_c=t(C)%*%INVB%*%C
      
      U=U+U_c
      UUtran=UUtran+UUtran_c
      Ustar=Ustar+Ustar_c
    }
    return(list(U=U,UUtran=UUtran,Ustar=Ustar))
  }
  
  #####################################################################################
  # MODULE: STAGEONE
  # Obtain Stage-one estimates for the correlation parameters
  
  # INPUT
  # Ustarold: Initial values for information matrix
  # alphaold: Initial values for the correlation parameters
  # y: Vector of continuous outcomes
  # X: Marginal mean covariates
  # n: Vector of cluster sample sizes (across periods)
  # m: Vector of cohort sizes
  # t: Number of periods
  # s2: Dispersion
  # alpadj: Bias adjustment to correlation estimating equations: 
  #         1 - no adjustment (Shults and Morrow, 2002)
  #         3 - matrix-adjustment described in the article      
  
  # OUTPUT
  # res: Stage-one estimates for the correlation parameters
  #####################################################################################
  
  STAGEONE=function(Ustarold, alphaold, beta, y, X, n, m, t, s2, alpadj){
    naiveold=ginv(Ustarold)
    locx=BEGINEND(n)
    
    # simplified estimating equation for alpha
    f=function(alpha){
      alpha0=alpha[1]
      alpha1=alpha[2]
      f0=f1=0
      
      # create inverse of AR1 correlation matrix
      IM=diag(1,t)
      C2=diag(c(0,rep(1,t-2),0),t)
      C1=diag(0,t)
      C1[cbind(2:t,1:(t-1))]=1
      C1[cbind(1:(t-1),(2:t))]=1
      INVF=IM+alpha1^2*C2-alpha1*C1 # omit proportionality constant
      
      for(i in 1:length(n)){
        X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
        y_c=y[locx[i,1]:locx[i,2]]
        mu_c=c(X_c%*%beta)
        C=X_c
        B=CREATEB(s2,alpha,m[i],t)
        INVB=INVBMAT(s2,alpha,m[i],t)
        
        # no adjustment
        if(alpadj==1){
          SQVARFUN=sqrt(rep(s2,n[i]))
          INVSQVAR=1/SQVARFUN
          RX=(y_c-mu_c)*INVSQVAR
          RM=tcrossprod(RX)
        }
        
        # matrix-adjusted correlation estimates
        if(alpadj==3){
          SQVARFUN=sqrt(rep(s2,n[i]))
          INVSQVAR=1/SQVARFUN
          CT=t(C)
          omega=C%*%naiveold%*%CT
          vminomega=B-omega
          Ci=B%*%ginv(vminomega)
          RX=(y_c-mu_c)*INVSQVAR
          Gi=tcrossprod(RX)
          RM=Ci%*%Gi
        }
        
        # alpha0 equation, omit proportionality constant
        Dqi1=alpha0^2*(m[i]-2)*(m[i]-1)+2*alpha0*(m[i]-1) 
        Dqi2=-(1+alpha0^2*(m[i]-1))                       
        DINVG=matrix(Dqi2,m[i],m[i])
        diag(DINVG)=Dqi1
        f0=f0+sum(diag(kronecker(DINVG,INVF)%*%RM))
        
        # alpha1 equation, omit proportionality constant
        qi1=1+(n[i]-2)*alpha0
        qi2=-alpha0
        INVG=matrix(qi2,m[i],m[i])
        diag(INVG)=qi1
        DINVF=2*alpha1*IM+2*alpha1*C2-(1+alpha1^2)*C1   
        f1=f1+sum(diag(kronecker(INVG,DINVF)%*%RM))
      }
      # return value
      return(c(f0,f1))
    }
    # optimization
    res<-multiroot(f,start=alphaold)$root
    return(res)
  }
  
  #####################################################################################
  # MODULE: STAGETWO
  # Obtain Stage-two estimates for the correlation parameters
  
  # INPUT
  # alpha: Stage-one estimates for the correlation parameters 
  
  # OUTPUT
  # tau: Estimate for within-period correlation, tau
  # rho: Estimate for decay parameter, rho
  #####################################################################################
  
  STAGETWO=function(alpha){
    alpha0=alpha[1]
    alpha1=alpha[2]
    tau_num=sum( m*(m-1)*alpha0*(2+(m-2)*alpha0) / (1+(m-1)*alpha0)^2 )
    tau_den=sum( m*(m-1)*(1+(m-1)*alpha0^2) / (1+(m-1)*alpha0)^2 )
    tau=tau_num/tau_den
    rho=2*alpha1/(1+alpha1^2)
    return(c(tau,rho))
  }
  
  #####################################################################################
  # MODULE: PHI
  # Estimate the dispersion parameter
  
  # INPUT
  # Ustar: Estimated information matrix
  # beta: Vector of marginal mean parameters
  # alpha: Vector of marginal correlation parameters
  # y: Vector of continuous outcomes
  # X: Marginal mean covariates
  # n: Vector of cluster sample sizes (across periods)
  # m: Vector of cohort sizes
  # t: Number of periods
  # p: Number of marginal mean parameters
  # s2: Dispersion estimate from the last update
  # alpadj: Bias adjustment to correlation estimating equations: 
  #         1 - no adjustment (Shults and Morrow, 2002)
  #         3 - matrix-adjustment described in the article      
  
  # OUTPUT
  # s2: Updated dispersion estimate
  # NPSDADJFLAG: Check whether correction factor is singular
  #             (1: singular, stop the program; 0: pass the check, progam resumes)
  #####################################################################################
  
  PHI=function(Ustar, beta, alpha, y, X, n, m, t, p, s2, alpadj){
    
    # initiate values
    RSSC=0
    NPSDADJFLAG=0
    naive=ginv(Ustar)
    locx=BEGINEND(n)
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      mu_c=c(X_c%*%beta)
      
      C=X_c
      B=CREATEB(s2,alpha,m[i],t)
      INVB=INVBMAT(s2,alpha,m[i],t)
      CtinvB=t(C)%*%INVB
      Hi1=C%*%naive%*%CtinvB
      
      if(alpadj==1){
        SQVARFUN=sqrt(rep(s2,n[i]))
        INVSQVAR=1/SQVARFUN
        RX=(y_c-mu_c)*INVSQVAR
        RM=tcrossprod(RX)
      }
      
      # Insert checks here to avoid singular correction factor
      if(alpadj==3){
        SQVARFUN=sqrt(rep(s2,n[i]))
        INVSQVAR=1/SQVARFUN
        CT=t(C)
        omega=C%*%naive%*%CT
        vminomega=B-omega
        psd_vmin=is_pos_def(vminomega)
        mineig=min(eigen(vminomega)$values)
        if(psd_vmin==1){
          # if(mineig<1e-12){cat("vminomega with psd_min =",psd_min,"mininum eigenvalue is",mineig,"\n")}
          Ci=B%*%ginv(vminomega)
          RX=(y_c-mu_c)*INVSQVAR
          Gi=tcrossprod(RX)
          RM=Ci%*%Gi
        } else{
          NPSDADJFLAG=1
          stop("vminomega is not positive definite")
        }
      }
      RSSC=RSSC+sum(s2*diag(RM))
    }
    # moment-based estimator
    s2=RSSC/(sum(n)-p)
    return(list(s2=s2,NPSDADJFLAG=NPSDADJFLAG))
  }
  
  #####################################################################################
  # MODULE: INITBETA
  # Generates initial values for beta. Linear regression using least squares
  
  # INPUT
  # y: Vector of continuous outcomes
  # X: Marginal mean covariates
  # n: Vector of cluster sample sizes (across periods)
  
  # OUTPUT
  # beta: Vector of marginal mean parameters
  # Ustar: approximate score vector
  # s2: Dispersion
  #####################################################################################
 
  INITBETA=function(y,X,n){
    beta=solve(t(X)%*%X,t(X)%*%y)
    u=c(X%*%beta)
    s2=sum((y-u)^2)/(sum(n)-ncol(X))
    Ustar=t(X)%*%X/s2
    return(list(beta=c(beta),Ustar=Ustar,s2=s2))
  }
  
  #####################################################################################
  # MODULE: INVBIG
  # compute (A - mm`)^{-1}c without performing the inverse directly
  # Computation details found in Preisser, Qaqish and Perin. (2008)
  # A note on deletion diagnostics for estimating equations. Biometrika
  
  # INPUT
  # ainvc: inverse of matrix A times vector c
  # ainvm: inverse of matrix A times matrix (with low no. of columns) M
  # M: matrix of eigen column vectors m1,m2, ..
  # c: vector 
  # start: of do loop
  # end:   of do loop, rank of X
  
  # OUTPUT
  # ainvc: inverse of matrix A times vector c
  #####################################################################################
  
  INVBIG=function(ainvc,ainvm,m,c,start,end){
    for(i in start:end){
      b=ainvm[,i]
      bt=t(b)
      btm=bt%*%m
      btmi=btm[,i]
      gam=1-btmi
      bg=b/gam
      ainvc=ainvc+bg%*%(bt%*%c)
      if(i<end){
        ainvm=ainvm+bg%*%btm
      }
    }
    return(ainvc)
  }
  
  #####################################################################################
  # MODULE: MAKEVAR
  # Creates covariance matrix of the marginal mean parameters
  
  # INPUT
  # beta: Vector of marginal mean parameters
  # alpha: Vector of marginal correlation parameters
  # y: Vector of continuous outcomes
  # X: Marginal mean covariates
  # n: Vector of cluster sample sizes (across periods)
  # m: Vector of cohort sizes
  # t: Number of periods
  # p: Number of marginal mean parameters
  # s2: Dispersion estimate from the last update
  
  # OUTPUT
  # robust: Robust covariance matrix for beta and alpha
  # naive: Naive (Model-Based) covariance matrix for beta
  # varMD: bias-corrected variance by Mancl and Derouen (2001)
  # varKC: bias-corrected variance by Kauermann and Carroll (2001)
  # varFG: bias-corrected variance by Fay and Graubard (2001)
  # ROBFLAG: Check singularity of the estimated covariate
  #         (1: singular; 0: pass the check)
  #####################################################################################
  
  MAKEVAR=function(beta, alpha, y, X, n, m, t, p, s2){
    
    # Run the score module
    SCORE_RES=SCORE(beta, alpha, y, X, n, m, t, p, s2)
    U=SCORE_RES$U
    UUtran=SCORE_RES$UUtran
    Ustar=SCORE_RES$Ustar
    ROBFLAG=0
    
    # Model-based estimator
    naive=ginv(Ustar)
    
    # new commands to compute INV(I - H1)
    eigenRES1=eigen(naive)
    evals1=eigenRES1$values
    evecs1=eigenRES1$vectors
    sqrevals1=sqrt(evals1)
    sqe1=evecs1%*%diag(sqrevals1)
    
    # Bias-corrected variance
    Ustar_c_array=UUtran_c_array=array(0,c(p,p,length(n)))
    UUtran=UUbc=UUbc2=UUbc3=Ustar=matrix(0,p,p)
    
    locx=BEGINEND(n)
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      mu_c=c(X_c%*%beta)
      
      U_i=U_c=rep(0,p)
      Ustar_c=matrix(0,p,p)
      
      # commands for beta
      C=X_c
      A=CREATEA(mu_c,y_c)
      B=CREATEB(s2,alpha,m[i],t)
      INVB=INVBMAT(s2,alpha,m[i],t)
      U_i=t(C)%*%INVB%*%A
      
      # commands for generalized inverse - beta
      ai1=INVB
      mm1=C%*%sqe1
      ai1A=ai1%*%A
      ai1m1=ai1%*%mm1
      ai1A=INVBIG(ai1A,ai1m1,mm1,A,1,p)
      U_c=t(C)%*%ai1A
      
      Ustar_c=t(C)%*%INVB%*%C
      Ustar=Ustar+Ustar_c
      UUtran_c=tcrossprod(U_i)
      UUtran=UUtran+UUtran_c
      UUbc_c=tcrossprod(U_c)
      UUbc=UUbc+UUbc_c
      UUbc_ic=tcrossprod(U_c,U_i)
      UUbc2=UUbc2+UUbc_ic
      
      Ustar_c_array[,,i]=Ustar_c
      UUtran_c_array[,,i]=UUtran_c
    }
    
    # Naive or Model-based estimator
    naive=ginv(Ustar)
    
    # BC0 or usual Sandwich estimator     
    robust=naive%*%UUtran%*%t(naive)
    
    # BC1 or Variance estimator due to Kauermann and Carroll (2001);
    varKC=naive%*%(UUbc2+t(UUbc2))%*%t(naive)/2
    
    # BC2 or Variance estimator due to Mancl and DeRouen (2001);
    varMD=naive%*%UUbc%*%t(naive)
    
    # calculating adjustment factor for BC3
    for(i in 1:length(n)){      
      Hi=diag(1/sqrt(1-pmin(0.75,c(diag(Ustar_c_array[,,i]%*%naive)))))
      UUbc3=UUbc3+Hi%*%UUtran_c_array[,,i]%*%Hi
    }
    
    # BC3 or Variance estimator due to Fay and Graubard (2001);
    varFG=naive%*%UUbc3%*%t(naive)
    
    if(min(diag(robust))<=0){ROBFLAG = 1}
    if(min(diag(varMD))<=0){ROBFLAG = 1}
    if(min(diag(varKC))<=0){ROBFLAG = 1}
    if(min(diag(varFG))<=0){ROBFLAG = 1}
    
    return(list(robust=robust,naive=naive,varMD=varMD,varKC=varKC,varFG=varFG,ROBFLAG=ROBFLAG))
  }
  
  #####################################################################################
  # MODULE: FITQLS
  # Performs the fitting algorithm of QLS/MAQLS
  
  # INPUT
  # y: Vector of continuous outcomes
  # X: Marginal mean covariates
  # n: Vector of cluster sample sizes (across periods)
  # m: Vector of cohort sizes
  # t: Number of periods
  # maxiter: Max number of iterations
  # epsilon: Tolerence for convergence
  # SINGFLAG: THE ALGORITHM terminated due to singular MB covariance matrix
  # ROBFLAG: THE ALGORITHM terminated due to singular robust variance matrix
  # ALPFLAG: THE ALGORITHM terminated due to infeasible correlation estimates
  # NPSDADJFLAG: Checks if the alpha adjustment factor matrix is positive definite, terminates if not
  
  # OUTPUT
  # beta: A p x 1 vector of marginal mean parameter estimates  
  # alpha: marginal correlation parameter estimates, tau and rho
  # phi: dispersion parameter
  # robust: Robust covariance matrix for beta
  # naive: Model-Based covariance matrix for beta
  # varMD: Bias-corrected variance due to Mancl and DeRouen (2001)
  # varKC: Bias-corrected variance due to Kauermann and Carroll (2001)
  # varFG: Bias-corrected variance due to Fay and Graubard (2001)
  # niter: Number of iterations required for convergence
  # converge: Did the algorithm converge (0 = no, 1 = yes)
  #####################################################################################
  FITQLS=function(y, X, n, m, t, maxiter, epsilon, 
                  SINGFLAG, ROBFLAG, ALPFLAG, NPSDADJFLAG){
    p=ncol(X)
    delta=rep(2*epsilon,p)
    max_modi=19
    
    # estimate beta
    alphaold=alphastart=rep(0.01,2)
    INITRES=INITBETA(y,X,n)
    beta=INITRES$beta
    Ustar=INITRES$Ustar
    s2=INITRES$s2
    
    # estimate alpha
    n_modi=0
    repeat{
      alpha=STAGEONE(Ustar, alphaold, beta, y, X, n, m, t, s2, alpadj=1)
      if((alpha[1]>-1/(max(m)-1)) & (alpha[1]<1) & (abs(alpha[2])<1)){
        ALPFLAG=0
      } else{
        ALPFLAG=1
      }
      if(ALPFLAG==1){
        alphaold=alphaold+c(0.05,0.05)
        n_modi=n_modi+1
      }
      if((n_modi>max_modi) | (ALPFLAG==0)){break}
    }
    alphanew=STAGETWO(alpha)
    if(ALPFLAG==1){
      stop("ALPFLAG==1 program terminated due to infeasible correlation estimates")
    }
    
    # Iterated updates
    niter=1
    converge=0
    while((niter<=maxiter) & (converge==0)){
      # beta estimating equations
      SCORE_RES=SCORE(beta, alphanew, y, X, n, m, t, p, s2)
      U=SCORE_RES$U
      UUtran=SCORE_RES$UUtran
      Ustar=SCORE_RES$Ustar
      psdustar=is_pos_def(Ustar)
      mineig=min(eigen(Ustar)$values)
      if(psdustar==TRUE){
        delta=solve(Ustar,U)
        beta=beta+delta
        converge=(max(abs(delta))<=epsilon)
      } else{
        SINGFLAG=1
      }
      
      # alpha estimating equations
      n_modi=0
      alphaold=alphastart
      repeat{
        alpha=STAGEONE(Ustar, alphaold, beta, y, X, n, m, t, s2, alpadj)
        if((alpha[1]>-1/(max(m)-1)) & (alpha[1]<1) & (abs(alpha[2])<1)){
          ALPFLAG=0
        } else{
          ALPFLAG=1
        }
        if(ALPFLAG==1){
          alphaold=alphaold+c(0.05,0.05)
          n_modi=n_modi+1
        }
        if((n_modi>max_modi) | (ALPFLAG==0)){break}
      }
      alphanew=STAGETWO(alpha)
      
      # dispersion parameter
      PHI_RES=PHI(Ustar, beta, alphanew, y, X, n, m, t, p, s2, alpadj)
      s2=PHI_RES$s2
      NPSDADJFLAG=PHI_RES$NPSDADJFLAG
      
      # counter
      niter<-niter+1
    }
    
    if(converge==1){
      MAKEVAR_RES=MAKEVAR(beta, alphanew, y, X, n, m, t, p, s2)
      robust=MAKEVAR_RES$robust
      naive=MAKEVAR_RES$naive
      varMD=MAKEVAR_RES$varMD
      varKC=MAKEVAR_RES$varKC
      varFG=MAKEVAR_RES$varFG
      ROBFLAG=MAKEVAR_RES$ROBFLAG
    }
    return(list(beta=beta,alpha=alphanew,phi=s2,robust=robust,naive=naive,varMD=varMD,
                varKC=varKC,varFG=varFG,niter=niter,converge=converge,SINGFLAG=SINGFLAG,
                ROBFLAG=ROBFLAG,ALPFLAG=ALPFLAG,NPSDADJFLAG=NPSDADJFLAG))
  }
  
  #####################################################################################
  # MODULE: RESULTS
  # Creates printed output to screen of parameters and other information
  
  # INPUT
  # beta: Vector of marginal mean parameters
  # alpha: Vector of marginal correlation Parameters 
  # phi: Dispersion
  # robust: Robust covariance matrix for beta and alpha
  # naive: Model-Based covariance matrix for marginal mean
  # varMD: Bias-corrected variance due to Mancl and DeRouen (2001)
  # varKC: Bias-corrected variance due to Kauermann and Carroll (2001)
  # varFG: Bias-corrected variance due to Fay and Graubard (2001)
  # niter: Number of iterations until convergence
  # n: Vector of cluster sample sizes (across periods)
  
  # OUTPUT
  # To Screen
  #####################################################################################
  
  RESULTS=function(beta, alpha, phi, robust, naive, varMD, varKC, varFG, niter, n){
    p=length(beta)
    K=length(n)
    
    # beta variances
    bSE=sqrt(diag(naive))
    bSEBC0=sqrt(diag(robust))
    bSEBC1=sqrt(diag(varKC))
    bSEBC2=sqrt(diag(varMD))
    bSEBC3=sqrt(diag(varFG))
    
    # to screen
    outbeta=cbind(beta,bSE,bSEBC0,bSEBC1,bSEBC2,bSEBC3)
    outalpha=cbind(alpha[1],alpha[2])
    colnames(outbeta)<-c("Estimate","stderr","BC0-stderr","BC1-stderr","BC2-stderr","BC3-stderr")
    colnames(outalpha)<-c("alpha0","alpha1")
    return(list(outbeta=outbeta,outalpha=outalpha,outphi=phi,niter=niter))
  }
  
  #####################################################################################
  # THE MAIN PROGRAM
  # SINGFLAG: THE ALGORITHM terminated due to singular MB covariance matrix
  # ROBFLAG: THE ALGORITHM terminated due to singular robust variance matrix
  # ALPFLAG: Check if the alpha estimates are within the plausible range
  # NPSDADJFLAG: Checks if the alpha adjustment factor matrix is positive definite, terminates if not
  #####################################################################################
  # reasons for non-results are identified and tallied
  SINGFLAG=0
  ROBFLAG=0
  ALPFLAG=0
  NPSDADJFLAG=0
  CONVFLAG=0
  
  # Fit the QLS Algorithm
  QLS_RES=FITQLS(y, X, n, m, t, maxiter, epsilon, SINGFLAG, ROBFLAG, ALPFLAG, NPSDADJFLAG)
  beta=QLS_RES$beta
  alpha=QLS_RES$alpha
  phi=QLS_RES$phi
  robust=QLS_RES$robust
  naive=QLS_RES$naive
  varMD=QLS_RES$varMD
  varKC=QLS_RES$varKC
  varFG=QLS_RES$varFG
  niter=QLS_RES$niter
  converge=QLS_RES$converge
  SINGFLAG=QLS_RES$SINGFLAG
  ROBFLAG=QLS_RES$ROBFLAG
  ALPFLAG=QLS_RES$ALPFLAG
  NPSDADJFLAG=QLS_RES$NPSDADJFLAG
  
  # Final Results
  if(converge==0 & SINGFLAG==0 & NPSDADJFLAG==0){CONVFLAG=1}
  if(converge==1 & ROBFLAG==0){
    return(RESULTS(beta, alpha, phi, robust, naive, varMD, varKC, varFG, niter, n))
  }
}
