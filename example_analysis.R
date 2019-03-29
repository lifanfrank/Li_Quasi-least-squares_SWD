#####################################################
# Illustrative use of the QLS/MAQLS program
# Analyzing a simulated trial data
#####################################################
set.seed(1234)
library(mvtnorm)

# change the dir to run the program
source(".../data_gen.R")
source(".../maqls_decay.R")

# design resources
n = 15               # number of clusters
t = 4                # number of periods
m = 20               # cohort size

# true parameter values
delta = 0.4          # effect size
alpha = c(0.1,0.8)  # values of tau=0.03, rho=0.2
beta=cumsum(c(0,0.1,0.1/2,0.1/(2^2)))   # period effect

# Create X matrix
X = NULL
trtSeq = matrix(0,t-1,t)
trtSeq[upper.tri(trtSeq)] = 1
g = n/(t-1)          # number of clusters per step
for(i in 1:(t-1)){
  for(j in 1:g){
    X = rbind(X,kronecker(rep(1,m),cbind(diag(t),trtSeq[i,])))}}

# Create cluster id for each obs
id = rep(1:n,each=m*t)
clsize = rep(m*t,n) # cluster size (across 4 periods)

# Generate outcome data
y = contGEN(n,m,t,delta,beta,alpha)
y = c(y)
cohortsize = clsize/t

# regular QLS
out1 = contQLS(y=y,X=X,id=id,n=clsize,m=cohortsize,t=t,maxiter=50,epsilon=0.001,alpadj=1)
#############################################################################
# Some output

# $outbeta
# Estimate     stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,] -0.11424762 0.09094165 0.07447470 0.07708864 0.07979432 0.07969879
# [2,] -0.02621135 0.09427154 0.08427490 0.08835641 0.09268996 0.09179680
# [3,]  0.04566648 0.10362113 0.09519600 0.09918412 0.10340060 0.10223850
# [4,]  0.01904666 0.11756289 0.11495690 0.12093539 0.12732893 0.12853151
# [5,]  0.47557729 0.07450269 0.08015067 0.08607153 0.09242982 0.08712111
# 
# $outalpha
# alpha0    alpha1
# [1,] 0.08376163 0.7948596
#############################################################################

out3 = contQLS(y=y,X=X,id=id,n=clsize,m=cohortsize,t=t,maxiter=50,epsilon=0.001,alpadj=3)
#############################################################################
# Some output

# $outbeta
# Estimate    stderr BC0-stderr BC1-stderr BC2-stderr BC3-stderr
# [1,] -0.11424762 0.0965101 0.07447470 0.07708864 0.07979432 0.07974604
# [2,] -0.02619479 0.1000070 0.08429865 0.08838006 0.09271351 0.09187804
# [3,]  0.04569960 0.1098318 0.09523978 0.09922933 0.10344721 0.10232256
# [4,]  0.01909633 0.1244953 0.11501014 0.12099158 0.12738809 0.12866898
# [5,]  0.47552762 0.0786440 0.08013356 0.08605364 0.09241111 0.08710925
# 
# $outalpha
# alpha0    alpha1
# [1,] 0.0994634 0.7969007
#############################################################################
