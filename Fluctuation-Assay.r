#### MSS Maximum Likelihood example from paper [2]
rDt <-c(0:7,9)
#[1] 0 1 2 3 4 5 6 7 9
cDt <-c(11,17,12,3,4,1,1,2,1)
#[1] 11 17 12 3 4 1 1 2 1

dt <- rep(0:9,c(11,17,12,3,4,1,1,2,0,1))
#[1]   0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[29] 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 4 4 4 4 5 6 7 7 9

conditionalProb_r = function(m,r) {
   #inputs: m, and r
   partials = numeric(r+1)
   partials[1] = exp(-1*m)
   if (r>0) {
       if (r==1) {return((m/r)*partials[1]/(r+1))}
           for (i in 1:(r)) {
               partialSum = 0
               for (j in 1:i) {
                   partialSum=partialSum+partials[j]/(i-(j-1)+1)
               }
               partials[i+1] = (m/i)*partialSum
           }
           return(partials[r+1])
   }
   return(partials[1])
}
  
loglik=function(m,rVec,rCard){ #input m, unique r and counts of r
   objfn = 0
   for (i in 1:length(rVec)){
       objfn = objfn + rCard[i]*log(conditionalProb_r(m,rVec[i]))
   }
   return(objfn)
}
  
f=function(z) {-1*loglik(z,rDt,cDt)}
  
m <- optimize(f,c(0,5))$minimum
m
#[1] 1.106733
  
#Generate CI with bootstrap
######################Confidence Intervals############
nBootstraps = 1000
  
##boot1 - draw with relacement from data
C=length(dt)
bReplicates1 <- numeric(nBootstraps)
  
set.seed(4) # for reproducibility
for (b in 1:nBootstraps) {
   indices <- sample(1:C,C,replace=TRUE)
   bRep <- dt[indices]
   # get unique r values
   rBRep <- as.numeric(names(table(bRep)))
    # get counts for each r value
   cBRep <- as.numeric(table(bRep))
   #build log likelihood function - negative for minimisation
   f=function(z) {-1*loglik(z,rBRep,cBRep)}
   #obtain maximum likleihood estimate for bootstrap replicate
   bReplicates1[b] <- optimize(f,c(0,5))$minimum
}
quantile(bReplicates1,c(0.025,0.975))
#               2.5%         97.5%
#     0.8984149 1.3579015
  
#Compare to normal approximation
m
#[1] 1.106733
sigma = 1.225*(m^-0.315)/sqrt(C)
u=log(m)+1.96*sigma*(exp(1)^(1.96*sigma))^-0.315
l=log(m)-1.96*sigma*(exp(1)^(1.96*sigma))^-0.315
exp(c(l,u))
#[[1] 0.8270193 1.4810507

# Assuming nt=320000 comput mu and CI for mu
nt=320000
mu=m/nt
mu
# 3.45854e-06
quantile(bReplicates1/nt,c(0.025,0.975))
#           2.5%        97.5% 
#   2.807547e-06 4.243442e-06
