
library(demogR)

#####################################################
#                                                   #
#                     EXAMPLE 1                     #
#  FINITE NUMBER OF POPULATION TRANSITION MATRICES  #
#             STOCHASTIC PROCESS:  IID              #
#                                                   #
#####################################################
# MATRIX 1
A1<- matrix(0,5,5)
## SURVIVALS
surv1<- c(0.5, 0.7, 0.95, 0.95)
for(i in 1:4){A1[i+1,i]<- surv1[i]}
## FERTILITIES
A1[1,]<- c(0, 0, 1.5, 2, 3)
## EIGEN ANALYSIS
ea1<- eigen.analysis(A1)  
l1<- ea1$lambda1

# MATRIX 2
A2<- matrix(0,5,5)
## SURVIVALS5
surv2<- c(0.35, 0.6, 0.95, 0.95)
for(i in 1:4){A2[i+1,i]<- surv2[i]}
## FERTILITIES
A2[1,]<- c(0, 0, 1.2, 1.75, 2.5)
## EIGEN ANALYSIS
ea2<- eigen.analysis(A2)  
l2<- ea2$lambda1

# MATRIX 3
A3<- matrix(0,5,5)
## SURVIVALS5
surv3<- c(0.2, 0.5, 0.95, 0.95)
for(i in 1:4){A3[i+1,i]<- surv3[i]}
## FERTILITIES
A3[1,]<- c(0, 0, 1, 1, 1.2)
### EIGEN ANALYSIS
ea3<- eigen.analysis(A3)  
l3<- ea3$lambda1

## COMBINE MATRICES INTO LIST
A_list<- list(A1=A1, A2=A2, A3=A3)

## FIXED MATRIX PROBABILITIES OF OCCURING
p1<- 0.3
p2<- 0.5
p3<- 1-p1-p2

## GEOMETRIC MEAN APPROXIMATION FOR GROWTH RATE
l_ap0<- l1^p1*l2^p2*l3^p3 

# AVERAGE POPULATION TRANSITION MATRIX
A<- p1*A1+p2*A2+p3*A3
ea<- eigen.analysis(A)
lA<- ea$lambda1

#########################################################
# METHOD 1: MULTIPLE SAMPLE PATHS; POPULATION NUMBERS   #
# NOTE: CONVERGENCE CHECK NEEDED BUT NOT PERFORMED HERE #
#########################################################
## NUMBER OF SAMPLE PATHS
M<- 1000
## SAMPLE PATH LENGTH
tau<- 500
## GENERATE SAMPLE PATHS
omega<- matrix(sample(1:3, tau*M, replace=TRUE, prob=c(p1,p2,p3)), M, tau)
## USE INITIAL AND SUBSEQUENT POPULATION NUMBERS TO ESTIMATE 
## STOCHASTIC GROWTH RATE
N0<- rep(500,5)
Ltau1<- sapply(1:M, function(w)
{
  A_t<- A_list[[omega[w,1]]]
  for(i in 2:tau)
  {
    A_t<- A_t %*% A_list[[omega[w,i]]]
  }
  N_tau<- A_t %*% N0
  L_t<- sum(N_tau)/sum(N0)
  return(L_t)
})
l_ap1<- exp(mean(log(Ltau1)/tau))

#########################################################
# METHOD 2: MULTIPLE SAMPLE PATHS; POPULATION STRUCTURE #
# NOTE: USES THE PATHS GENERATED IN METHOD 1 ABOVE      #
#########################################################
## 2A: USE INITIAL AND SUBSEQUENT POPULATION STRUCTURE TO ESTIMATE
## STOCHASTIC GROWTH RATE
Y0<- N0/sum(N0)
Ltau2<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega[w,1]]]%*%Y0/sum(A_list[[omega[w,1]]]%*%Y0)
  L_t<- sum(A_list[[omega[w,1]]]%*%Y0)
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega[w,i]]]%*%Y_t/sum(A_list[[omega[w,i]]]%*%Y_t)
    L_t<- sum(A_list[[omega[w,i]]]%*%Y_t)*L_t
  }
  return(L_t)
})
l_ap2a<- exp(mean(log(Ltau2)/tau))

## 2B: SAME AS 2A BUT CHANGES THE ORDER OF THE PRODUCT-LOG TO LOG-SUM
logLtau<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega[w,1]]]%*%Y0/sum(A_list[[omega[w,1]]]%*%Y0)
  loglt_sum<- log(sum(A_list[[omega[w,1]]]%*%Y0))
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega[w,i]]]%*%Y_t/sum(A_list[[omega[w,i]]]%*%Y_t)
    loglt_sum<- log(sum(A_list[[omega[w,i]]]%*%Y_t))+loglt_sum
  }
  return(loglt_sum)
})
l_ap2b<- exp(mean(logLtau/tau))


### TEST CONVERGENCE
#### CHECK log(Ltau2) IS APPROXIMATELY NORMALLY DISTRIBUTED
# 2A:
hist(log(Ltau2))
qqnorm(log(Ltau2))
qqline(log(Ltau2), col="red", lty=3)
# 2B: 
hist(logLtau)
qqnorm(logLtau)
qqline(logLtau, col="red", lty=3)

#### CHECK THAT THE QUANTILES OF EARLIER AND LATER RUNS ARE SIMILAR
#### AND DISTRIBUTION IS NOT CHANGING -- MAY NEED LARGER TAU
##### RERUN ABOVE BUT STORE LOG OF ANNUAL GROWTH RATES 
##### (DONE WITH METHOD 2B BUT CAN ALSO BE DONE WITH 2A)
logl_t<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega[w,1]]]%*%Y0/sum(A_list[[omega[w,1]]]%*%Y0)
  loglt<- log(sum(A_list[[omega[w,1]]]%*%Y0))
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega[w,i]]]%*%Y_t/sum(A_list[[omega[w,i]]]%*%Y_t)
    loglt<- c(loglt, log(sum(A_list[[omega[w,i]]]%*%Y_t)))
  }
  return(loglt)
})
##### COMPARE QUANTILES FROM t=450 and t=500
qqplot(colSums(logl_t[1:450,]), colSums(logl_t))
abline(0,1, col="red", lty=3)

##### NOTE THAT THE GROWTH RATE OF A GIVEN RUN HAS CONVERGED
x<- sample(1:M, 1)
plot(1:tau, cumsum(logl_t[,x])/1:tau, type = "l")
abline(sum(logl_t[,x])/tau,0,col="red", lty=3)

##### LOG GROWTH CONVERGES TO INCREASING LINEARLY -- LONGER TAU MAY BE 
##### USEFUL FOR SEEING THIS/NEEDED FOR CONVERGENCE (SEE NEXT SECTION)
plot(1:tau, cumsum(logl_t[,x]), type = "l")
points(1:tau, logl_t[1,x]+sum(logl_t[,x])/tau*1:tau, type="l", col="red", lty=3)

 
## METHOD 2B MAY BE BETTER WHEN A LONG TIME STEP IS NEEDED FOR 
## CONVERGENCE IN SAMPLE PATHS
## SAMPLE PATH LENGTH
tau2<- 100000
## GENERATE SAMPLE PATHS
omega2<- sample(1:3, tau2, replace=TRUE, prob=c(p1,p2,p3))
## APPLY METHOD 2A
Y_t<- A_list[[omega2[1]]]%*%Y0/sum(A_list[[omega2[1]]]%*%Y0)
L_t<- sum(A_list[[omega2[1]]]%*%Y0)
for(i in 2:tau2)
{
  Y_t<- A_list[[omega2[i]]]%*%Y_t/sum(A_list[[omega2[i]]]%*%Y_t)
  L_t<- sum(A_list[[omega2[i]]]%*%Y_t)*L_t
}
L_t
exp(log(L_t)/tau2)
## APPLY METHOD 2B
Y_t<- A_list[[omega2[1]]]%*%Y0/sum(A_list[[omega2[1]]]%*%Y0)
logL_t<- log(sum(A_list[[omega2[1]]]%*%Y0))
for(i in 2:tau2)
{
  Y_t<- A_list[[omega2[i]]]%*%Y_t/sum(A_list[[omega2[i]]]%*%Y_t)
  logL_t<- log(sum(A_list[[omega2[i]]]%*%Y_t))+logL_t
}
logL_t
exp(logL_t/tau2)

## METHOD 2B WITH logL_t STORED AND ABOVE PLOT REVISITED
Y_t<- A_list[[omega2[1]]]%*%Y0/sum(A_list[[omega2[1]]]%*%Y0)
logL_t<- log(sum(A_list[[omega2[1]]]%*%Y0))
for(i in 2:tau2)
{
  Y_t<- A_list[[omega2[i]]]%*%Y_t/sum(A_list[[omega2[i]]]%*%Y_t)
  logL_t<- c(logL_t, log(sum(A_list[[omega2[i]]]%*%Y_t))+logL_t[i-1])
}
plot(1:tau2, logL_t, type = "l")
points(1:tau2, logL_t[1]+logL_t[tau2]/tau2*1:tau2, type="l", col="red", lty=3)
  # SLOPES SHOULD BE APPROXIMATELY THE SAME (BUT INTERCEPT MAY BE OFF)

##################################################
# METHOD 3: APPROXIMATION FORMULA                #
#           (TULJAPURKAR & CASWELL)              #
# NOTE: ONLY IS A GOOD APPROXIMATION WHEN VITAL  # 
#       RATE CVs<<1                              #
##################################################
## DEVIATIONS FROM AVERAGE POPULATION TRANSITION MATRIX
H1<- A1-A
H2<- A2-A
H3<- A3-A
# PAIRWISE COVARIANCES
## NOTE:  THE APPROACH USED BELOW WORKS ONLY IF MATRICES H1, H2, AND 
##        H3 HAVE NON-ZERO ELEMENTS IN ALL THE SAME PLACES (AND A IS
##        NON-ZERO HERE TOO). THIS WILL OFTEN BE THE CASE, AND THE 
##        APPROACH CAN BE MODIFIED TO ACCOUNT FOR DIFFERENCES.
indx<- which(H1!=0) 
COV1<- sapply(indx, function(x)
{
  out<- sapply(indx, function(y)
  {
    H1[x]/A[x]*H1[y]/A[y]
  })
  return(out)
})
COV2<- sapply(indx, function(x)
{
  out<- sapply(indx, function(y)
  {
    H2[x]/A[x]*H2[y]/A[y]
  })
  return(out)
})
COV3<- sapply(indx, function(x)
{
  out<- sapply(indx, function(y)
  {
    H3[x]/A[x]*H3[y]/A[y]
  })
  return(out)
})
COV<- p1*COV1+p2*COV2+p3*COV3
# ELASTICITY PRODUCTS
E<- sapply(indx, function(x)
{
  out<- sapply(indx, function(y)
  {
    ea$elasticities[x]*ea$elasticities[y]
  })
  return(out)
})
## COMPUTE THE APPROXIMATION
l_ap3<- exp(log(lA)-1/2*sum(E*COV))

## VITAL RATE CVs
sigs<- sqrt((A1-A)^2*p1+(A2-A)^2*p2+(A3-A)^2*p3)
CV<- sigs/A
CV[is.na(CV)]<- 0
max(CV) #0.2877
## VARIATION IN VITAL RATES COULD BE TOO LARGE TO USE THIS APPROXIMATION
## (ALTHOUGH IT APPEARS TO BE DOING OKAY)


#####################################################
#                                                   #
#                     EXAMPLE 2                     #
#  FINITE NUMBER OF POPULATION TRANSITION MATRICES  #
#         STOCHASTIC PROCESS:  MARKOV CHAIN         #
#                                                   #
#####################################################
## USE A1, A2, AND A3 FROM EXAMPLE 1 ABOVE
## TRANSITION MATRIX
P<- matrix(c(0.4, 0.5, 0.1, 0.3, 0.5, 0.2, 0.05, 0.55, 0.4), 3, 3)
## FIND STATIONARY PROBABILITIES FOR MATRICES A1, A2, AND A3
eaP<- eigen.analysis(P)
q<- eaP$stable.age
l_ap0<- l1^q[1]*l2^q[2]*l3^q[3]

#########################################################
# METHOD 1: MULTIPLE SAMPLE PATHS; POPULATION NUMBERS   #
# NOTE: CONVERGENCE CHECK NEEDED BUT NOT PERFORMED HERE #
#########################################################
## NUMBER OF SAMPLE PATHS
M<- 500
## SAMPLE PATH LENGTH
tau<- 5000
## GENERATE SAMPLE PATHS
omega<- t(sapply(1:M, function(x)
{
  w<-sample(1:3,1)
  for(i in 2:tau)
  {
    w<- c(w,sample(1:3, 1, prob=P[,w[i-1]]))
  }
  return(w)
}))
## USE INITIAL AND SUBSEQUENT POPULATION NUMBERS TO ESTIMATE 
## STOCHASTIC GROWTH RATE
N0<- rep(500,5)
Ltau1<- sapply(1:M, function(w)
{
  A_t<- A_list[[omega[w,1]]]
  for(i in 2:tau)
  {
    A_t<- A_t %*% A_list[[omega[w,i]]]
  }
  N_tau<- A_t %*% N0
  L_t<- sum(N_tau)/sum(N0)
  return(L_t)
})
l_ap1<- exp(mean(log(Ltau1)/tau))

#########################################################
# METHOD 2: MULTIPLE SAMPLE PATHS; POPULATION STRUCTURE #
# NOTE: USES THE PATHS GENERATED IN METHOD 1 ABOVE      #
#########################################################
## 2A: USE INITIAL AND SUBSEQUENT POPULATION STRUCTURE TO ESTIMATE
## STOCHASTIC GROWTH RATE
Y0<- N0/sum(N0)
Ltau2<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega[w,1]]]%*%Y0/sum(A_list[[omega[w,1]]]%*%Y0)
  L_t<- sum(A_list[[omega[w,1]]]%*%Y0)
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega[w,i]]]%*%Y_t/sum(A_list[[omega[w,i]]]%*%Y_t)
    L_t<- sum(A_list[[omega[w,i]]]%*%Y_t)*L_t
  }
  return(L_t)
})
l_ap2a<- exp(mean(log(Ltau2)/tau))

## 2B: SAME AS 2A BUT CHANGES THE ORDER OF THE PRODUCT-LOG TO LOG-SUM
logLtau<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega[w,1]]]%*%Y0/sum(A_list[[omega[w,1]]]%*%Y0)
  loglt_sum<- log(sum(A_list[[omega[w,1]]]%*%Y0))
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega[w,i]]]%*%Y_t/sum(A_list[[omega[w,i]]]%*%Y_t)
    loglt_sum<- log(sum(A_list[[omega[w,i]]]%*%Y_t))+loglt_sum
  }
  return(loglt_sum)
})
l_ap2b<- exp(mean(logLtau/tau))


### TEST CONVERGENCE
#### CHECK log(Ltau2) IS APPROXIMATELY NORMALLY DISTRIBUTED
# 2A:
hist(log(Ltau2))
qqnorm(log(Ltau2))
qqline(log(Ltau2), col="red", lty=3)
# 2B: 
hist(logLtau)
qqnorm(logLtau)
qqline(logLtau, col="red", lty=3)

#### CHECK THAT THE QUANTILES OF EARLIER AND LATER RUNS ARE SIMILAR
#### AND DISTRIBUTION IS NOT CHANGING -- MAY NEED LARGER TAU
##### RERUN ABOVE BUT STORE LOG OF ANNUAL GROWTH RATES 
##### (DONE WITH METHOD 2B BUT CAN ALSO BE DONE WITH 2A)
logl_t<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega[w,1]]]%*%Y0/sum(A_list[[omega[w,1]]]%*%Y0)
  loglt<- log(sum(A_list[[omega[w,1]]]%*%Y0))
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega[w,i]]]%*%Y_t/sum(A_list[[omega[w,i]]]%*%Y_t)
    loglt<- c(loglt, log(sum(A_list[[omega[w,i]]]%*%Y_t)))
  }
  return(loglt)
})
##### COMPARE QUANTILES FROM t=4900 and t=5000
qqplot(colSums(logl_t[1:4900,]), colSums(logl_t))
abline(0,1, col="red", lty=3)

##### COMPARE QUANTILES FROM t=4500 and t=5000
qqplot(colSums(logl_t[1:4500,]), colSums(logl_t))
abline(0,1, col="red", lty=3)

##### COMPARE QUANTILES FROM t=4000 and t=1000
qqplot(colSums(logl_t[1:4000,]), colSums(logl_t))
abline(0,1, col="red", lty=3)

##### NOTE THAT THE GROWTH RATE OF A GIVEN RUN HAS CONVERGED
x<- sample(1:M, 1)
plot(1:tau, cumsum(logl_t[,x])/1:tau, type = "l")
abline(sum(logl_t[,x])/tau,0,col="red", lty=3)

##### WHAT ABOUT LOG GROWTH?
## SAMPLE PATH LENGTH
tau2<- 100000
## GENERATE SAMPLE PATHS
omega2<- sample(1:3,1)
for(i in 2:tau2)
{
  omega2<- c(omega2,sample(1:3, 1, prob=P[,omega2[i-1]]))
}

## METHOD 2B WITH logL_t STORED AND LOG GROWTH PLOT
Y_t<- A_list[[omega2[1]]]%*%Y0/sum(A_list[[omega2[1]]]%*%Y0)
logL_t<- log(sum(A_list[[omega2[1]]]%*%Y0))
for(i in 2:tau2)
{
  Y_t<- A_list[[omega2[i]]]%*%Y_t/sum(A_list[[omega2[i]]]%*%Y_t)
  logL_t<- c(logL_t, log(sum(A_list[[omega2[i]]]%*%Y_t))+logL_t[i-1])
}
plot(1:tau2, logL_t, type = "l")
points(1:tau2, logL_t[1]+logL_t[tau2]/tau2*1:tau2, type="l", col="red", lty=3)
# SLOPE AT LATER TIMES SHOULD BE SIMILAR TO SLOPE OF DOTTED RED LINE
# GIVEN ENOUGH TIME HAS PASSED

##### 2C: CONSIDER 2B RAN LIKE EXAMPLE 2 USING STABLE MATRIX PROBABILITIES
omega2<- matrix(sample(1:3, replace=TRUE, M*tau, prob=q), M, tau)
logLtau2<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega2[w,1]]]%*%Y0/sum(A_list[[omega2[w,1]]]%*%Y0)
  loglt_sum<- log(sum(A_list[[omega2[w,1]]]%*%Y0))
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega2[w,i]]]%*%Y_t/sum(A_list[[omega2[w,i]]]%*%Y_t)
    loglt_sum<- log(sum(A_list[[omega2[w,i]]]%*%Y_t))+loglt_sum
  }
  return(loglt_sum)
})
l_ap2c<- exp(mean(logLtau2/tau))

### TEST CONVERGENCE
#### CHECK log(Ltau2) IS APPROXIMATELY NORMALLY DISTRIBUTED
# 2C:
hist(logLtau2)
qqnorm(logLtau2)
qqline(logLtau2, col="red", lty=3)

#### CHECK THAT THE QUANTILES OF EARLIER AND LATER RUNS ARE SIMILAR
#### AND DISTRIBUTION IS NOT CHANGING -- MAY NEED LARGER TAU
logl_t2<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega2[w,1]]]%*%Y0/sum(A_list[[omega2[w,1]]]%*%Y0)
  loglt<- log(sum(A_list[[omega2[w,1]]]%*%Y0))
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega2[w,i]]]%*%Y_t/sum(A_list[[omega2[w,i]]]%*%Y_t)
    loglt<- c(loglt, log(sum(A_list[[omega2[w,i]]]%*%Y_t)))
  }
  return(loglt)
})
##### COMPARE QUANTILES FROM t=4900 and t=5000
qqplot(colSums(logl_t2[1:4900,]), colSums(logl_t2))
abline(0,1, col="red", lty=3)

##### COMPARE QUANTILES FROM t=4500 and t=5000
qqplot(colSums(logl_t2[1:4500,]), colSums(logl_t2))
abline(0,1, col="red", lty=3)

##### COMPARE QUANTILES FROM t=4000 and t=1000
qqplot(colSums(logl_t2[1:4000,]), colSums(logl_t2))
abline(0,1, col="red", lty=3)

##### NOTE THAT THE GROWTH RATE OF A GIVEN RUN HAS CONVERGED
x<- sample(1:M, 1)
plot(1:tau, cumsum(logl_t2[,x])/1:tau, type = "l")
abline(sum(logl_t2[,x])/tau,0,col="red", lty=3)

##### LOG GROWTH
## SAMPLE PATH LENGTH
tau2<- 100000
## GENERATE SAMPLE PATHS
omega3<- sample(1:3, replace=TRUE, tau2, prob=q)

## METHOD 2B WITH logL_t STORED AND LOG GROWTH PLOT
Y_t<- A_list[[omega3[1]]]%*%Y0/sum(A_list[[omega3[1]]]%*%Y0)
logL_t<- log(sum(A_list[[omega3[1]]]%*%Y0))
for(i in 2:tau2)
{
  Y_t<- A_list[[omega3[i]]]%*%Y_t/sum(A_list[[omega3[i]]]%*%Y_t)
  logL_t<- c(logL_t, log(sum(A_list[[omega3[i]]]%*%Y_t))+logL_t[i-1])
}
plot(1:tau2, logL_t, type = "l")
points(1:tau2, logL_t[1]+logL_t[tau2]/tau2*1:tau2, type="l", col="red", lty=3)


##################################################
# METHOD 3: AN APPROXIMATION FORMULA EXISTS      # 
#           UNDER CERTAIN CONDITIONS; NEED TO    #
#            REFER TO TULJAPURKAR 1990, CH 12    #   
##################################################


#######################################################
#                                                     #
#  INFINITE NUMBER OF POPULATION TRANSITION MATRICES  #
#             STOCHASTIC PROCESS:  IID                #
#                                                     #
#######################################################
# METHODS 1-3 APPLY TO THIS CASE AS WELL
