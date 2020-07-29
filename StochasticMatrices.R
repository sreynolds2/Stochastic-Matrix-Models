
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

##################################################
# METHOD 1: MULTIPLE SAMPLE PATHS                #
# NOTE: NEED TO CHECK CONVERGENCE FOR ALL CASES  #
##################################################
## NUMBER OF SAMPLE PATHS
M<- 1000
## SAMPLE PATH LENGTH
tau<- 500
## GENERATE SAMPLE PATHS
omega<- matrix(sample(1:3, tau*M, replace=TRUE, prob=c(p1,p2,p3)), M, tau)
## 1A: USE INITIAL AND SUBSEQUENT POPULATION NUMBERS
N0<- rep(500,5)
Ltau_A<- sapply(1:M, function(w)
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
l_ap1a<- exp(mean(log(Ltau_A)/tau))

## 1B: USE INITIAL AND SUBSEQUENT POPULATION STRUCTURE
Y0<- N0/sum(N0)
Ltau_B<- sapply(1:M, function(w)
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
l_ap1b<- exp(mean(log(Ltau_B)/tau))

## 1C: SAME AS 1B BUT CHANGES THE ORDER OF THE PRODUCT-LOG TO LOG-SUM
logltau_sum<- sapply(1:M, function(w)
{
  Y_t<- A_list[[omega[w,1]]]%*%Y0/sum(A_list[[omega[w,1]]]%*%Y0)
  logl_t<- log(sum(A_list[[omega[w,1]]]%*%Y0))
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega[w,i]]]%*%Y_t/sum(A_list[[omega[w,i]]]%*%Y_t)
    logl_t<- log(sum(A_list[[omega[w,i]]]%*%Y_t))+logl_t
  }
  return(logl_t)
})
l_ap1c<- exp(mean(logltau_sum/tau))

##################################################
# METHOD 2: SINGLE SAMPLE PATH                   #
# ABOVE SHOULD BE MORE RELIABLE BUT IF YOU       #
# NEEDED TO CUT OUT REPLICATES FOR SOME REASON   #
# YOU COULD TRY THE FOLLOWING                    #
##################################################
tau<- 20000
test<- sapply(1:100, function(x)
{
  omega<- sample(1:3, tau, replace=TRUE, prob=c(p1,p2,p3))
  Y_t<- A_list[[omega[1]]]%*%Y0/sum(A_list[[omega[1]]]%*%Y0)
  L_t<- sum(A_list[[omega[1]]]%*%Y0)
  for(i in 2:tau)
  {
    Y_t<- A_list[[omega[i]]]%*%Y_t/sum(A_list[[omega[i]]]%*%Y_t)
    L_t<- c(L_t, sum(A_list[[omega[i]]]%*%Y_t))
  }
  return(log(prod(L_t))/tau)
})

plot(1:tau, log(cumprod(L_t)), type = "l")
points(1:tau, 8*log(L_t[1])+log(prod(L_t))/tau*1:tau, type="l", col="red", lty=3)

plot(1:tau, log(cumprod(L_t))/1:tau, type = "l")
abline(log(prod(L_t))/tau,0,col="red", lty=3)

##################################################
# METHOD 2: APPROXIMATION FORMULA                #
# NOTE: ONLY IS A GOOD APPROXIMATION WHEN VITAL  # 
#       RATE CVs<<1                              #
##################################################
## DEVIATIONS FROM AVERAGE POPULATION TRANSITION MATRIX
H1<- A1-A
H2<- A2-A
H3<- A3-A
# PAIRWISE COVARIANCES
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
## (ALTHOUGH IT APPEARS TO BE DOING WELL)



#####################################################
#                                                   #
#  FINITE NUMBER OF POPULATION TRANSITION MATRICES  #
#         STOCHASTIC PROCESS:  MARKOV CHAIN         #
#                                                   #
#####################################################

## USE A1, A2, AND A3 ABOVE
## TRANSITION MATRIX



#######################################################
#                                                     #
#  INFINITE NUMBER OF POPULATION TRANSITION MATRICES  #
#             STOCHASTIC PROCESS:  IID                #
#                                                     #
#######################################################

