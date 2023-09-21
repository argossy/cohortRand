# cohortRand, a randomization package for cohort randomization

## Introduction

For Randomized Controlled Trial of intervention studies, the participants are usually randomized as group or cohort to a few number of intervention arms. In this case, instead of randomizing subject by subject, the RCT randomization is done by group or cohort. Minimization method tends to minimize the differences among groups while controlling for a number of prognostic factors. cohortRand is designed to provide user friendly R implementation to perform cohort randomization using existing miniRand R package. It also provides summary statistics to track how each cohort has been randomized.


## Key features
1. Randomize RCT subject by cohort using minization method
2. Provide summary statistics for randomized subjets

## Install GitHub Version
To install `cohortRand` directly from GitHub, run

```r
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("argossy/cohortRand@main")
```


## Basic Usage
### Randomizing cohort 1

```r
s## simulation code for random assignments

## add cohort by cohort in each chunk, then save the output by the date
## Cohort 1

library(devtools)
library(gtsummary)
library(knitr)
library(Minirand) # minization method
library(cohortRand)


### run example

#fileout = sprintf('random_sequence_cohort_%s.csv', co_id)

# define parameters for the randomization


co_id = 1              ## cohort number
set.seed(2023 + co_id) ## random see for this cohort
method='Range'         ## minization method
ntrt <- 3              ## number of treatment types
ratio <- rep(1,ntrt)   ## assignment ratio
#trtseq <- 1:ntrt # treatment assignment number

rand.mat <- NULL



## simulating the information of this cohort
nsample <- round(runif(1,18,25)) # sample size of this cohort
  
## create random subject with covariates
covmat <- random_cov(nsample)

## check strata of this random cohort  
sum_cov(covmat)
  
## create input matrix with ID, covariates and empty randomization
rand.mat <- create_rand_mat(rand.mat, covmat)

## start cohort 1 randomization using minization method
rand.mat <- cohortRand(rand.mat, ntrt = ntrt, ratio = ratio,method=method)
#write.csv(rand.mat, file = fileout,row.names = FALSE,quote=FALSE)
  

### the randomization assignment
#cat(sprintf("#### Summary of all randomized subjects"))
sum_random(rand.mat,full.tab = TRUE,sum.tab = FALSE)

cat(sprintf("#### Summary of all randomized subjects"))
sum_random(rand.mat,full.tab = FALSE,sum.tab = TRUE)

```

### Randomizing cohort 2
```r
# define parameters for the randomization


co_id = 2              ## cohort number
set.seed(2023 + co_id) ## random seed for this cohort
method='Range'         ## minization method
ntrt <- 3              ## number of treatment types
ratio <- rep(1,ntrt)   ## assignment ratio
#trtseq <- 1:ntrt # treatment assignment number



## simulating the information of this cohort
nsample <- round(runif(1,18,25)) # sample size of this cohort
  
## create random subject with covariates
covmat <- random_cov(nsample, rand.mat$ID)

## check strata of this random cohort  
sum_cov(covmat)
  
## create input matrix with ID, covariates and empty randomization
rand.mat <- create_rand_mat(rand.mat, covmat)


## start cohort 1 randomization using minization method
rand.mat <- cohortRand(rand.mat, ntrt = ntrt, ratio = ratio,method=method)
#write.csv(rand.mat, file = fileout,row.names = FALSE,quote=FALSE)
  
rand.mat.co <- rand.mat[rand.mat$cohort == co_id,]
### the randomization assignment of cohort 2
sum_random(rand.mat.co,full.tab=TRUE)

cat(sprintf("#### Summary of all randomized subjects"))
sum_random(rand.mat,full.tab = FALSE,sum.tab = TRUE)







```

### Randomizing cohort 3



```r

co_id = 3              ## cohort number
set.seed(2023 + co_id) ## random seed for this cohort
method='Range'         ## minization method
ntrt <- 3              ## number of treatment types
ratio <- rep(1,ntrt)   ## assignment ratio
#trtseq <- 1:ntrt # treatment assignment number



## simulating the information of this cohort
nsample <- round(runif(1,18,25)) # sample size of this cohort
  
## create random subject with covariates
covmat <- random_cov(nsample, rand.mat$ID)

## check strata of this random cohort  
sum_cov(covmat)
  
## create input matrix with ID, covariates and empty randomization
rand.mat <- create_rand_mat(rand.mat, covmat)


## start cohort 1 randomization using minization method
rand.mat <- cohortRand(rand.mat, ntrt = ntrt, ratio = ratio,method=method)
#write.csv(rand.mat, file = fileout,row.names = FALSE,quote=FALSE)
  
rand.mat.co <- rand.mat[rand.mat$cohort == co_id,]
### the randomization assignment of cohort 2
sum_random(rand.mat.co,full.tab=TRUE)

cat(sprintf("#### Summary of all randomized subjects"))
sum_random(rand.mat,full.tab = FALSE,sum.tab = TRUE)


```







## Contact
Kai Xia: kxia@med.unc.edu

