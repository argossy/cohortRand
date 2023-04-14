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

The package includes reference manual, sample data and a Vignette.

## Basic Usage

```r
library(cohortRand)


## simulation code for randomization assignments

## randomizating subject by cohort, then save the output by the date

## Cohort 1




```


## Contact
Kai Xia: kxia@med.unc.edu

