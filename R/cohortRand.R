### features of cohort randomization
## input
#### existing and new subjects
####     with covariates/strata and assigned and unassigned randomization
####.   row names as subject ID
#### ratio of trt
#### number of treatment groups
#### cohort #


### procedure
## 1) when no subject has been assigned, randomly assign the first subject
## 2) the future subject will be sequentially added to the assignment by minirand
## 3) create a randomization ID as 'Random_XXX'
## 4) output current imbalance and summary table of randomization by cohort



cohortRand <- function(rand.mat,ntrt=2, ratio=rep(1,ntrt), trtseq= 1:ntrt,method='Range') {
  # ratio = c(1,1);ntrt = 2;
  #trtseq <- 1:ntrt
  # if(is.null(trtNames)){
  #   trtNames <- LETTERS[trtseq]
  # }

  ## assume first column is ID and coded as 'Name_XXX'
  id_data <- rand.mat[,1]



  ##


  ## assume covariates has CV in names
  ## extract strata
  name_cv = grep('CV',colnames(rand.mat),value=TRUE)
  covmat <- rand.mat[,name_cv]
  colnames(covmat) = gsub('CV.','', colnames(covmat))

  # use balance weight
  wi <- 1/length(name_cv)
  covwt <- rep(wi, length(name_cv))

  ## extract current randomization info
  #res <- rand.mat$random


  ## check whether this is the first cohort
  #cond1 <- sum(!is.na(rand.mat$cohort)) == 0
  cond2 <- sum(!is.na(rand.mat$random)) == 0

  if(cond2){
    # if this is the 1st cohort with no subject has been assigned
    # randomly assign the first subject
    n_cur_cohort <- dim(covmat)[1]
    cohort_id = 1
    cur_idx = 1
    rand.mat$random[cur_idx] <- sample(trtseq, 1, replace = TRUE, prob = ratio/sum(ratio))
    rand.mat$cohort[cur_idx] <- cohort_id

    # for 2 to N subjects, use minirand to randomize
    if(n_cur_cohort > 1){
      res <- rand.mat$random
      res[is.na(res)] = 100
      cur_idx <- cur_idx + 1
      id_start = cur_idx

      for( cur_idx in id_start:n_cur_cohort){
        res[cur_idx] <- Minirand(covmat=covmat, cur_idx, covwt=covwt, ratio=ratio,
                           ntrt=ntrt, trtseq=trtseq, method=method, result=res, p = 0.9)
        rand.mat$cohort[cur_idx] = cohort_id
      }
      res[res == 100] = NA
      rand.mat$random = res
    }
      ## randomization ID
      random.id <- paste(1:n_cur_cohort)
      rand.mat$random.id = random.id


  } else{
    ## if this is not first cohort
    # compute the new subject
    n_cur_cohort <- sum(is.na(rand.mat$cohort))

    ## previous cohort + current cohort
    n_total <- dim(rand.mat)[1]

    # calculate the cohort ID of this new cohort
    cohort_id = max(rand.mat$cohort, na.rm = TRUE) + 1

    ## assign the starting id as the first subject without cohort id
    id_start = which(is.na(rand.mat$cohort))[1]

    res <- rand.mat$random
    res[is.na(res)] = 100

    for( cur_idx in id_start:n_total){
      res[cur_idx] <- Minirand(covmat=covmat, cur_idx, covwt=covwt, ratio=ratio,
                               ntrt=ntrt, trtseq=trtseq, method="Range", result=res, p = 0.9)
      rand.mat$cohort[cur_idx] = cohort_id
      rand.mat$random.id[cur_idx] = cur_idx
    }
    res[res == 100] = NA

    rand.mat$random = res


  }

  return(rand.mat)

}




### generate new rand.mat using new covariates
### assume covmat has ID and covariates
create_rand_mat <- function(rand.mat, covmat){
  ## add columns for covmat
  ## create input matrix with both covartes and randomization
  #cohort = rep(NA, dim(covmat)[1])
  #random = cohort
  #random.id = cohort


  ## add randomization information
  dat.mat <- data.frame(covmat, cohort=NA, random=NA, random.id=NA)

  rand.mat <- rbind(rand.mat, dat.mat)

  # #covmat <- data.frame()
  # if(is.null(rand.mat)){
  #   rand.mat <- covmat
  # } else{
  #   rand.mat <- rbind(rand.mat, covmat)
  # }
  #
  return(rand.mat)
}


##
sum_random <- function(rand.mat){
  rand.mat.out <- rand.mat

  rand.mat.out$cohort = as.factor(rand.mat.out$cohort)
  colnames(rand.mat.out) <- gsub('CV.','', colnames(rand.mat.out))

  trtNames <- LETTERS[1:ntrt]
  rand.mat.out$random <- trtNames[rand.mat.out$random ]

  # output file
  write.csv(rand.mat, file = fileout,row.names = FALSE,quote=FALSE)

  ### output into PDF
  print(kable(rand.mat.out,caption=sprintf('Cohort %s randomization asignment',co_id)))
  #cat('\\pagebreak')

  tb1 <- rand.mat.out %>%
    select(-c(ID,random.id)) %>%
    tbl_summary(by=random) %>%
    modify_header(label ~ "**Strata**") %>%
    modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment Received**") %>%
    modify_caption("**Summary table of randomization assignment**") %>%
    bold_labels()

  tb1 %>% as_gt()
}


sum_cov <- function(covmat){
  rand.mat.out <- covmat

  #rand.mat.out$cohort = as.factor(rand.mat.out$cohort)
  colnames(rand.mat.out) <- gsub('CV.','', colnames(rand.mat.out))

  #trtNames <- LETTERS[1:ntrt]
  #rand.mat.out$random <- trtNames[rand.mat.out$random ]

  # output file
  #write.csv(rand.mat, file = fileout,row.names = FALSE,quote=FALSE)

  ### output into PDF
  #print(kable(rand.mat.out,caption=sprintf('Cohort %s randomization asignment',co_id)))
  #cat('\\pagebreak')

  tb1 <- rand.mat.out %>%
    select(-c(ID)) %>%
    tbl_summary() %>%
    modify_header(label ~ "**Strata**") %>%
    #modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment Received**") %>%
    modify_caption("**Summary table of strata in this cohort**") %>%
    bold_labels()

  tb1 %>% as_gt()
}


#
#
# ### run example
# library(Minirand)
#
# set.seed(2023)
# ntrt <- 2
# nsample <- 12
# trtseq <- 1:ntrt
# ratio <- rep(1,ntrt)
# c1 <- sample(c('M', 'F'), nsample, replace = TRUE, prob = c(0.4, 0.6))
# c2 <- sample(c('H','L'), nsample, replace = TRUE, prob = c(0.3, 0.7))
# covmat <- cbind(c1, c2) # generate the matrix of covariate factors for the subjects
# # label of the covariates
# colnames(covmat) = c("CV.Gender", "CV.Risk")
#
# ## create input matrix with both covartes and randomization
# cohort = rep(NA, dim(covmat)[1])
# random = rep(NA, dim(covmat)[1])
# ID = paste("ID_",1:nsample,sep='')
# dat.mat <- data.frame(ID,covmat, cohort, random)
#
# ## start cohort 1 randomization using minization method
# rand.mat <- cohortRand(dat.mat, ntrt = 2)
#
# trt1 <- rand.mat$random
# balance1 <- randbalance(trt1, covmat, ntrt, trtseq)
# balance1
#
# covwt <- c(1/2,1/2)
# totimbal(trt = trt1, covmat = covmat, covwt = covwt,
#          ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range")
#
#
# ## cohort 2
# nsample = 10
# c1 <- sample(c('M', 'F'), nsample, replace = TRUE, prob = c(0.4, 0.6))
# c2 <- sample(c('H','L'), nsample, replace = TRUE, prob = c(0.3, 0.7))
# covmat <- cbind(c1, c2) # generate the matrix of covariate factors for the subjects
# # label of the covariates
# colnames(covmat) = c("CV.Gender", "CV.Risk")
#
# ## create input matrix with both covartes and randomization
# cohort = rep(NA, dim(covmat)[1])
# random = cohort
# random.id = cohort
# id_start = dim(rand.mat)[1] + 1
# id_end = id_start + nsample - 1
# ID = paste('ID_', id_start:id_end, sep='')
#
# dat.mat <- data.frame(ID,covmat, cohort, random,random.id)
#
# rand.mat <- rbind(rand.mat, dat.mat)
# rand.mat <- cohortRand(rand.mat, ntrt = 2)
#
# trt1 <- rand.mat$random
# covmat <- rand.mat[,2:3]
# balance1 <- randbalance(trt1, covmat, ntrt, trtseq)
# balance1
#
# covwt <- c(1/2,1/2)
# totimbal(trt = trt1, covmat = covmat, covwt = covwt,
#          ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range")
#
#
#
#
# ## cohort 3
# nsample = 15
# c1 <- sample(c('M', 'F'), nsample, replace = TRUE, prob = c(0.4, 0.6))
# c2 <- sample(c('H','L'), nsample, replace = TRUE, prob = c(0.3, 0.7))
# covmat <- cbind(c1, c2) # generate the matrix of covariate factors for the subjects
# # label of the covariates
# colnames(covmat) = c("CV.Gender", "CV.Risk")
#
# ## create input matrix with both covartes and randomization
# cohort = rep(NA, dim(covmat)[1])
# random = cohort
# random.id = cohort
# id_start = dim(rand.mat)[1] + 1
# id_end = id_start + nsample - 1
# ID = paste('ID_', id_start:id_end, sep='')
#
# dat.mat <- data.frame(ID,covmat, cohort, random,random.id)
#
# rand.mat <- rbind(rand.mat, dat.mat)
# rand.mat <- cohortRand(rand.mat, ntrt = 2)
#
# trt1 <- rand.mat$random
# covmat <- rand.mat[,2:3]
# balance1 <- randbalance(trt1, covmat, ntrt, trtseq)
# balance1
#
# covwt <- c(1/2,1/2)
# totimbal(trt = trt1, covmat = covmat, covwt = covwt,
#          ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range")
#
#
#
#
#
#
#
#
#
#
#
# ### output
# library(gtsummary)
#
# rand.mat.out <- rand.mat
#
# ## remove ID column
# idx_rm = which(colnames(rand.mat.out) %in% c('ID','random.id'))
# rand.mat.out = rand.mat.out[,-idx_rm]
#
# colnames(rand.mat.out) <- gsub('CV.','', colnames(rand.mat.out))
# #trtNames <- paste('Treatment',LETTERS[1:ntrt])
# trtNames <- LETTERS[1:ntrt]
# rand.mat.out$random <- trtNames[rand.mat.out$random ]
#
# rand.mat.out %>%
#   tbl_summary(by=random) %>%
#   modify_header(label ~ "**Strata**") %>%
#   modify_spanning_header(c("stat_1", "stat_2") ~ "**Treatment Received**") %>%
#   modify_caption("**Table 1. Randomization assignment**") %>%
#   bold_labels()
#
