#### Construct ROC curves to compare Cox and regularized regression models for clinical, radiomic and all
#### Only cluster features
# load libraries
library(survminer)
library(dplyr)
library(pec)
library(caret)
library(c060)
library(peperr)
library(survAUC)
library(ggplot2)

library(survival)
library(survivalROC)
library(glmnet)
library(nnet)
library(randomForestSRC)
library(survivalsvm)
library(birk)
library(bnnSurvival)
library(CoxBoost)

setwd("/Users/ani/Dropbox/2018Radiomics/RcodeFINAL/")
features <- read.csv("data431_surv.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

set.seed(123)
cutoff <- 365*3
ts <- seq(1, cutoff, 5)
nk <- 10

### Get PatientIDs and patients with bad masks
#patients = feature.dat$PatientID
#bad.patients = c("LUNG1-034", "LUNG1-040", "LUNG1-044", "LUNG1-068",
#                 "LUNG1-083", "LUNG1-084", "LUNG1-094", "LUNG1-096",
#                 "LUNG1-110", "LUNG1-137", "LUNG1-143", "LUNG1-144",
#                 "LUNG1-146", "LUNG1-164", "LUNG1-166", "LUNG1-167",
#                 "LUNG1-191", "LUNG1-192", "LUNG1-208", "LUNG1-212",
#                 "LUNG1-214", "LUNG1-222", "LUNG1-227")
#bad.ind.ft = which(feature.dat$PatientID %in% bad.patients)
#feature.dat.new = feature.dat.new[-bad.ind.ft,]

features <- features %>% 
  filter(!is.na(F1)) %>%
  filter(!is.na(Histology)) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(Overall.Stage)) %>%
  filter(!is.na(clinical.T.Stage))

id <- features$PatientID
features <- features %>% dplyr::select(-c(PatientID))

# categorical variables
# one-hot-encoding
features_cat <- features %>%
  dplyr::select(Overall.Stage, Histology, gender)
dmy <-  dummyVars("~.", data =features_cat) 
features_cat_ohe <- data.frame(predict(dmy, newdata = features_cat))
# other clinical variables
features_clc <- features %>%
  dplyr::select(clinical.T.Stage, Clinical.N.Stage, Clinical.M.Stage, age) %>%
  mutate(Tstage = as.numeric(clinical.T.Stage), Nstage = as.numeric(Clinical.N.Stage), Mstage = as.numeric(Clinical.M.Stage)) %>%
  dplyr::select(age, Tstage, Nstage, Mstage)

surv.time <- features$Survival.time
status <- features$deadstatus.event
surv_obj <- Surv(surv.time, status)
features_full = cbind(surv.time, status, features_clc, features_cat)

# matrix for regularized cox
features_mat <- as.matrix(cbind(features_clc, features_cat_ohe))

# create folds
index.kfold <- createFolds(seq(1, nrow(features)), k = nk, 
                           list = TRUE, returnTrain = FALSE)

auc.clinical.cox3.cv = c()
auc.clinical.lasso3.cv = c()
auc.clinical.eln3.cv = c()
auc.clinical.ridge3.cv = c()
auc.clinical.rf3.cv  = c()
auc.clinical.svm3.cv  = c()
auc.clinical.boost3.cv  = c()
auc.clinical.bnn3.cv = c()
auc.clinical.nnet3.cv = c()

for(cv.i in 1:10){
  
  ### Define training and test sets
  tr.ind = (1:dim(features)[1])[-index.kfold[[cv.i]]]
  
  ####### ********Model 1: Cox******** ########
  ### Run model with clinical variables only
  
  train.cox = features_full[tr.ind,]
  test.cox = features_full[-tr.ind,]
  fit.surv.clinical <- coxph(Surv(surv.time, status) ~ ., data=train.cox)
  ### Calculate the predicted values for the test set and add to the plot of ROC
  pr.clinical.test.cox = predict(fit.surv.clinical, newdata = test.cox)
  
  n.obs = dim(test.cox)[1]
  rc.clinical.cox3 = survivalROC(Stime = test.cox$surv.time, 
                                 status = test.cox$status, 
                                 marker = pr.clinical.test.cox, 
                                 predict.time = cutoff,
                                 span = 0.25*n.obs^(-0.20))
  auc.clinical.cox3.cv = c(auc.clinical.cox3.cv, rc.clinical.cox3$AUC)
  

  ####### Model 1: Clinical: Elnet#########
  ### Extract X, Y, and C from the training set
  X.train <- features_mat[tr.ind,]
  Y.train <- surv.time[tr.ind]
  C.train <- status[tr.ind]
  X.test <- features_mat[-tr.ind,]
  C.test <- status[-tr.ind]
  
  #### Cross validation for choosing the value of lambda
  cv.fit <- cv.glmnet(X.train, Surv(Y.train,C.train), alpha = 0.5, family = "cox")
  #### Elastic net
  fit <- glmnet(X.train, Surv(Y.train,C.train), family = "cox", alpha = 0.5, lambda = cv.fit$lambda.min)
  pr.clinical.test.elnet = predict(fit, X.test, type = "link", s = cv.fit$lambda.min)
  
  n.obs = dim(X.test)[1]
  rc.clinical.eln3 = survivalROC(Stime = surv.time[-tr.ind], 
                                 status = C.test, 
                                 marker = pr.clinical.test.elnet, 
                                 predict.time = cutoff,
                                 span = 0.25*n.obs^(-0.20))
  auc.clinical.eln3.cv = c(auc.clinical.eln3.cv, rc.clinical.eln3$AUC)

  
  #### Lasso
  cv.fit <- cv.glmnet(X.train, Surv(Y.train,C.train), alpha = 1, family = "cox")
  fit <- glmnet(X.train, Surv(Y.train,C.train), family = "cox", alpha = 1, lambda = cv.fit$lambda.min)
  
  #### Plot the ROC curve 
  pr.clinical.test.lasso = predict(fit, X.test, type = "link", s = cv.fit$lambda.min)
  n.obs = dim(X.test)[1]
  rc.clinical.lasso3 = survivalROC(Stime = surv.time[-tr.ind], 
                                   status = C.test, 
                                   marker = pr.clinical.test.lasso, 
                                   predict.time = cutoff,
                                   span = 0.25*n.obs^(-0.20))
  auc.clinical.lasso3.cv = c(auc.clinical.lasso3.cv, rc.clinical.lasso3$AUC)
  
  
  #### Ridge
  cv.fit <- cv.glmnet(X.train, Surv(Y.train,C.train), alpha = 0, family = "cox")
  fit <- glmnet(X.train, Surv(Y.train,C.train), family = "cox", alpha = 0, lambda = cv.fit$lambda.min)
  #### Plot the ROC curve 
  pr.clinical.test.ridge = predict(fit, X.test, type = "link", s = cv.fit$lambda.min)
  n.obs = dim(X.test)[1]
  rc.clinical.ridge3 = survivalROC(Stime = surv.time[-tr.ind], 
                                   status = C.test, 
                                   marker = pr.clinical.test.ridge, 
                                   predict.time = cutoff,
                                   span = 0.25*n.obs^(-0.20))
  auc.clinical.ridge3.cv = c(auc.clinical.ridge3.cv, rc.clinical.ridge3$AUC)
  

  ##### Model 1, random forest #####
  train.rf = as.data.frame(cbind(surv.time[tr.ind], status[tr.ind], features_mat[tr.ind,]))
  names(train.rf)[1] = "surv.time"
  names(train.rf)[2] = "status"
  test.rf = as.data.frame(cbind(surv.time[-tr.ind], status[-tr.ind], features_mat[-tr.ind,]))
  names(test.rf)[1] = "surv.time"
  names(test.rf)[2] = "status"
  
  fit.rf.clinical <- rfsrc(Surv(surv.time, status) ~ ., data = train.rf, nsplit = 8) 
  
  ### Calculate the predicted values for the test set and add to the plot of ROC
  pr.clinical.test.rf = predict(fit.rf.clinical, newdata = test.rf)
  n.obs = dim(test.rf)[1]
  rc.clinical.rf3 = survivalROC(Stime = test.rf$surv.time, 
                                status = test.rf$status, 
                                marker = pr.clinical.test.rf$predicted, 
                                predict.time = cutoff,
                                span = 0.25*n.obs^(-0.20))
  auc.clinical.rf3.cv = c(auc.clinical.rf3.cv, rc.clinical.rf3$AUC)
  
  
  ##### Model 1, SVM #####
  train.svm = train.rf
  test.svm = train.rf
  fit.svm.clinical <- survivalsvm(Surv(surv.time, status) ~ ., 
                                  data = train.svm, 
                                  type = "regression", gamma.mu = 1, 
                                  opt.meth = "quadprog", kernel = "lin_kernel") 
  
  ### Calculate the predicted values for the test set and add to the plot of ROC
  pr.clinical.test.svm = predict(fit.svm.clinical, newdata = test.svm)
  n.obs = dim(test.svm)[1]
  rc.clinical.svm3 = survivalROC(Stime = test.svm$surv.time, 
                                 status = test.svm$status, 
                                 marker = pr.clinical.test.svm$predicted, 
                                 predict.time = cutoff,
                                 span = 0.25*n.obs^(-0.20))
  auc.clinical.svm3.cv = c(auc.clinical.svm3.cv, rc.clinical.svm3$AUC)
  
  
  ##### Model 1: boosting #########
  x = data.matrix(features_mat[tr.ind,])
  cbfit <- CoxBoost(time = surv.time[tr.ind],
                    status = status[tr.ind],unpen.index=c(1,2),
                    x=x, stepno=100, penalty=100) 
  
  x.test = data.matrix(features_mat[-tr.ind,])
  pr.clinical.test.boost = predict(cbfit, x.test)
  
  n.obs = dim(x.test)[1]
  rc.clinical.boost3 = survivalROC(Stime = surv.time[-tr.ind], 
                                   status = status[-tr.ind], 
                                   marker = pr.clinical.test.boost, 
                                   predict.time = cutoff,
                                   span = 0.25*n.obs^(-0.20))
  auc.clinical.boost3.cv = c(auc.clinical.boost3.cv, rc.clinical.boost3$AUC)
  
  ##### Model 1: bnn #########
  train.bnn = train.rf
  test.bnn = test.rf
  fit.bnn.clinical <- bnnSurvival(Surv(surv.time, status) ~ ., train.bnn,
                                  k = 30, num_base_learners = 10, num_features_per_base_learner = 3)
  ### Calculate the predicted values for the test set and add to the plot of ROC
  pr.clinical.test.bnn <- predict(fit.bnn.clinical, test.bnn)
  n.obs = dim(test.bnn)[1]
  ind.bnn = which.closest(timepoints(pr.clinical.test.bnn), cutoff)
  rc.clinical.bnn3 = survivalROC(Stime = test.bnn$surv.time, 
                                 status = test.bnn$status, 
                                 marker = predictions(pr.clinical.test.bnn)[,ind.bnn], 
                                 predict.time = cutoff,
                                 span = 0.25*n.obs^(-0.20))
  auc.clinical.bnn3.cv = c(auc.clinical.bnn3.cv, rc.clinical.bnn3$AUC)
  
  ##### Model 2: nnet #########
  train.nnet = train.rf
  test.nnet = test.rf
  fit.nnet.clinical <- nnet(Surv(surv.time, status) ~ ., train.nnet,
                            size=1, maxit=500)
  ### Calculate the predicted values for the test set and add to the plot of ROC
  pr.clinical.test.nnet <- predict(fit.nnet.clinical, test.nnet)
  n.obs = dim(test.rf)[1]
  rc.clinical.nnet3 = survivalROC(Stime = test.nnet$surv.time, 
                                  status = test.nnet$status, 
                                  marker = pr.clinical.test.nnet[,2], 
                                  predict.time = cutoff,
                                  span = 0.25*n.obs^(-0.20))
  auc.clinical.nnet3.cv = c(auc.clinical.nnet3.cv, rc.clinical.nnet3$AUC)
  
print(cv.i)
}

### Both feature.dat.new and feature.dat.new.el have clinial and radiomic features
### The categorical clinical features are saved as numeric in feature.dat.new.el 

results.clinical3 = rbind(auc.clinical.cox3.cv, auc.clinical.lasso3.cv, 
                      auc.clinical.ridge3.cv, auc.clinical.eln3.cv,
                      auc.clinical.rf3.cv, auc.clinical.svm3.cv,
                      auc.clinical.boost3.cv, auc.clinical.bnn3.cv)
round(apply(results.clinical3, 1, mean),3)
round(apply(results.clinical3, 1, sd),3)


