library(readxl)
library(survival)
library(glmnet)
library(survminer)
library(dplyr)
library(rms)

#Define functions
gm_mean = function(x, na.rm=TRUE){
  #Calculates geometric mean
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
calc_BagRiskScore = function(ns.trim, all.coefs) {
  #Function calculates bagged ridge recalibrated risk scores. 
  #ns.trim: positive control normalized counts exported from nSolver
  #all.coefs: list containing 500 ridge regression submodels comprising the bootstrap aggregated model
  wanted.genes.housekeeping = c('ACTB','CASC3','GUSB','KIAA1715','MRPL19','POP4','TTC21B')
  wanted.genes = c('ARID1B','CCL21','CCN1','CCND2','CD3E','CDC20','CDK6','CDKN2A','CDKN2C','CHEK1','CKS2','COL1A1','ESR1','EZH2','FBLIM1','FGFR4','GAS1','IFNGR1','IGF2','KDR','KIF20A','KRT14','LINC02593','MDM4','MMP9','MUTYH','MYBL1','PGK1','PGR','PIM1','SPOP','TAGLN','TMEM30B','USF1')
  preds = c()
  for (i in seq(1,length(all.coefs))) {
    preds = cbind(preds,predict(all.coefs[[i]], newx = log2(data.matrix((ns.trim[,wanted.genes]))/apply(ns.trim[,wanted.genes.housekeeping],1,gm_mean))))
  }
  newy=rowMeans(preds)
  return(newy)
}

#Data indexing
dat = dat[dat$HK.7best.geom.mean/dat$Pos.ctrl.geom.mean > 0.20,] #Remove low signal
traindat = subset(dat,dat$Site == "UCSF") #Training data
dat.mat = data.matrix(traindat[,4:103]) #Indexing 100 genes
dat.lffr = data.matrix(traindat[,c('LF','LFFR')]) #Indexing outcomes

#Lasso model
set.seed(1)
colnames(dat.lffr) = c("status","time")
fit <- glmnet(dat.mat, dat.lffr, family = "cox")
cvfit <- cv.glmnet(dat.mat, dat.lffr, family = "cox", type.measure = "C")
plot(cvfit)
traindat$Risk.Score = dat.mat%*%coef(fit, s = cvfit$lambda.1se) #Multiply gene values by model coefficients
traindat$Risk.Score = (traindat$Risk.Score - (-4))/(4-(-4)) #Linearly rescale by upper and lower bounds

#Nested cut off determination
#LF
opt.cut.1 = surv_cutpoint(traindat, time = "LFFR", event = "LF",variables = c("Risk.Score"))
upper = subset(traindat,Risk.Score>opt.cut.1)
lower = subset(traindat,Risk.Score<=opt.cut.1)
opt.cut.upper.LF = surv_cutpoint(upper, time = "LFFR", event = "LF",variables = c("Risk.Score"))
opt.cut.lower.LF = surv_cutpoint(lower, time = "LFFR", event = "LF",variables = c("Risk.Score"))
#OS
opt.cut.1 = surv_cutpoint(traindat, time = "OS", event = "VS",variables = c("Risk.Score"))
upper = subset(traindat,Risk.Score>opt.cut.1)
lower = subset(traindat,Risk.Score<=opt.cut.1)
opt.cut.upper = surv_cutpoint(upper, time = "OS", event = "VS",variables = c("Risk.Score"))
opt.cut.lower = surv_cutpoint(lower, time = "OS", event = "VS",variables = c("Risk.Score"))

#Ridge recalibration
wanted.genes.housekeeping = c('ACTB','CASC3','GUSB','KIAA1715','MRPL19','POP4','TTC21B')
wanted.genes = c('ARID1B','CCL21','CCN1','CCND2','CD3E','CDC20','CDK6','CDKN2A','CDKN2C','CHEK1','CKS2','COL1A1','ESR1','EZH2','FBLIM1','FGFR4','GAS1','IFNGR1','IGF2','KDR','KIF20A','KRT14','LINC02593','MDM4','MMP9','MUTYH','MYBL1','PGK1','PGR','PIM1','SPOP','TAGLN','TMEM30B','USF1')
nb = 500 #Number of submodels
for (i in seq(1,nb)){
    #Bootstrap resampling with replacement
    #ns.trim is the positive control normalized data matrix depending on context
    train_ind <- sample(seq_len(nrow(ns.trim)), size = dim(ns.trim,1),replace = TRUE)
    x_var = log2(data.matrix((ns.trim[train_ind,wanted.genes]))/apply(ns.trim[train_ind,wanted.genes.housekeeping],1,gm_mean))
    y_var <- traindat$Risk.Score[match(ns.trim[train_ind,]$ID,traindat$ID)]
    lambda_seq <- 10^seq(2, -2, by = -.1)
    #Using glmnet function to build  ridge regression
    fit <- glmnet(x_var, y_var, alpha = 0, lambda  = lambda_seq)
    ridge_cv <- cv.glmnet(x_var, y_var, alpha = 0, lambda = lambda_seq)
    best_lambda <- ridge_cv$lambda.min
    best_ridge <- glmnet(x_var, y_var, alpha = 0, lambda = best_lambda)
    all.coefs[[i]] = best_ridge #Resulting list containing 500 submodels
}

#Validation risk score calculation
frozen.risks = calc_BagRiskScore(subset(ns.trim.val,Tissue_Type=='Frozen'),all.coefs.frozen)
ffpe.risks = calc_BagRiskScore(subset(ns.trim.val,Tissue_Type=='FFPE'),all.coefs.ffpe)
val.RiskScores = c(frozen.risks,ffpe.risks)
ns.trim.val = rbind(subset(ns.trim.val,Tissue_Type=='Frozen'),subset(ns.trim.val,Tissue_Type=='FFPE'))
ns.trim.val$RiskScores = val.RiskScores

