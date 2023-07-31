#©2023 The Regents of the University of California. All Rights Reserved

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

#Generate pairwise likelihood-ratio test heatmap comparing classification systems

#Grade = WHO 2016 grade, WHO.2021 = WHO 2021 grade, Chen = targeted gene expression risk score, Maas = DNA methylation based "integrated score", Sahm = DNA methylation families, Driver = "integrated grade", Olar.Lasso = methylation probe based scoree, Patel = "gene expresion types", 4_subgroups = DNA subgroups, DNA_methylation_groups = DNA groups.

models = c("Chen","WHO.2021","Sahm","Driver","Olar.Lasso","`4_subgroups`","Maas","Patel","DNA_methylation_groups","Grade") 

df.heatmap = subset(df,!is.na(df$Sahm)&!is.na(df$Chen)&!is.na(df$Driver)&!is.na(df$Patel)&!is.na(df$WHO.2021)) #Obtain complete cases with data available for all comparison classification systems, N=290
confusion_pvals = matrix(nrow=10,ncol=10)

for (i in seq(1,10)){
  for (j in seq(1,10)){
    modeli_j = as.formula(paste("Surv(LFFR,LF)~",models[i],"+",models[j],sep=""))
    modeli = as.formula(paste("Surv(LFFR,LF)~",models[i],sep=""))
    model1 = coxph(modeli_j,data=df.heatmap)
    model2 = coxph(modeli,data=df.heatmap)
    confusion_pvals[i,j] = lrtest(model1, model2)$`Pr(>Chisq)`[2] #Likelihood ratio test
  }
}

pheatmap(-log2(matrix(p.adjust((confusion_pvals),method = 'hochberg'),nrow = 10,ncol=10)),cluster_rows=FALSE, cluster_cols=FALSE,scale = 'none',cellwidth = 16,cellheight = 16) #heatmap

#Brier error scores
mod0=coxph(Surv(LFFR,LF)~Chen,data=df.heatmap,x=TRUE)
mod1=coxph(Surv(LFFR,LF)~Chen.groups,data=df.heatmap,x=TRUE)
mod2=coxph(Surv(LFFR,LF)~Driver.continuous,data=df.heatmap,x=TRUE)
mod3=coxph(Surv(LFFR,LF)~Driver.group,data=df.heatmap,x=TRUE)
mod4=coxph(Surv(LFFR,LF)~Maas,data=df.heatmap,x=TRUE)
mod5=coxph(Surv(LFFR,LF)~Sahm,data=df.heatmap,x=TRUE)
mod6=coxph(Surv(LFFR,LF)~Olar.Lasso,data=df.heatmap,x=TRUE)
mod7=coxph(Surv(LFFR,LF)~DNA_methylation_groups,data=df.heatmap,x=TRUE)
mod8=coxph(Surv(LFFR,LF)~`4_subgroups`,data=df.heatmap,x=TRUE)
mod9=coxph(Surv(LFFR,LF)~WHO.2021,data=df.heatmap,x=TRUE)
mod10=coxph(Surv(LFFR,LF)~Grade,data=df.heatmap,x=TRUE)
mod11=coxph(Surv(LFFR,LF)~Patel,data=df.heatmap,x=TRUE)

brier <- pec(list("Chen" = mod0, "Chen discrete"=mod1, "Driver"=mod2, "Driver discrete" = mod13, "Maas" = mod4, "Sahm"=mod5, "Olar"=mod6, "Choudhury"=mod7,"Choudhury-Chen"=mod8,"Patel"=mod11,"WHO 2021"=mod9,"WHO 2016"=mod10),data=df.heatmap,formula=Surv(LFFR,LF)~1,splitMethod = 'Boot632',keep.matrix = TRUE,B=1000,M=250)
print(brier)
plot(brier,xlim=c(0,10),ylim = c(0,0.25),legend = FALSE)

#Bootstrap delta-AUC, variable N depending on pairwise complete cases

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(Driver)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$Driver,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(WHO.2021)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$WHO.2021,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(Grade)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$Grade,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(Sahm)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$Sahm,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(Maas)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$Maas,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(Patel)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$Patel,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(Olar.Lasso)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$Olar.Lasso,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(`4_subgroups`)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$`4_subgroups`,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

b = 1000 #resample 1000 times
deltas = c()
df.2 = subset(df,!is.na(DNA_methylation_groups)&!is.na(Chen))
for (i in seq(1,b)){
  samp_inds = sample(seq(1:length(df.2$ID)),replace = TRUE)
  dt = df.2[samp_inds,]
  s1 = survivalROC(dt$LFFR,dt$LF,dt$Chen,method='KM',predict.time=5)
  d = survivalROC(dt$LFFR,dt$LF,dt$DNA_methylation_groups,method='KM',predict.time=5)
  delta = s1$AUC - d$AUC
  deltas = c(deltas,delta)
}

median(deltas) 
quantile(deltas,0.025) 
quantile(deltas,0.975)
sum(deltas<0)/1000 #bootstrap P value

#Nomograms and calibration curves
data=df

data$Extent.Of.Resection = as.factor(data$EOR)
data$Adjuvant.Radiotherapy = as.factor(data$Adj_RT)
data$Setting = as.factor(data$Recurrent)
data$Gene.Risk.Score = (data$Chen)

data$Extent.Of.Resection = revalue(data$Extent.Of.Resection, c("GTR"="Gross-total", "STR"="Sub-total"))
data$Adjuvant.Radiotherapy = revalue(data$Adjuvant.Radiotherapy, c("0"="No", "1"="Yes"))
data$Setting = revalue(data$Setting, c("0"="Primary", "1"="Recurrent"))

data$WHO.Grade.2021 = data$WHO.2021
data$WHO.Grade.2016 = data$Grade

data = data[,c('Age','LFFR','LF','OS','VS','Extent.Of.Resection','Adjuvant.Radiotherapy','Setting','Gene.Risk.Score','WHO.Grade.2021','WHO.Grade.2016')]

dd = datadist(data)
options(datadist='dd')
units(data$LFFR)<- 'Years'

f2 <- cph(Surv(LFFR,LF)~Gene.Risk.Score + Extent.Of.Resection +  WHO.Grade.2021 + Setting,data=data,x=TRUE,y=TRUE,surv=TRUE)
t2=Survival(f2)
nom=nomogram(f2,fun=list(function(x)t2(5,x)),funlabel = '5-year LFFR')
plot(nom)

f2 <- cph(Surv(OS,VS)~Age + Gene.Risk.Score + Extent.Of.Resection + WHO.Grade.2021 + Setting,data=data,x=TRUE,y=TRUE,surv=TRUE)
t2=Survival(f2)
nom=nomogram(f2,fun=list(function(x)t2(5,x)),funlabel = '5-year OS')
plot(nom)

f2 <- cph(Surv(LFFR,LF)~Gene.Risk.Score + Extent.Of.Resection +  WHO.Grade.2016 + Setting,data=data,x=TRUE,y=TRUE,surv=TRUE)
t2=Survival(f2)
nom=nomogram(f2,fun=list(function(x)t2(5,x)),funlabel = '5-year LFFR')
plot(nom)

f2 <- cph(Surv(OS,VS)~Age + Gene.Risk.Score + Extent.Of.Resection + WHO.Grade.2016 + Setting,data=data,x=TRUE,y=TRUE,surv=TRUE)
t2=Survival(f2)
nom=nomogram(f2,fun=list(function(x)t2(5,x)),funlabel = '5-year OS')
plot(nom)

# C-index and calibration plots for methylation model
#Optimism corrected c-index
mod1=as.formula(Surv(LFFR,LF)~Gene.Risk.Score+WHO.Grade.2021+Extent.Of.Resection+Setting)
f <- cph(mod1,data=data,x=TRUE,y=TRUE)
calib = rms::calibrate(cph(mod1,data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=5),u=5,conf.int=TRUE,B=1000,cmethod='KM',m=75)
plot(calib,xlab='Predicted 5-year LFFR',ylab='Observed 5-year LFFR',xlim=c(0,1))

# C-index and calibration plots for methylation model
#Optimism corrected c-index
mod1=as.formula(Surv(LFFR,LF)~Gene.Risk.Score+WHO.Grade.2016+Extent.Of.Resection+Setting)
f <- cph(mod1,data=data,x=TRUE,y=TRUE)
calib = rms::calibrate(cph(mod1,data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=5),u=5,conf.int=TRUE,B=1000,cmethod='KM',m=150)
plot(calib,xlab='Predicted 5-year LFFR',ylab='Observed 5-year LFFR',xlim=c(0,1))

# C-index and calibration plots for methylation model
#Optimism corrected c-index
mod1=as.formula(Surv(OS,VS)~Age+Gene.Risk.Score+WHO.Grade.2021+Extent.Of.Resection+Setting)
f <- cph(mod1,data=data,x=TRUE,y=TRUE)
calib = rms::calibrate(cph(mod1,data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=5),u=5,conf.int=TRUE,B=1000,cmethod='KM',m=75)
plot(calib,xlab='Predicted 5-year OS',ylab='Observed 5-year OS',xlim=c(0,1))

# C-index and calibration plots for methylation model
#Optimism corrected c-index
mod1=as.formula(Surv(OS,VS)~Age+Gene.Risk.Score+WHO.Grade.2016+Extent.Of.Resection+Setting)
f <- cph(mod1,data=data,x=TRUE,y=TRUE)
calib = rms::calibrate(cph(mod1,data=data,x=TRUE,y=TRUE,surv=TRUE,time.inc=5),u=5,conf.int=TRUE,B=1000,cmethod='KM',m=150)
plot(calib,xlab='Predicted 5-year OS',ylab='Observed 5-year OS',xlim=c(0,1))

#Propensity Matching
tmp4 = rbind(cbind(subset(data,Site!='UCSF')$Adj_RT,subset(data,Site!='UCSF')$LFFR,subset(data,Site!='UCSF')$LF,subset(data,Site!='UCSF')$Chen,subset(data,Site!='UCSF')$Grade,subset(data,Site!='UCSF')$EOR,subset(data,Site!='UCSF')$ID,subset(data,Site!='UCSF')$Chen.Groups))
tmp4 = tmp4[!is.na(tmp4[,4]),]
colnames(tmp4) = c('RT','LFFR','LF','RiskScore','Grade','EOR','ID','RiskGroup')
tmp5 = as.data.frame(tmp4[complete.cases(as.data.frame(tmp4)),])
tmp5$RiskScore = as.numeric(as.character(tmp5$RiskScore))
tmp5$LFFR = as.numeric(as.character(tmp5$LFFR))
tmp5$LF = as.numeric(as.character(tmp5$LF))

tmp5.low = subset(tmp5,!((RiskGroup=='Intermediate'&EOR=='STR')|RiskGroup=='High'))
tmp5.high = subset(tmp5,(RiskGroup=='Intermediate'&EOR=='STR')|RiskGroup=='High')

#Unfavorable strata match
modmatch = matchit(RT~RiskScore+Grade+EOR,method='nearest',caliper=0.2,ratio = 3,data=tmp5.high,replace = FALSE)
summary(modmatch)
matcheddata = match.data(modmatch)
ggsurvplot(survfit(Surv(LFFR,LF)~RT,data=matcheddata),pval=TRUE,risk.table = TRUE,xlim=c(0,10))

#Favorable strata match
modmatch = matchit(RT~RiskScore+Grade+EOR,method='nearest',caliper=0.2,ratio = 3,data=tmp5.low,replace = FALSE)
summary(modmatch)
matcheddata = match.data(modmatch)
ggsurvplot(survfit(Surv(LFFR,LF)~RT,data=matcheddata),pval=TRUE,risk.table = TRUE,xlim=c(0,10))

