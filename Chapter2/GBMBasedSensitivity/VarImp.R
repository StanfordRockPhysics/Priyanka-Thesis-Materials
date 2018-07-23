#setwd("C:/Work/Post Berkeley/Priyanka Thesis GBM Tuning")

library(R.matlab)
library(beepr)
library(gbm)
library(caret)

# reading in data, for data format reference the file HS_dfp.mat has been included in folder
dat = readMat('HS_dfp.mat')
n = round(2*(dim(dat$c33)[1])/3)
m = dim(dat$c33)[1]

train.data = data.frame(dfp = dat$dfp[1:n],c33 = dat$c33[1:n],c44 = dat$c44[1:n],eps = dat$eps[1:n],gam = dat$gam[1:n],del = dat$del[1:n])
test.data = data.frame(dfp = dat$dfp[(n+1):m],c33 = dat$c33[(n+1):m],c44 = dat$c44[(n+1):m],eps = dat$eps[(n+1):m],gam = dat$gam[(n+1):m],del = dat$del[(n+1):m])


# Sample code for using Caret for GBM parameter tuning, can take a while to run, we used results of the tuning to run GBM in the following section
id<-c(4,5)
s<-c(0.05,0.1,0.3)
trees<-c(5000,10000,15000)

h.params<-expand.grid(id,s,trees)
caret_Grid<-expand.grid(n.trees=(1:3)*5000,interaction.depth=c(3,4,5),shrinkage=c(0.001,0.10,0.05),n.minobsinnode=10)

gbm.caret<-train(dfp~.,data=train.data,method="gbm",train.fraction=0.90,metric="RMSE",trControl=trainControl(method="cv",number=10),verbose=FALSE,tuneGrid=caret_Grid)

print(gbm.caret)


# Running GBM training
set.seed(444)
gbm.train_p<-gbm(dfp~.,data=train.data,train.fraction=0.90,interaction.depth=5,shrinkage=.1,n.trees=5000,bag.fraction=0.5,cv.folds=5,verbose=T)
summary(gbm.train_p)

# Predictions using trained GBM
pred.gbm_p <<- predict(gbm.train_p,newdata = test.data,n.trees = 5000)

# Plots and errors
plot(test.data$dfp,pred.gbm_p,col = 'gray83')
abline(0,1, col = 'blue')
squared.errors<-sqrt(sum((test.data$dfp-pred.gbm_p)^2)/sum(test.data$dfp^2))