library(e1071)
library(dplyr)
##数据读入
setwd("/Users/xuanyuanqiao/Desktop/homework4/homework4-data")
xtraindata=read.csv("xtrain.csv",header=TRUE)
xtrain=xtraindata[,-1]
ytraindata=read.csv("ytrain.csv",header=TRUE)
ytrain=as.factor(ytraindata[,-1])
##数据输入
dat=data.frame(x=xtrain,y=ytrain)
##构建模型
out=svm(ytrain~., data=dat, kernel="linear", cost =10)
summary(out)

##用训练集查看预测情况（应该是100%）
train_table<-table(out$fitted,ytrain)
train_table

accuracy1<-sum(diag(train_table))/sum(train_table)
accuracy1

##用该模型预测测试集并查看预测情况
##读数据
xtest=read.csv("xtest.csv",header=TRUE)
xtest=xtest[,-1]
ytest=read.csv("ytest.csv",header=TRUE)
ytest=as.factor((ytest)[,-1])

data.te=data.frame(x=xtest,y=ytest)
##进行预测
model_predict<-predict(out,newdata=data.te)
train_table2<-table(model_predict, data.te$y)
train_table2
##计算准确率
accuracy2<-sum(diag(train_table2))/sum(train_table2)
accuracy2
