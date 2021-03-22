set.seed(1000)
memory.limit(9999999999)

load("C:/Users/qiangli/Desktop/Drug Response Project/DR_Pred/trainFrame.RData")

library(caret)
library(pRRophetic)
library(ridge)
library(MASS)
library(dplyr)
library(caret)
library(glmnet)
