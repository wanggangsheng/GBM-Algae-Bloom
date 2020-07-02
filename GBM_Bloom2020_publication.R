
## source("GBM_Bloom2020_publication.R")

#Clean up Environment
rm(list = ls())
setwd("/Users/wgs/ownCloud/Rcode/Algae")
# ensure the results are repeatable
set.seed(2019)

# load the library
library(mlbench)
library(dplyr)
library(lavaan)
library(corrplot)
library(caret)
library(rpart)
library(gbm) ## relative.influence
library(rattle)
library(ggplot2)
library(ggsignif)
library(ggpubr) ##stat_regline_equation
library(tidyr)     ## gather: combine multiple columns into a single column 
library(agricolae)  ##LSD1, tapply.stat
library(psych) ## describe.by
library(matrixStats) # colQuantiles



path <- getwd()

group.colors2 <- c("yellowgreen","orange")
group.colors  <- c("cyan","orange","blue","red")

file_par0 <- "Algae_input_data.txt"
##------------------------------------------------
# path_io <- paste0(path, "/", dir_io)
path_io <- path
file_par <- paste0(path,"/", file_par0)
prefix_io <- "fig"
path_fig <- paste0(path_io, "/",prefix_io)
ifelse(!dir.exists(file.path(path_io, prefix_io)), dir.create(file.path(path_io, prefix_io)), FALSE)

par0 <- read.table(file_par, sep="\t", header=TRUE)
par1 <- par0[par0$Month < 9,]
write.csv(par1,paste0(path_io, "/Rdata_GBM.csv"))
# par1 <- par0[par0$Month < 9 & par0$Year > 2003,]
par1 <- par1[,4:ncol(par0)]
par1$Bloom <- ifelse(par1$Algae >= 1, "presence","absence")

# mean1 <- aggregate(par1,)
##---------------------------------------------------------

describe.by(par1,par1$Bloom)

prob <- c(0, 2.5, 5, 25, 50, 75, 95, 97.5, 100)
group_names <- unique(par1$Bloom) 
for(i in 1:length(group_names)) {
  sub_par1 <- par1[par1$Bloom == group_names[i],]
  qt <- colQuantiles(data.matrix(sub_par1),  probs = prob/100)
  print(group_names[i])
  print(qt)
}

par01 <- par0[,which(colnames(par0)%in%c("HRWL","HRV","HRQ","YRWL","YRV","YRQ"))]
library("PerformanceAnalytics")
chart.Correlation(par01, histogram=TRUE, pch=19)

par4 <- par3[,which(!colnames(par3)%in%c("HRQ","YRQ"))]

## with predictors from previous 10-day
par5 <- par4[,which(substr(colnames(par4),1,1) == "p")]

## with predictors from current 10-day
# par5 <- par4[,which(!substr(colnames(par4),1,1) == "p")]

par5$Bloom <- par4$Bloom
names(par5)
##---------------------------------------------------------

## Booted Regression Tree: two classes
## with predictors from current and previous 10-day
# par <- par4

## with predictors from current or previous 10-day
par <-par5
npar <- ncol(par) - 1
par_name <- names(par)[1:npar]
group <- names(par)[npar + 1]
par_log <- vector(mode="integer",length=npar)
par_log[1:npar] <- 0

xlab1 <- "Algal Bloom"
par_label <- par_name
color_use <- group.colors2

par$Bloom <- as.factor(par$Bloom)  ## otherwise Error: All group levels must be finite


##------------------------------------------------------------------------------------------------

method1 <- 'gbm' #this is where model type is specified, could be rf, gbm, rpart
method2 <- 'rpart'

outcomeName <- "Bloom"

par[,outcomeName] <- as.factor(par[,outcomeName])
predictorNames <- names(par)[names(par) != outcomeName]
predictorNames_sort <- predictorNames[order(predictorNames)]
npredictor <- length(predictorNames)
objControl    <- trainControl(method='cv', number=10 , returnResamp='none', 
                              summaryFunction = twoClassSummary, classProbs = TRUE) #10 fold cross-validation

desiredOutput   <- "raw"

summary(par[,outcomeName])


##----------------------------------------------------------------------------
## model with all data
# ensure the results are repeatable
set.seed(2019)

nLoop_max <- 10000
var_imp_nrun <- matrix(nrow=nLoop_max,ncol=npredictor+2) ## n predictors + 2 criteia (Accuracy & Kappa)
var_imp_nrun <- as.data.frame(var_imp_nrun)
colnames(var_imp_nrun) <- c(predictorNames_sort,"Accuracy", "Kappa")
nrun <- 0
nLoop <- 0
Kappa_train0 <- 0
Kappa_train <- Kappa_train0
model_list <- vector("list", length=nLoop_max)
# while (nLoop <= nLoop_max & (Kappa_train0 < 1.0))
while (nLoop < nLoop_max)
{
  nLoop <- nLoop + 1
  model <- train(par[,predictorNames], par[,outcomeName],
                      method=method1, 
                      trControl=objControl,
                      verbose = FALSE,
                      metric = "ROC") ## ROC or Kappa


  # print("The most important predictors")
  variable_importance <- varImp(model,scale=TRUE)
  # print(varImp(model,scale=T))
  # plot(variable_importance)
  # filterVarImp(par[,predictorNames], as.factor(par[,outcomeName]))

  predictions <- predict(object=model, par[,predictorNames], type='raw')
  # print("Results on the training data")
  obs <- as.factor(ifelse(par[,outcomeName]=="presence",1,0))
  sim <- as.factor(ifelse(predictions=="presence",1,0))
  levels(sim) <- c(0,1)
  confusion_Matrix <- confusionMatrix(obs, sim)
  Kappa_train <- confusion_Matrix[3]$overall[2]
  Accuracy_train <- confusion_Matrix[3]$overall[1]
  print(paste0("Loop = ",nLoop,"; Kappa = ", Kappa_train))

  vimp1 <- variable_importance$importance
  var_imp_nrun[nLoop,1:npredictor] <- vimp1[predictorNames_sort,]
  var_imp_nrun[nLoop,npredictor + 1] <- Accuracy_train
  var_imp_nrun[nLoop,npredictor + 2] <- Kappa_train

    if (Kappa_train == 1.0) {
      nrun <- nrun  + 1
      results <- list(model=model,varImp=vimp1,Accuracy=Accuracy_train,Kappa=Kappa_train)
      model_list[[nrun]] <- results

    }


} ## end while

# model_list_best <- model_list[1:nrun]
##-------------------------------------------------------
## nrun = 1950
save(model_list, file = paste0(path_io, "/model_list.Rdata"))

Algae_prob <- matrix(nrow=nrow(par),ncol=nrun) 
for (i in 1:nrun){
  predictions_prob <- predict(object=model_list[[i]]$model, par[,predictorNames], type='prob')
  Algae_prob[,i] <- predictions_prob[,2]
}

df1 <- cbind(par, Algae_prob)
fn1 <- paste0(path_io, "/Algae_prob.txt")
write.table(df1,fn1,sep="\t",row.names=FALSE,col.names = TRUE)


var_imp_nrun2 <- as.data.frame(var_imp_nrun[var_imp_nrun$Kappa==1,1:npredictor])
var_imp_mean <- colMeans(var_imp_nrun2)
var_imp_mean

st=format(Sys.time(), "%y%m%d-%H%M%S")
fn2 <- paste0(path_out, "/varImp_nrun_",st,".txt")
write.table(round(var_imp_nrun,2),fn2,sep="\t",row.names=FALSE)



