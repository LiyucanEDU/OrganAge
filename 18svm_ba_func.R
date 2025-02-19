library(e1071)
library(caret)
library(randomForest)
svm_ba_func<-function(dat,organ_statname,organ_name,costmax=5,PCA=F,RFE=T,healthyid=F){
  #organ_statname<-statname_list[[8]];organname=names(statname_list)[8]
  organ_variable<-dat[,c('age',organ_statname)]
  #有完全数据的样本ID
  ID_complete<-which(cal_individua_NA_rate(organ_variable)==0 )
  organ_variable<-organ_variable[ID_complete,]
  target <- organ_variable[,'age']
  organ_variable[,-1]<-apply(organ_variable[,-1],2,scale)
  data <- organ_variable # 数据已经准备好

  if(length(which(isFALSE(healthyid)))==0){

    organ_variable_healthy<-dat[,c('ID','age',organ_statname)]
    organ_variable_healthy<-organ_variable_healthy %>% filter(ID %in% healthyid)
    #有完全数据的样本ID
    ID_complete_healthy<-which(cal_individua_NA_rate(organ_variable_healthy)==0 )
    organ_variable_healthy<-organ_variable_healthy[ID_complete_healthy,]
    target_healthy <- organ_variable_healthy[,'age']
    organ_variable_healthy[,-c(1:2)]<-apply(organ_variable_healthy[,-c(1:2)],2,scale)
    data_healthy <- organ_variable_healthy[,-1] # 数据已经准备好
  }

  if(PCA & ncol(data)>50){#特征大于50个
    datX<-organ_variable[,-1]
    datX_pca<-prcomp(datX)
    t<-datX_pca %>% summary();t

    pcs_id<-which(t$importance[3,]>0.95)[1]
    data <- data.frame(datX_pca$x[,1:pcs_id])
    #print(pcs_id)
    data$age <-organ_variable$age
  }

  optVariables<-organ_statname
if(RFE&length(optVariables)>5){
  # 设置训练控制
  ctrl <- rfeControl(functions=rfFuncs, method="cv", number=10)

  svmFuncs <- rfFuncs

  # 使用SVM代替随机森林
  svmFuncs$fit <- function(x, y, first, last, ...) {
    library(e1071)
    out <- svm(x, y, ...)
    out
  }

  # 使用SVM的预测函数
  svmFuncs$pred <- function(object, x)  {
    library(e1071)
    predict(object, x)
  }

  # 使用SVM的等级重要性
  svmFuncs$rank <- function(object, x, y) {
    library(e1071)
    imp <- varImp(object, scale = FALSE)
    imp
  }
  predictors<-data %>% dplyr::select(!age)
  outcome<-data$age

  s<-ifelse(ncol(predictors)>80,80,ncol(predictors))
  result <- rfe(predictors, outcome, sizes=c(2:s), rfeControl=ctrl, functions=svmFuncs)

  optVariables<-result$optVariables
  data<-data[,c('age',optVariables)]
}


  #Tune the SVM model


  if(length(which(isFALSE(healthyid)))==0){
    OptModelsvm=tune(svm, age~., data=data_healthy)#,ranges=list(elsilon=seq(0,0.5,0.1), cost=1:costmax)
    BstModel=OptModelsvm$best.model
  }else{
    OptModelsvm=tune(svm, age~., data=data)#,ranges=list(elsilon=seq(0,0.5,0.1), cost=1:costmax)
    BstModel=OptModelsvm$best.model
  }

  #Predict Y using best model
  PredYBst=predict(BstModel,data)
  svm_best.para<-OptModelsvm$best.parameters

  fit=BstModel
  #feature importance
  w <- t(fit$coefs) %*% fit$SV                 # weight vectors
  w <- apply(w, 2, function(v){sqrt(sum(v^2))})  # weight
  w <- sort(w, decreasing = T)
  #print(w)
  #install.packages('rminer')
  #library(rminer)
  #class(fit)
  #Importance(fit,data)#

  organ_res<-data.frame(dat[ID_complete,c('ID','sex','age')])
  organ_res$BA<-PredYBst
  organ_res$BAgaplm<-lm(BA~age,data=organ_res) %>% resid
  organ_res$BAgap<-with(data=organ_res,(BA-age))
  organ_res$BAgap_adjust<-lm(BAgap~age,data=organ_res) %>% resid
  organ_age_cor=correlation(organ_res,select = c('age','BA'))
  gap_range=range(organ_res$BAgap)

  colnames(organ_res)[4:7]<-c(paste(organ_name,'_BA',sep=''),paste(organ_name,'_BA_Gap_lm',sep=''),
                              paste(organ_name,'_BA_Gap',sep=''),paste(organ_name,'_adj',sep=''))

  checkcor<-data.frame(BA_cor=t(cor(PredYBst,data)),CA_cor=t(cor(data$age,data)))
  checkcor<-checkcor[order(checkcor$BA_cor,decreasing = T),]
  checkcor$BAcor_multiply_CAcor<-with(checkcor,ifelse(BA_cor*CA_cor>0,T,F))


  organ_res_name<-c(paste(organ_name,'_res',sep=''))
  final_list<-list(organ_res=organ_res,checkcor=checkcor,
                   organ_age_cor=organ_age_cor,
                   gap_range=gap_range,svm_best.para=svm_best.para,optVariables=optVariables)
  names(final_list)[1]<-organ_res_name
  return(final_list)
}
