library(e1071)
library(caret)
library(randomForest)
svm_ba_func<-function(dat,organ_statname,organ_name,costmax=5){
  #organ_statname<-statname_list[[8]];organname=names(statname_list)[8]
  organ_variable<-dat[,c('age',organ_statname)]
  #有完全数据的样本ID
  ID_complete<-which(cal_individua_NA_rate(organ_variable)==0 )
  organ_variable<-organ_variable[ID_complete,]
  target <- organ_variable[,'age']
  organ_variable[,-1]<-apply(organ_variable[,-1],2,scale)
  data <- organ_variable # 数据已经准备好
  optVariables<-organ_statname

  #Tune the SVM model
  OptModelsvm=tune(svm, age~., data=data,ranges=list(elsilon=seq(0,0.5,0.1), cost=1:costmax))#
  BstModel=OptModelsvm$best.model

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


