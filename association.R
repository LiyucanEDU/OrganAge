
#1.加载R包====
library(compareGroups)
library(readxl)
library(dplyr)
#2.读取ba数据====
bagaplm<-svm_res_0120[[1]][[4]]
bagap<-svm_res_0120[[1]][[6]]

organ_complete<-c('Cardiovascular','Bone','Gut','Kidney','Metabolism','Brain','Cognition','Body')#','Body.Brain'
organ_simple<-c('Cardiac','Bone','Gut','Kidney','Metab.','Brain','Cognition','Body')#,'Body.Brain'


ba_descr<-descrTable(~.,data = ba,show.n = T,show.all = T)
bagap_descr<-descrTable(sex~.,data = bagap,show.n = T,show.all = T)
bagaplm_descr<-descrTable(sex~.,data = bagaplm,show.n = T,show.all = T)



ba_predicted<-bagap %>% mutate(Cardiovascular=`Cardiac`+age,
                            Bone=`Bone`+age,
                            Gut=`Gut`+age,
                            Brain=`Brain`+age,
                            Kidney=`Kidney`+age,
                            Metabolism=`Metab.`+age,
                            Cognition=`Cognition`+age,
                            Body=`Body`+age) %>%
  dplyr::select(all_of(c('ID','sex','age',organ_complete)))

ba_predicted<-bagap %>% mutate(Cardiac=`Cardiac`+age,
                               Bone=`Bone`+age,
                               Gut=`Gut`+age,
                               Brain=`Brain`+age,
                               Kidney=`Kidney`+age,
                               Metab.=`Metab.`+age,
                               Cognition=`Cognition`+age,
                               Body=`Body`+age) %>%
  dplyr::select(all_of(c('ID','sex','age',organ_simple)))


correlation::correlation(ba,select = 'age',select2 =organ_complete,
                         p_adjust = 'none',include_factors = T)

ba_predicted_bodysepe<-full_join(ba_predicted[c('ID','sex','age','Cardiovascular',
                                       'Bone','Gut','Kidney','Metabolism',
                                       'Brain','Cognition')],svm_res_1111_body_BA)

save(ba_predicted_bodysepe,file='./06newoutput/09SvmRes/BodySexSepa_ba.Rdata')

colnames(ba_predicted)

##画BA CA线性相关点图====
source("25BA-CA_BioAge.R")

filename='./06newoutput/14BACAdot/baca_1113.pdf'
filename='./06newoutput/14BACAdot/f2_baca_0311.pdf'
plot_BACA(ba_predicted,organ_order = organ_simple,
          nrow =1,width = 13, height =4,legend.position='bottom')
ggsave(filename = filename, plot=p,
       device='pdf',width = width, height = height, dpi = 600, units = 'in')


ba_man<-ba_predicted %>% filter(sex==1)
ba_woman<-ba_predicted %>% filter(sex==2)


filename='./06newoutput/14BACAdot/f2_bacaman_0305.pdf'
plot_BACA(ba_man,filename,organ_order = organ_simple,
          nrow =1,width = 13, height =3,legend.position='bottom')
filename='./06newoutput/14BACAdot/f2_bacawoman_0305.pdf'
plot_BACA(ba_woman,filename,organ_order = organ_simple,
          nrow =1,width = 13, height =3,legend.position='bottom')


#3.ba间的偏相关性====
#source("D:/work/01衰老/02标志物/20230417/22interplay_func.R")
load("D:/work/01衰老/02标志物/20230417/03output/svm_res_0602.Rdata")
##3.1 gap之间的偏相关性


bagaplm_pcor<-ba_gap_pcor_func(bagaplm,organBA_name=organ_simple)
bagaplm_man_pcor<-ba_gap_pcor_func(bagaplm_man,organBA_name=organ_simple)
bagaplm_woman_pcor<-ba_gap_pcor_func(bagaplm_woman,organBA_name=organ_simple)

#source("D:/work/01衰老/02标志物/20230417/24corheatmap_Heatmap.R")
filename='./06newoutput/14BACAdot/gaplmheat0120.pdf'
cor_plot_withPvalue_heatmap(bagaplm_pcor[,c(1,2,3,9)],filename=filename,
                            organ_order =organ_simple)
filename='./06newoutput/14BACAdot/gaplm_manheat0120.pdf'
cor_plot_withPvalue_heatmap(bagaplm_man_pcor[,c(1,2,3,9)],filename=filename,
                            organ_order =organ_simple)
filename='./06newoutput/14BACAdot/gaplm_womanheat0120.pdf'
cor_plot_withPvalue_heatmap(bagaplm_woman_pcor[,c(1,2,3,9)],filename=filename,
                            organ_order =organ_simple)


##3.2年龄之间的相关热图及网络====
balm<-bagaplm %>% mutate(Cardiac=`Cardiac`+age,
                            Bone=`Bone`+age,
                            Gut=`Gut`+age,
                            Brain=`Brain`+age,
                            Kidney=`Kidney`+age,
                            Metab.=`Metab.`+age,
                            Cognition=`Cognition`+age,
                            Body=`Body`+age) %>%
  dplyr::select(all_of(c('ID','sex','age',organ_simple)))
balm_man<-balm %>% filter(sex==1)
balm_woman<-balm %>% filter(sex==2)

balm_descr<-descrTable(sex~.,data = balm,show.n = T,show.all = T,digits = 3)
bagaplm_descr<-descrTable(sex~.,data = bagaplm,show.n = T,show.all = T,digits = 3)
export2csv(balm_descr, file='./06newoutput/16table/f2_balm0120_descr.csv')
export2csv(bagaplm_descr, file='./06newoutput/16table/f2_bagaplm0120_descr.csv')


ba_pcor<-ba_pcor_func(balm,organ_name=organ_simple)#
ba_woman_pcor<-ba_pcor_func(balm_woman,organ_name=organ_simple)
ba_man_pcor<-ba_pcor_func(balm_man,organ_name=organ_simple)

stand_ba_pcor<-function(ba_pcor){
ba_pcor<-ba_pcor %>% dplyr::select(v1,v2,n,estimate,p.value,p.fdr)
colnames(ba_pcor)<-c('organ1','organ2','n','r','p.value','p.fdr')
ba_pcor$r<-sprintf("%.3f", ba_pcor$r)
ba_pcor$p.value<-format(signif(ba_pcor$p.value, digits = 4), scientific = TRUE)
ba_pcor$p.fdr<-format(signif(ba_pcor$p.fdr, digits = 4), scientific = TRUE)
return(ba_pcor)}

stand_ba_pcor<-stand_ba_pcor(ba_pcor)
stand_ba_woman_pcor<-stand_ba_pcor(ba_woman_pcor)
stand_ba_man_pcor<-stand_ba_pcor(ba_man_pcor)

stand_bagaplm_pcor<-stand_ba_pcor(bagaplm_pcor)
stand_bagaplm_woman_pcor<-stand_ba_pcor(bagaplm_woman_pcor)
stand_bagaplm_man_pcor<-stand_ba_pcor(bagaplm_man_pcor)


write.csv(stand_ba_pcor,file = './06newoutput/16table/f2_ba_pcor_0120.csv',row.names = F)
write.csv(stand_ba_woman_pcor,file = './06newoutput/16table/f2_ba_pcor_woman_0120.csv',row.names = F)
write.csv(stand_ba_man_pcor,file = './06newoutput/16table/f2_ba_pcor_man_0120.csv',row.names = F)
write.csv(stand_bagaplm_pcor,file = './06newoutput/16table/f2_bagaplm_pcor_0120.csv',row.names = F)
write.csv(stand_bagaplm_woman_pcor,file = './06newoutput/16table/f2_bagaplm_woman_pcor_0120.csv',row.names = F)
write.csv(stand_bagaplm_man_pcor,file = './06newoutput/16table/f2_bagaplm_man_pcor_0120.csv',row.names = F)



filename='./03output/07PcorHeatmap/organbaheat0712.pdf'
filename='./06newoutput/14BACAdot/baheat0120.pdf'
cor_plot_withPvalue_heatmap(ba_pcor[,c(1,2,3,9)],filename=filename,
                            organ_order = organ_simple)
filename='./06newoutput/14BACAdot/ba_woman_heat0120.pdf'
cor_plot_withPvalue_heatmap(ba_woman_pcor[,c(1,2,3,9)],filename=filename,
                            organ_order = organ_simple)
filename='./06newoutput/14BACAdot/ba_man_heat0120.pdf'
cor_plot_withPvalue_heatmap(ba_man_pcor[,c(1,2,3,9)],filename=filename,
                            organ_order = organ_simple)




##由ba相关性构建网络====
filename='./06newoutput/14BACAdot/banet0120.pdf'
ba_cor_net_func(ba_pcor,filename = filename)
filename='./06newoutput/14BACAdot/banet_woman_1128.pdf'
ba_cor_net_func(ba_woman_pcor,filename = filename)
filename='./06newoutput/14BACAdot/baGaplmNet0310.pdf'
ba_cor_net_func(bagaplm_pcor,filename = filename)



#散点相关性图
##bioage包
axis_type<-rep('float',length(organ_complete))
names(axis_type)<-organ_complete
plot_baa(bagap, agevar=organ_simple,
         label=organ_complete,
         axis_type=axis_type)
library(psych)
pairs.panels(bagap[,-c(1:3)],
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)



#5.ba和生活方式====
#source("D:/work/01衰老/02标志物/20230417/22interplay_func.R")
##5.1整理好的生活方式====
load(file='./03output/08lifestyle/lifestyle.Rdata')
load(file='./03output/08lifestyle/lifestyle_use.Rdata')
##5.2植物性饮食指数====
pdi<-read_excel("01input_data/植物性饮食指数.xlsx")
pdi$hPDI2<-as.numeric(as.factor(pdi$hPDI2))
pdi$PDI2<-as.numeric(as.factor(pdi$PDI2))
pdi$uPDI2<-as.numeric(as.factor(pdi$uPDI2))
dput(colnames(pdi))#c("ID", "hPDI1", "hPDI2", "PDI1", "PDI2", "uPDI1", "uPDI2")
colnames(pdi)<-c("ID", "hPDI", "hPDI_cate", "PDI", "PDI_cate", "uPDI", "uPDI_cate")
pdi_test<-full_join(bagap,pdi)
write.csv(pdi,file = '02derived_data/PDI.csv')

#连续型====
correlation(pdi_test,select = organ_simple,select2 = c('hPDI1','PDI1','uPDI1'),
            p_adjust = 'none')##相反的====
correlation(pdi_test,select = organ_simple,select2 = c('uPDI2'),#'hPDI2','PDI2',
            p_adjust = 'none')
correlation(pdi_test,select = organ_simple,select2 = c('hPDI2','PDI2','uPDI2'),#,
            p_adjust = 'fdr')

lifestyle<-left_join(lifestyle_dat_correction,pdi)
colnames(lifestyle_use)
lifestyle_use<-rbind(lifestyle_use,
                     data.frame(name=c('hPDI1','PDI1','uPDI1','hPDI2','PDI2','uPDI2','Cognition'),
                                              label=c('healthy','healthy','unhealthy',
                                                      'healthy','healthy','unhealthy','healthy'),
                                reproductive_label=c(F,F,F,F,F,F,F)))
##5.3健康状态自评互评====
healthstatus<-dat2023_0905 %>% dplyr::select(ID,healthstatus_comparation,healthstatus_self)
healthstatus$healthstatus_comparation[which(healthstatus$healthstatus_comparation==4)]<-NA
healthstatus$healthstatus_comparation<-multireplace(healthstatus$healthstatus_comparation,1:3,3:1)
healthstatus$healthstatus_self<-multireplace(healthstatus$healthstatus_self,1:4,4:1)
healthstatus_test<-full_join(bagap,healthstatus)
colnames(healthstatus_test)
correlation(healthstatus_test,select =organ_simple,select2 = c("healthstatus_comparation","healthstatus_self"),
            p_adjust = 'none')

colnames(ba_scale1)<-multireplace(colnames(ba_scale1),organ_simple,organ_complete)
organBA_name<-organ_complete
##5.4睡眠状况====
dat2023_0905<-read.csv('./02derived_data/dat2023_0905.csv')
lifestyle_dat_correction$notclear_duringdays %>% head#0是清醒，1是不清醒
dat2023_0905$notclear_duringdays %>% head

lifestyle_dat_correction$sleeplate30_not %>% head#0是没入睡困难有，1是有
dat2023_0905$sleeplate30_not %>% head

sleep<-dat2023_0905 %>% dplyr::select(ID,notclear_duringdays,sleeplate30_not,
                                 sleep_time_hrs,wakeearly_not,snore_not)
##
sleep$sleep_time_hrs_cate<-ifelse(sleep$sleep_time_hrs < 6, 1, ifelse(sleep$sleep_time_hrs <= 9, 2, 1))
sleep$snore_not_cate<-ifelse(sleep$snore_not<3,1,2)#3是不知道或者否，
complete.cases(sleep) %>% which %>% length()
sleep<-sleep[complete.cases(sleep),]
sleep$sleep_score<-sleep %>% dplyr::select(notclear_duringdays,sleeplate30_not,
                           sleep_time_hrs_cate,wakeearly_not,snore_not_cate) %>% rowSums()
table(sleep$sleep_score)
sleep$sleep_cate<-ifelse(sleep$sleep_score<7,1,ifelse(sleep$sleep_score <= 9, 2, 3))
sleep_test<-full_join(bagap,sleep)
correlation(sleep_test,select = organ_simple,select2 = c('notclear_duringdays','sleeplate30_not',
                                                         'sleep_time_hrs',
                                                         'wakeearly_not','snore_not'),
            p_adjust = 'fdr')
correlation(sleep_test,select = organ_simple,select2 = c('score'),
            p_adjust = 'fdr')
correlation(sleep_test,select = organ_simple,select2 = c('score'),
            p_adjust = 'none')
correlation(sleep_test,select = organ_simple,select2 = c('score_cate'),
            p_adjust = 'none')

write.csv(sleep,file = './02derived_data/sleep0916.csv',row.names = F)
ba_lifestyle<-left_join(ba_scale1,lifestyle_dat_correction)
ba_lifestyle$exercise_frequency %>% table
##5.5筛选生活方式====
life_new<-c('smoke', 'drink', 'tea','exercise_frequency','sleep_cate',
'Cognition', 'Identification','TUG', 'Tinetti',
 "FreshFruit",'Dairy','Dessert','Staple',#"FreshVegetable",
"gum_bleeding",  "stool_scale",
#"bad_breath", "abdominal_distension","fart",  "diarrhea", "tenesmus","borborygmus",
# "upset_intestinal_problems", "odor_fart", "antibiotic_usual",
#"sleep_drughelp", "notclear_duringdays", "sleep_time_hrs", "sleeplate30_not","wakeearly_not",
'menophania_age','menses_regular','menses_interval1','menopause_age')
life_new_label<-multireplace(life_new,lifestyle_use$name,lifestyle_use$label)
life_new_reproductive_label<-multireplace(life_new,lifestyle_use$name,lifestyle_use$reproductive_label)
life_new_reproductive_label<-as.logical(life_new_reproductive_label)
life_new_reproductive_label[c(4:10)]<-FALSE
life_new_name<-c('Smoke','Alcohol','Tea','Exercise Frequency','Sleep score',
                 'Cognition age','Olfactory identification', 'Tinetti','TUG',
                  "Fresh Fruit",'Dairy',"Fresh Vegetable",'Staple',#'Dessert',
                 "Gum bleeding",  "Stool scale",
                # "Bad breath","Abdominal distension","Fart", "Diarrhea", "Tenesmus","Borborygmus",
                # "upset intestinal problems", "Odor fart", "Antibiotic usual",
                #"Sleep drughelp", "Not clear during days", "Sleep time hrs", "Sleeplate30 not","Wakeearly not",
                'Age at menarche ','Menses regular','Menstrual cycle','Age at menopause')


##5.6热图====
organ_complete<-c('Cardiovascular','Bone','Gut','Kidney','Metabolism','Brain','Body')#,'Cognition','Body.Brain'
organ_simple<-c('Cardiac','Bone','Gut','Kidney','Metab.','Brain','Body')
bagap<-svm_res_0915[[4]]
ba_lifestyle<-left_join(bagap,lifestyle_dat_correction)
ba_lifestyle<-left_join(ba_lifestyle,sleep[,c('ID','sleep_cate')])
colnames(ba_lifestyle)<-gsub("\nGap", "", colnames(ba_lifestyle), fixed = TRUE)
colnames(ba_lifestyle)<-multireplace(colnames(ba_lifestyle),organ_simple,organ_complete)
colnames(ba_lifestyle)<-multireplace(colnames(ba_lifestyle),life_new,life_new_name)

ba_lifestyle_1<-ba_lifestyle %>% filter(sex=='1')
ba_lifestyle_2<-ba_lifestyle %>% filter(sex=='2')


filename='./03output/08lifestyle/0722_fdr_max.png'
filename='./06newoutput/11lifestyle/0912.png'
ba_lifestyle_cor<-ba_lifestyle_cor_plot(ba_lifestyle,organBA_name,lifestyle_use$name,
                                        lifestyle_use$label,lifestyle_use$reproductive_label,
                                        filename,filter_p ='p' )#p.fdr
dput(lifestyle_use$name)
setdiff(colnames(ba_lifestyle),lifestyle_use$name)
setdiff(lifestyle_use$name,colnames(ba_lifestyle))

lifestyle_name_select<-setdiff(lifestyle_use$name,c('rhythm','variability','pace',
                                                    'hPDI1','PDI1','uPDI1','hPDI2',
                                                    'PDI2','uPDI2',))
filename='./06newoutput/11lifestyle/0914.png'
filename='./06newoutput/11lifestyle/0914_fdr1.pdf'
dput(organ_complete)
organ_complete %in% colnames(ba_lifestyle)
life_new_name %in% colnames(ba_lifestyle)
ba_lifestyle$`Stool scale` [ba_lifestyle$`Stool scale`==8]<-NA #8是不固定
ba_lifestyle_cor_all<-ba_lifestyle_cor_plot(ba_lifestyle,organBA_name=c("Cardiovascular", "Bone", "Gut", "Kidney", "Metabolism", "Brain",
                                                                        "Body"),
                                            life_new_name,
                                        life_new_label,life_new_reproductive_label,
                                        filename,filter_p ='p.fdr' )#
write.csv(ba_lifestyle_cor_all,file = './06newoutput/11lifestyle/0914_fdr1.csv')
#man
filename='./06newoutput/11lifestyle/0914_fdr_MAN1.pdf'
ba_lifestyle_cor_man<-ba_lifestyle_cor_plot(ba_lifestyle_1,organBA_name=c("Cardiovascular", "Bone", "Gut", "Kidney", "Metabolism", "Brain",
                                                                          "Body"),
                                            life_new_name,
                                        life_new_label,life_new_reproductive_label,
                                        filename,filter_p ='p.fdr' )
write.csv(ba_lifestyle_cor_man,file = './06newoutput/11lifestyle/0914_fdr_MAN1.csv')
#woman
filename='./06newoutput/11lifestyle/0914_fdr_WOMAN1.pdf'
ba_lifestyle_cor_woman<-ba_lifestyle_cor_plot(ba_lifestyle_2,organBA_name=c("Bone", "Gut", "Cardiovascular", "Kidney", "Metabolism", "Brain",
                                                                            "Body"),
                                        life_new_name,
                                        life_new_label,life_new_reproductive_label,
                                        filename,filter_p ='p.fdr',t_mat = F )
write.csv(ba_lifestyle_cor_woman,file = './06newoutput/11lifestyle/0914_fdr_WOMAN1.csv')

##5.7所有生活方式，不矫正====
load(file='./03output/08lifestyle/lifestyle.Rdata')
load(file='./03output/08lifestyle/lifestyle_use.Rdata')
pdi<-read.csv('02derived_data/PDI.csv')
sleep<-read.csv('./02derived_data/sleep0916.csv')
dput(colnames(pdi))
dput(colnames(sleep))
lifestyle_use<-rbind(lifestyle_use,
                     data.frame(name=c("hPDI", "hPDI_cate", "PDI", "PDI_cate", "uPDI", "uPDI_cate",
                                        "snore_not", "sleep_score", "sleep_cate"),
                                label=c('healthy','healthy','unhealthy','healthy','healthy','unhealthy',
                                        'healthy','healthy','healthy'),
                                reproductive_label=c(F,F,F,F,F,F,F,F,F)))
lifestyle<-full_join(lifestyle_dat_correction,sleep[,c('ID','snore_not',"sleep_score", "sleep_cate")])
lifestyle<-left_join(lifestyle,pdi)
dput(colnames(lifestyle))
save(lifestyle_use,file = './03output/08lifestyle/lifestyle_use0920.Rdata')
save(ba_lifestyle,file = './03output/08lifestyle/ba_lifestyle0920.Rdata')
save(lifestyle,file = './03output/08lifestyle/lifestyle0920.Rdata')
save(lifestyle,file = './03output/08lifestyle/lifestyle1111.Rdata')


organ_complete<-c('Cardiovascular','Bone','Gut','Kidney','Metabolism','Brain','Cognition','Body')
organ_simple<-c('Cardiac','Bone','Gut','Kidney','Metab.','Brain','Cognition','Body')
bagap<-svm_res_0915[[4]]
ba_lifestyle<-left_join(svm_res_0121[[4]],lifestyle_dat_correction)
ba_lifestyle<-left_join(ba_lifestyle,sleep)
ba_lifestyle<-left_join(ba_lifestyle,pdi)
colnames(ba_lifestyle)<-gsub("\nGap", "", colnames(ba_lifestyle), fixed = TRUE)
colnames(ba_lifestyle)<-multireplace(colnames(ba_lifestyle),organ_simple,organ_complete)


ba_lifestyle_1<-ba_lifestyle %>% filter(sex=='1')
ba_lifestyle_2<-ba_lifestyle %>% filter(sex=='2')

organ_complete %in% colnames(ba_lifestyle)
lifestyle_use$name %in% colnames(ba_lifestyle)
dput(lifestyle_use$name)

filename='./06newoutput/11lifestyle/0918_all.pdf'
filename='./06newoutput/11lifestyle/0918_all_fdr.pdf'
filename='./06newoutput/11lifestyle/1129.pdf'
ba_lifestyle_cor_all<-ba_lifestyle_pcor_plot(ba_lifestyle,organBA_name=organ_complete,
                                            lifestyle_use$name,
                                            lifestyle_use$label,lifestyle_use$reproductive_label,
                                            filename,filter_p ='p' ,
                                            disease_dat_num = disease_num,
                                            width = 20)#.fdr
#
lifestyle
Identification_cor<-correlation(lifestyle,select = 'Identification',select2 =lifestyle_use$name,
            p_adjust = 'fdr' )
lifestyle
##5.8继承5.6====
load('./03output/08lifestyle/lifestyle_use0920.Rdata')
load('./03output/08lifestyle/lifestyle0920.Rdata')
load('./06newoutput/09SvmRes/svm_res_0121.Rdata')
bagaplm<-svm_res_0121[[1]][[4]]
lifestyle$family_num[lifestyle$family_num==0]<-NA
lifestyle$income_average <- ifelse(is.na(lifestyle$total_income) | is.na(lifestyle$family_num), NA, lifestyle$total_income/lifestyle$family_num)
save(lifestyle,file = './03output/08lifestyle/lifestyle1129.Rdata')
organ_complete<-c('Cardiovascular','Bone','Gut','Kidney','Metabolism','Brain','Cognition','Body')#','Body.Brain'
organ_simple<-c('Cardiac','Bone','Gut','Kidney','Metab.','Brain','Cognition','Body')

life_new<-c( "edu_y", "family_num",'income_average',#"total_income",
             'smoke', 'drink', 'tea','exercise_frequency','sleep_cate',#
            'hPDI','uPDI','PDI',
            'TUG', 'Tinetti','Identification',
            'menophania_age','menses_regular','menses_interval1','menopause_age'
            #"FreshFruit",'Dairy','Dessert','Staple',#"FreshVegetable",
            #"gum_bleeding",  "stool_scale",
            #"bad_breath", "abdominal_distension","fart",  "diarrhea", "tenesmus","borborygmus",
            # "upset_intestinal_problems", "odor_fart", "antibiotic_usual",
            #"sleep_drughelp", "notclear_duringdays", "sleep_time_hrs", "sleeplate30_not","wakeearly_not",
           )
life_new_label<-multireplace(life_new,lifestyle_use$name,lifestyle_use$label)
life_new_reproductive_label<-multireplace(life_new,lifestyle_use$name,lifestyle_use$reproductive_label)

life_new_reproductive_label<-rep(FALSE,length(life_new_reproductive_label))
life_new_reproductive_label[(length(life_new_reproductive_label)-3):length(life_new_reproductive_label)]<-TRUE
life_new_reproductive_label<-as.logical(life_new_reproductive_label)
life_new_name<-c("Education years","No.of family members",'Per capita income',#"Total income",
                 'Smoke','Alcohol','Tea','Exercise frequency','Sleep score',
                 'hPDI','uPDI','PDI',
                  'TUG','Tinetti','Olfactory identification',
                 'Age at menarche ','Menses regular','Menstrual cycle length','Age at menopause'
                 #"Fresh Fruit",'Dairy',"Fresh Vegetable",'Staple',#'Dessert',
                 #"Gum bleeding",  "Stool scale",
                 # "Bad breath","Abdominal distension","Fart", "Diarrhea", "Tenesmus","Borborygmus",
                 # "upset intestinal problems", "Odor fart", "Antibiotic usual",
                 #"Sleep drughelp", "Not clear during days", "Sleep time hrs", "Sleeplate30 not","Wakeearly not",
                 )

ba_lifestyle<-left_join(svm_res_0121[[1]][[4]],lifestyle)
#ba_lifestyle<-left_join(ba_lifestyle,sleep[,c('ID','sleep_cate')])
colnames(ba_lifestyle)<-gsub("\nGap", "", colnames(ba_lifestyle), fixed = TRUE)
colnames(ba_lifestyle)<-multireplace(colnames(ba_lifestyle),organ_simple,organ_complete)
colnames(ba_lifestyle)<-multireplace(colnames(ba_lifestyle),organ_complete,organ_simple)
colnames(ba_lifestyle)<-multireplace(colnames(ba_lifestyle),life_new,life_new_name)
dput(colnames(ba_lifestyle))
ba_lifestyle$sex
ba_lifestyle_1<-ba_lifestyle %>% filter(sex=='1')
ba_lifestyle_2<-ba_lifestyle %>% filter(sex=='2')


filename='./06newoutput/11lifestyle/1129_gaplm_fdr.pdf'
filename='./06newoutput/11lifestyle/1129_gaplm_p.pdf'
filename='./06newoutput/11lifestyle/1129_gaplm_fdr.pdf'
filename='./06newoutput/11lifestyle/1130_gaplm_fdr.pdf'
filename='./06newoutput/11lifestyle/1205_gaplm_fdr_long.pdf'
filename='./06newoutput/11lifestyle/0120_gaplm_fdr_long.pdf'
filename='./06newoutput/11lifestyle/0323_gaplm_fdr_long.pdf'
dput(organ_complete)
organ_simple %in% colnames(ba_lifestyle)
life_new_name %in% colnames(ba_lifestyle)

#ba_lifestyle$`Stool scale` [ba_lifestyle$`Stool scale`==8]<-NA #8是不固定
disease_fac<-read.csv('./06newoutput/12violinplot/disease_fac_1129.csv')
disease_num<-disease_fac[,c('ID','MetabolicHealthsum')]
disease_num$MetabolicHealthsum<-as.numeric(disease_num$MetabolicHealthsum)

ba_lifestyle_cor_all<-ba_lifestyle_pcor_plot(ba_lifestyle,
                                            organBA_name=organ_simple,
                                            life_new_name,
                                            life_new_label,life_new_reproductive_label,
                                            filename,filter_p ='p.fdr',t_mat=F,#p
                                           # disease_dat_num = disease_num[,1],
                                            pcor_var = c('age','sex'),
      width=5.5,height=6,#width=7,height=5.5,
      )#,'MetabolicHealthsum'

ba_lifestyle_cor_all_stand<-ba_lifestyle_cor_all %>% dplyr::select(Parameter1,Parameter2,
                                                                   n,r,p,p.fdr)
colnames(ba_lifestyle_cor_all_stand)[1:2]<-c('Agegap','Lifestyle')
ba_lifestyle_cor_all_stand$r<-sprintf("%.3f", ba_lifestyle_cor_all_stand$r)
ba_lifestyle_cor_all_stand$p<- format(signif(ba_lifestyle_cor_all_stand$p, digits = 4), scientific = TRUE)
ba_lifestyle_cor_all_stand$p.fdr<- format(signif(ba_lifestyle_cor_all_stand$p.fdr, digits = 4), scientific = TRUE)


openxlsx::write.xlsx(ba_lifestyle_cor_all_stand,file = './06newoutput/16table/f3_BALifestyleCor0123.xlsx')

#6.ba和慢性病====


#6.1先准备疾病df====
BaseFollowDiseaseBA <- read.csv("02derived_data/BaseFollowDiseaseBA.csv")
disease_fac<-BaseFollowDiseaseBA %>% dplyr::select(ID,Hypertension,Hyperlipidemia,Diabetes,Osteoporosis,CVD,Dementia)
#disease_fac<-full_join(BaseFollowDiseaseBA,cog_set)
disease_fac$MetabolicHealthsum<-disease_fac %>% dplyr::select(Hypertension,Hyperlipidemia,Diabetes) %>% rowSums()
disease_fac$MetabolicHealthsum<-as.factor(disease_fac$MetabolicHealthsum)


disease_fac$MetabolicHealth_2fac2<-ifelse(disease_fac$MetabolicHealthsum %in% c(0),0,
                                         ifelse(disease_fac$MetabolicHealthsum %in% c(1,2,3),1,NA))
disease_fac$MetabolicHealth_3fac<-ifelse(disease_fac$MetabolicHealthsum %in% c(0),0,
                                         ifelse(disease_fac$MetabolicHealthsum %in% c(1),1,
                                                ifelse(disease_fac$MetabolicHealthsum %in% c(2,3),2,NA)))

disease_fac$MetabolicHealth_3fac_ml<-ifelse(disease_fac$MetabolicHealth_3fac %in% c(0),0,
                                         ifelse(disease_fac$MetabolicHealth_3fac%in% c(1),1,NA))
disease_fac$MetabolicHealth_3fac_hl<-ifelse(disease_fac$MetabolicHealth_3fac %in% c(0),0,
                                            ifelse(disease_fac$MetabolicHealth_3fac %in% c(2),1,NA))
disease_fac$UnhealthMetabolic<-factor(disease_fac$MetabolicHealth_2fac2,levels = c(0,1),labels  = c('NonUnhealthMetabolic','UnhealthMetabolic'))
disease_fac$Moderate<-factor(disease_fac$MetabolicHealth_3fac_ml,levels = c(0,1),labels  = c('NonModerate','Moderate'))
disease_fac$High<-factor(disease_fac$MetabolicHealth_3fac_hl,levels = c(0,1),labels  = c('NonHigh','High'))
disease_fac$MetabolicHealth<-factor(disease_fac$MetabolicHealth_3fac,levels = c(0,1,2),labels  = c('healthy','moderately\nunhealthy',
                                                                                       'highly\nunhealthy'))
# disease_fac$Osteoporosis<-factor(disease_fac$Osteoporosis,levels = c(0,1,2),labels = c('healthy','osteopenia','osteoporosis'))
# disease_fac$CognitiveDecline<-factor(disease_fac$,levels = c(0,1,2),labels = c('healthy','mild','severe'))


table(disease_fac$UnhealthMetabolic)
table(disease_fac$MetabolicHealth_3fac)
table(disease_fac$Moderate)
table(disease_fac$High)
length(which(is.na(disease_fac$MetabolicHealth_3fac)))

disease_fac$Hyperlipidemia<-factor(disease_fac$Hyperlipidemia,levels = c(0,1),labels  = c('NonHyperlipidemia','Hyperlipidemia'))
disease_fac$Hypertension<-factor(disease_fac$Hypertension,levels = c(0,1),labels  = c('NonHypertension','Hypertension'))
disease_fac$Diabetes<-factor(disease_fac$Diabetes,levels = c(0,1),labels  = c('NonDiabetes','Diabetes'))
disease_fac$Osteoporosis<-ifelse(disease_fac$Osteoporosis>1,1,0)
disease_fac$Osteoporosis<-factor(disease_fac$Osteoporosis,levels = c(0,1),labels  = c('NonOsteoporosis','Osteoporosis'))
disease_fac$CVD<-factor(disease_fac$CVD,levels = c(0,1),labels  = c('NonCVD','CVD'))
disease_fac$Dementia<-factor(disease_fac$Dementia,levels = c(0,1),labels  = c('NonDementia','Dementia'))

write.csv(disease_fac,file = './06newoutput/12violinplot/disease_fac_1129.csv',row.names = F)

#disease_fac$MMSE<-multireplace(disease_fac$ID,dat2023_0905$ID,dat2023_0905$mmse_level)
#disease_fac$MMSE<-factor(disease_fac$MMSE,levels = c(0,1),labels  = c('NonDementia','Dementia'))
#disease_fac$Osteopenia<-factor(disease_fac$Osteopenia,levels = c(0,1,2),labels  = c('Normal','Osteopenia','Osteoporosis'))
#disease_fac$CognitiveDisorder<-factor(disease_fac$CognitiveDisorder,levels = c(0,1,2),labels  = c('Normal','MCI','Dementia'))
#forest_plot(cont_disease_res_sig)

bagap<-svm_res_0908[[4]]
bagap<-svm_res_0911[[4]]
bagap<-svm_res_0915[[4]]
bagap<-svm_res_1111[[4]]
bagap<-svm_res_1116[[4]]


descrTable(svm_res_0915[[4]])
descrTable(svm_res_0911[[4]])
descrTable(bagap_bodysepe)

bagap<-bagap %>% filter(sex=='1')
bagap<-bagap[complete.cases(bagap),]

#6.1 split volin====

#ba_disease$sex<-factor(ba_disease$sex,levels = c(1,2),labels  = c('Men','Women'))
# descrTable(MetabolicHealth_3fac~Cardiac+Bone+Gut+Kidney+Metab.+Brain+Body,data = ba_disease)
# descrTable(MetabolicHealth_2fac1~Cardiac+Bone+Gut+Kidney+Metab.+Brain+Body,data = ba_disease)
# descrTable(MetabolicHealth_2fac2~Cardiac+Bone+Gut+Kidney+Metab.+Brain+Body,data = ba_disease)
#

disease_fac<-read.csv('./06newoutput/12violinplot/disease_fac_1129.csv')
ba_disease <- left_join(bagaplm, disease_fac)
descrTable(MetabolicHealth~age,data = ba_disease)
correlation(data=ba_disease,select = 'Gut',select2 = 'age')
table(ba_disease$UnhealthMetabolic)
range(ba_disease$age[ba_disease$UnhealthMetabolic=='NonUnhealthMetabolic'],na.rm = T)
length()
hist(ba_disease$age[ba_disease$UnhealthMetabolic=='NonUnhealthMetabolic'])

colnames(ba_disease)
ba_disease_sum<-ba_disease %>%
  group_by(UnhealthMetabolic, CVD,Dementia,Osteoporosis) %>%
  summarize(n = n())
healthy_id<-ba_disease %>% filter(UnhealthMetabolic=='NonUnhealthMetabolic',
                                 # CVD=='NonCVD',
                                 # Dementia=='NonDementia',
                                 # Osteoporosis=='NonOsteoporosis'
                                 ) %>%
  dplyr::select(ID,age,sex)
healthyid<-healthy_id$ID
table(healthy_id$sex)
hist(healthy_id$age)
range(healthy_id$age,na.rm = T)

get_dat_long <- function(bagaplm, disease_fac) {
  library(tidyr)
  library(dplyr)
  ba_disease <- left_join(bagaplm, disease_fac)
  dat_long = ba_disease %>% pivot_longer(all_of(organ_simple),
                                         names_to = "method",
                                         values_to = "measure") %>%
    mutate(method = factor(method, levels = organ_simple, labels = organ_simple))
  #t检验之前先进行正态检验
  # 创建一个空的正态性检验结果数据框
  dat_long$measure[dat_long$method == 'Metab.'] %>% max(na.rm = T)
  dat_long_test <- dat_long %>%
    mutate(measure = ifelse(method == 'Metab.' &
                              as.numeric(measure) > 5, NA, measure))
  colnames(dat_long_test)[colnames(dat_long_test) == 'method'] <-
    'organ'
  return(dat_long_test)
}

##6.2代谢状态====
###6.2.0函数====
metabolichealth_res <- function(dat_long_test) {
  library(ggpubr)
  library(rstatix)
  dat_long_test$MetabolicHealth<-factor(dat_long_test$MetabolicHealth,
                                        levels = c('healthy',
                                                   'moderately\nunhealthy',
                                                   'highly\nunhealthy'))
  dat_long_test$measure_ind <-
    ifelse(is.na(dat_long_test$measure), F, T)
  wilcoxtest <-
    compare_means(
      data = dat_long_test,
      formula = measure ~ MetabolicHealth,
      group.by = 'organ' ,
      method = 'wilcox.test',
      p.adjust.method = 'none',
      alternative='greater',
    )#,ref.group = 'healthy'

  group_info <- dat_long_test %>%
    group_by(MetabolicHealth, organ, measure_ind) %>%
    summarize(n = n(),
              mean = mean(measure, na.rm = TRUE),
              var=var(measure,na.rm = T))
  group_info <-
    group_info %>% filter(measure_ind == T, !is.na(MetabolicHealth))
  group_info$measure_ind <- NULL

  metabolichealth_all <-
    left_join(group_info,
              wilcoxtest,
              by = c("MetabolicHealth" = "group1", 'organ' = 'organ'))
  colnames(group_info)[c(1, 3:5)] <- c('MetabolicHealth2', 'n2', 'mean2','var2')

  metabolichealth_all <-
    left_join(
      group_info,
      metabolichealth_all,
      by = c("MetabolicHealth2" = "group2", 'organ' = 'organ')
    )
  metabolichealth_all <- na.omit(metabolichealth_all)#metabolichealth_all[-(25:32), ]
  metabolichealth_select <-
    metabolichealth_all %>% dplyr::select(organ, MetabolicHealth,
                                          MetabolicHealth2, n, n2, mean, mean2, p)
}

metabolichealth_simple <- function(dat_long_test) {
  library(ggpubr)
  library(rstatix)

  wilcoxtest <-
    compare_means(
      data = dat_long_test,
      formula = measure ~ MetabolicHealth,
      group.by = 'organ' ,
      method = 'wilcox.test',
      p.adjust.method = 'none',
      alternative='greater',
    )#,ref.group = 'healthy'

  wilcox_effsize_res<-dat_long_test%>%
    group_by(organ) %>%
    wilcox_effsize(measure ~ MetabolicHealth,
                   p.adjust.method = "fdr",
                   alternative='less')

  wilcox_all<-full_join(wilcoxtest,wilcox_effsize_res)
  wilcox_all<-wilcox_all %>% mutate(
    p.adj.fdr=p.adjust(p,method='fdr')
  )
  wilcox_all %>% dplyr::select("organ","group1", "group2", "n1",
                               "n2", "p","p.adj.fdr","effsize")
}

multidisease_simple <- function(dat_long_test,disease=c("UnhealthMetabolic",#
                                                        "Dementia", "Osteoporosis")) {
  wilcox_res<-c()
  library(ggpubr)
  library(rstatix)
  for(i in disease){
    formula=as.formula(paste0('measure','~',i))
  wilcoxtest <-
    compare_means(
      data = dat_long_test,
      formula =formula ,
      group.by = 'organ' ,
      method = 'wilcox.test',
      p.adjust.method = 'none',
    )#,ref.group = 'healthy'

  wilcox_effsize_res<-dat_long_test%>%
    group_by(organ) %>%
    wilcox_effsize(formula,
                   p.adjust.method = "fdr")

  wilcox_all<-full_join(wilcoxtest,wilcox_effsize_res)
  wilcox_all<-wilcox_all %>% mutate(
    p.adj.fdr=p.adjust(p,method='fdr')
  )
  wilcox_all <-wilcox_all %>% dplyr::select("organ","group1", "group2", "n1",
                               "n2", "p","p.adj.fdr","effsize")
  wilcox_res<-rbind(wilcox_res,wilcox_all)
  }
return(wilcox_res)
}


bagaplm<-svm_res_0120[[1]][[4]]
bagap<-svm_res_0120[[1]][[6]]

bagaplm_man <-bagaplm %>% filter(sex==1)
bagaplm_woman <-bagaplm %>% filter(sex==2)

###6.2.1 gap和疾病====
dat_long_test<-get_dat_long(bagaplm, disease_fac)
dat_long_test$MetabolicHealth<-factor(dat_long_test$MetabolicHealth,
                                      levels = c('healthy',
                                                 'moderately\nunhealthy',
                                                 'highly\nunhealthy'))

###6.2.2代谢健康====
metabolichealth_all<-metabolichealth_res(dat_long_test)#有均值信息
metabolichealth_all$p.adj<-p.adjust(metabolichealth_all$p,method = 'fdr')

metabolichealth_simple_df<-metabolichealth_simple(dat_long_test)
metabolichealth_simple_output<-metabolichealth_simple_df %>% mutate(
  p=format(signif(p, digits = 4), scientific = TRUE) ,
  p.adj.fdr=format(signif(p.adj.fdr, digits = 4), scientific = TRUE),
  effsize=sprintf("%.3f", effsize)
)

dat_long_man<-get_dat_long(bagaplm_man, disease_fac)
metabolichealth_man<-metabolichealth_res(dat_long_man)

dat_long_woman<-get_dat_long(bagaplm_woman, disease_fac)
metabolichealth_woman<-metabolichealth_res(dat_long_woman)


write.xlsx(metabolichealth_all,file='./06newoutput/16table/f4_metabolichealth_all_0120.xlsx')
write.csv(metabolichealth_man,file='./06newoutput/16table/f4_metabolichealth_man.csv',row.names = F)
write.csv(metabolichealth_woman,file='./06newoutput/16table/f4_metabolichealth_woman.csv',row.names = F)

metabolichealth_list<-list(all=metabolichealth_all,man=metabolichealth_man,woman=metabolichealth_woman)
openxlsx::write.xlsx(metabolichealth_simple_output,file='./06newoutput/16table/f4_metabolichealth_0124.xlsx')


filename='./06newoutput/12violinplot/boxplot1128.pdf'
filename='./06newoutput/12violinplot/boxplot0120.pdf'
MetabolicHealth_boxplot_func(dat_long_test,filename)

##所有疾病====
dat_long_test$UnhealthMetabolic<-factor(dat_long_test$UnhealthMetabolic,
                                        levels = c('NonUnhealthMetabolic',
                                                   'UnhealthMetabolic'
                                        ))
dat_long_test$Dementia<-factor(dat_long_test$Dementia,
                               levels = c('NonDementia',
                                          'Dementia'),
                               labels = c('NonCognitiveImpairment',
                                          'CognitiveImpairment'))
dat_long_test$Osteoporosis<-factor(dat_long_test$Osteoporosis,
                                   levels = c('NonOsteoporosis',
                                              'Osteoporosis'))

multidisease_simple_df<-multidisease_simple(dat_long_test)

multidisease_simple_output<-multidisease_simple_df %>% mutate(
  p=format(signif(p, digits = 4), scientific = TRUE) ,
  p.adj.fdr=format(signif(p.adj.fdr, digits = 4), scientific = TRUE),
  effsize=sprintf("%.3f", effsize)
)
write.xlsx(multidisease_simple_output,file='./06newoutput/16table/f4_3disease_0124.xlsx')

disease<-c("Hypertension", "Diabetes", "Hyperlipidemia",#
           "Dementia", "Osteoporosis", "CVD")#
disease<-c("UnhealthMetabolic",#
           "Dementia", "Osteoporosis")
#disease<-c('Moderate','High',"UnhealthMetabolic","Dementia", "Osteoporosis", "CVD")
organs<-c('Cardiac','Bone','Gut','Kidney','Metab.','Brain','Cognition','Body')#,'Body.Brain','gut_BA_Gap'
normal_test_res<-normal_test_df_func(dat=dat_long_test,disease = disease,organs = organ_simple)
normal_man_res<-normal_test_df_func(dat=dat_long_man,disease = disease,organs = organ_simple)
normal_woman_res<-normal_test_df_func(dat=dat_long_woman,disease = disease,organs = organ_simple)
dput(colnames(normal_test_res))

output_name<-c("disease", "Organ", "normal_n", "abnormal_n", "mean_normal",
               "mean_abnormal","wilcoxtest_p","wilcoxtest_fdr", "Cohen.s.d", "Cohen.s.d.scientific")

disease4_gap_list<-list(all=normal_test_res[,output_name],
                    man=normal_man_res[,output_name],
                    woman=normal_woman_res[,output_name])

openxlsx::write.xlsx(disease4_gap_list,file='./06newoutput/16table/f4_disease4_gap_list_0120.xlsx')


filename='./06newoutput/12violinplot/violin1111_all_wilcoxtest.pdf'
filename='./06newoutput/12violinplot/violin1111_all_bodysepe_wilcoxtest.pdf'
filename='./06newoutput/12violinplot/violin1111_all_bodysepe_wilcoxtest_origin.pdf'

filename='./06newoutput/12violinplot/violin1116_other.pdf'
disease_name = c('Osteoporosis',"Dementia", 'CVD')##,,'Dementia'
violinmerge_func(dat_long, disease_name,filename,plabel='p.adj.label.star')#'p.label.star'

filename='./06newoutput/12violinplot/violin0120.pdf'
disease_name = c("UnhealthMetabolic",'Osteoporosis',"Dementia")##,,'Dementia', 'CVD'
violinmerge_func(dat_long_test, disease_name,filename,plabel='p.adj.label.star')#'p.label.star'


filename='./06newoutput/12violinplot/violin1007_all_wilcoxtest.pdf'
disease_name = c('Hypertension','Diabetes','Osteoporosis',"Dementia", 'CVD')##'Hyperlipidemia',,'Dementia'
violinmerge_func(dat_long, disease_name,filename,plabel='p.adj.label.star')#'p.label.star'


filename='./06newoutput/12violinplot/violin1007_all_wilcoxtest_pvalue.pdf'
violinmerge_func(dat_long, disease_name,filename,plabel='p.adj.scientific')#'p.label.star'


filename='./06newoutput/12violinplot/violin1007_intersect_wilcoxtest_fdr_test.pdf'
disease_name = c('Hypertension','Diabetes','Osteoporosis',"Dementia", 'CVD')##'Hyperlipidemia',,'Dementia'
violinmerge_func(dat_long, disease_name,filename,plabel='p.adj.label.star')

filename='./06newoutput/12violinplot/violin1007-0911_intersect_wilcoxtest_fdr_test.pdf'
disease_name = c('Hypertension','Diabetes','Osteoporosis',"Dementia", 'CVD')##'Hyperlipidemia',,'Dementia'
violinmerge_func(dat_long, disease_name,filename,plabel='p.adj.label.star')

filename='./06newoutput/12violinplot/violin0911_intersect_wilcoxtest_fdr_test.pdf'
disease_name = c('Hypertension','Diabetes','Osteoporosis',"Dementia")#, 'CVD'#'Hyperlipidemia',,'Dementia'
violinmerge_func(dat_long, disease_name,filename,plabel='p.adj.label.star')

filename='./06newoutput/12violinplot/violin0911_intersect_wilcoxtest_fdr_pvalue.pdf'
disease_name = c('Hypertension','Diabetes','Osteoporosis',"Dementia", 'CVD')#'Hyperlipidemia',,'Dementia'
violinmerge_func(dat_long, disease_name,filename,plabel='p.adj.label.num')


filename='./03output/09violinplot/violin0718_psignif.pdf'
violinmerge_func(dat_long, disease_name[1:6],filename)
filename='./03output/09violinplot/violin0718_pvalue.pdf'
violinmerge_func(dat_long, disease_name[1:6],filename,plabel = 'p.label.num')#"cohensd.label"#"p.label.num","p.signif"
filename='./03output/09violinplot/violin0718_cohend.pdf'
violinmerge_func(dat_long, disease_name[1:6],filename,plabel = 'cohensd.label')

filename='./03output/09violinplot/violin0722_max.pdf'
disease_name = c('Hyperlipidemia','Hypertension','Diabetes','Osteoporosis','Dementia','CVD')
violinmerge_func(dat_long, disease_name[1:6],filename)


filename='./03output/09violinplot/violin0714_2.pdf'
violinmerge_func(dat_long, disease_name[4:6],filename)

data_test<-dat_long
colnames(data_test)[10]<-'m'
data_test$Hypertension<-factor(data_test$Hypertension,
                               levels = c('Hypertension','NonHypertension'),
                               labels = c('normal','abnormal'))
wilcoxtest<-compare_means(data = data_test,formula = measure~Hypertension,group.by ='m' ,
                          method = 'wilcox.test',ref.group = 'normal')#
filter(data_test, !is.na(Hypertension)) %>%
  ggplot(., aes(x = m, y = measure,fill =Hypertension)) +#,
  introdataviz::geom_split_violin(alpha = .4,inherit.aes = T)+#,aes(x = m, y = measure,fill =Hypertension)
  geom_pwc(method = "wilcox_test", label = "p = {format(signif(p, digits = 3), scientific = TRUE)}{p.signif}",
           #aes(x = m, y = measure,fill =Hypertension),
           hide.ns=T)#
  stat_pvalue_manual(wilcoxtest, label = "p.adj",y.position = 3,x='m')




#7.中介分析====
##7.1选取变量====
dat2022_1017 <- read.csv("01input_data/dat2022_1017.csv")
  load(file = './03output/08lifestyle/lifestyle0920.Rdata')
dat2022_1017$BMI<-dat2022_1017$Weight_now/(dat2022_1017$Height_now*0.01)^2
dat2022_1017<-full_join(dat2022_1017,BaseFollowDiseaseBA[,c('ID','CVD')])
dat2022_1017<-full_join(dat2022_1017,lifestyle[,c('ID','sleep_cate')])

med_dat_origin<-dat2022_1017 %>% dplyr::select(ID,smoke,Identification,HTD2,
                                               DM2,drink,BMI,CVD,sleep_cate,Staple)

med_dat<-left_join(bagaplm,med_dat_origin)
colnames(med_dat)

##7.2协变量插补====
#取smoke,Identification,brain完全的数据集，其他的插补
# 创建一个新的数据框用于存储插补后的数据
imputed_df <- data.frame()
df<-med_dat
# 选择自变量、因变量和中介变量完整的样本
###7.2.1大脑====
complete_cases <- complete.cases(df$smoke, df$`Identification`, df$`Brain`) #Gap
complete_brain <- df[complete_cases, ] %>% dplyr::select(age,smoke,Brain,Identification,sex,HTD2,
                                                         DM2,drink,BMI,CVD,Staple,Cognition)#,sleep_cate

# 对其他存在缺失的变量进行中值插补
for (col in colnames(complete_brain)) {
  if (col != "smoke" && col != "Identification" && col != "Brain") {
    imputed_value <- median(complete_brain[[col]], na.rm = TRUE)  # 使用中值进行插补
    complete_brain[[col]] <- ifelse(is.na(complete_brain[[col]]), imputed_value, complete_brain[[col]])
  }
}
###7.2.2肾脏====
complete_cases <- complete.cases(df$Smoke, df$`Identification`, df$`Kidney`) #Gap
complete_Kidney <- df[complete_cases, ] %>% dplyr::select(age,smoke,Kidney,Identification,sex,HTD2,
                                                         DM2,drink,BMI,CVD,Staple,sleep_cate)

# 对其他存在缺失的变量进行中值插补
for (col in colnames(complete_Kidney)) {
  if (col != "smoke" && col != "Identification" && col != "Kidney") {
    imputed_value <- median(complete_Kidney[[col]], na.rm = TRUE)  # 使用中值进行插补
    complete_Kidney[[col]] <- ifelse(is.na(complete_Kidney[[col]]), imputed_value, complete_Kidney[[col]])
  }
}
###7.2.3身体====
complete_cases <- complete.cases(df$Smoke, df$`Identification`, df$`Body`) #Gap
which(complete_cases) %>% length
complete_Body <- df[complete_cases, ] %>% dplyr::select(age,smoke,Body,Identification,sex,HTD2,
                                                       DM2,drink,BMI,CVD,Cognition)

# 对其他存在缺失的变量进行中值插补
for (col in colnames(complete_Body)) {
  if (col != "smoke" && col != "Identification" && col != "Body") {
    imputed_value <- median(complete_Body[[col]], na.rm = TRUE)  # 使用中值进行插补
    complete_Body[[col]] <- ifelse(is.na(complete_Body[[col]]), imputed_value, complete_Body[[col]])
  }
}


#library(bruceR)

##7.3中介分析====
###7.3.1 brain====
dput(colnames(complete_brain))
complete_brain
Brain_PROCESS<-PROCESS(complete_brain, y = "Identification", x = "smoke",
        meds = c("Brain"),covs = c('sex','drink','HTD2',#'edu_y',
                                   'DM2','age','Cognition','Staple','BMI','CVD'),#,'total_income'
        ci = "boot", nsim = 1000, seed = 1)##
Brain_PROCESS

Brain_model <- lm(Identification ~ sex+smoke+Brain+HTD2+DM2+age+drink+BMI+CVD+Cognition+Staple, data=complete_brain)
summary(Brain_model) # 因变量和自变量关系
confint(Brain_model, 'Brain', level=0.95)

Brain_mediator <- lm(Brain ~ sex+smoke+HTD2+DM2+age+drink+BMI+CVD+Cognition+Staple, data=complete_brain)
summary(Brain_mediator) # 中介变量和自变量关系
confint(Brain_mediator, 'smoke', level=0.95)

Brain_mediation_model <-mediation::mediate(Brain_mediator,Brain_model,  treat = "smoke", mediator = "Brain",boot = T)
summary(Brain_mediation_model)

Brain_mediation_model$n0.ci#Prop. Mediated CI

###7.3.2 Kidney====
dput(colnames(complete_Kidney))
library(bruceR)
Kidney_PROCESS<-PROCESS(complete_Kidney, y = "Identification", x = "smoke",
                       meds = c("Kidney"),covs = c('sex','drink','HTD2',#'edu_y',
                                                  'DM2','age','BMI','CVD','sleep_cate','Staple'),#,'total_income'
                       ci = "boot", nsim = 1000, seed = 1)##
Kidney_PROCESS

Kidney_model <- lm(Identification ~ sex+smoke+Kidney+HTD2+DM2+age+drink+BMI+CVD+Staple+sleep_cate, data=complete_Kidney)
summary(Kidney_model) # 因变量和自变量关系
confint(Kidney_model, 'Kidney', level=0.95)

Kidney_mediator <- lm(Kidney ~ sex+smoke+HTD2+DM2+age+drink+BMI+CVD+Staple+sleep_cate, data=complete_Kidney)
summary(Kidney_mediator) # 中介变量和自变量关系
confint(Kidney_mediator, 'smoke', level=0.95)

Kidney_mediation_model <-mediation::mediate(Kidney_mediator,Kidney_model,  treat = "smoke", mediator = "Kidney",boot = T)
summary(Kidney_mediation_model)
Kidney_mediation_model$n0.ci#Prop. Mediated CI

###7.3.3 body====
dput(colnames(complete_Body))
Body_PROCESS<-PROCESS(complete_Body, y = "Identification", x = "smoke",
                       meds = c("Body"),covs = c('sex','drink','HTD2',#'edu_y',
                                                  'DM2','age','Cognition','BMI','CVD'),#,'total_income'
                       ci = "boot", nsim = 1000, seed = 1)##
Body_PROCESS

Body_model <- lm(Identification ~ sex+smoke+Body+HTD2+DM2+age+drink+Cognition+BMI+CVD, data=complete_Body)
summary(Body_model) # 因变量和自变量关系
confint(Body_model, 'Body', level=0.95)

Body_mediator <- lm(Body ~ sex+smoke+HTD2+DM2+age+drink+Cognition+BMI+CVD, data=complete_Body)
summary(Body_mediator) # 中介变量和自变量关系
confint(Body_mediator, 'smoke', level=0.95)

Body_mediation_model <-mediation::mediate(Body_mediator,Body_model,  treat = "smoke", mediator = "Body",boot = T)
summary(Body_mediation_model)
Body_mediation_model$n0.ci#Prop. Mediated CI
##★7.4通用中介函数====
#7.4.1先准备全部所需数据====
correlation(data = ba_lifestyle,select = 'uPDI',select2 = 'TUG')#测试，显著
load('./03output/08lifestyle/lifestyle1129.Rdata')
ba_lifestyle<-full_join(bagaplm,lifestyle)

disease_fac<-read.csv('./06newoutput/12violinplot/disease_fac_1129.csv')
disease_num<-disease_fac[,c('ID','MetabolicHealthsum')]
disease_num$MetabolicHealthsum<-as.numeric(disease_num$MetabolicHealthsum)

dat2022_1017 <- read.csv("01input_data/dat2022_1017.csv")
dat2022_1017$BMI<-dat2022_1017$Weight_now/(dat2022_1017$Height_now*0.01)^2

med_dat<-full_join(ba_lifestyle,disease_num)
med_dat<-full_join(med_dat,dat2022_1017[,c('ID','BMI')])

 #选择自变量、因变量和中介变量完整的样本
life_med_func <-
  function(data = med_dat,
           xvar = "smoke",
           mvar = 'Brain',
           yvar = 'Identification',
           covar = c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink')) {
    library(bruceR)
    library(mediation)
    c(xvar,mvar,yvar,covar) %in% colnames(med_dat)


    complete_cases <-
      complete.cases(med_dat[,c(xvar,mvar,yvar)]) #Gap


    complete_organ <-
      med_dat[complete_cases,] %>% dplyr::select(all_of(c(xvar,mvar,yvar,covar)))#,sleep_cate

    # 对其他存在缺失的变量进行中值插补
    for (col in colnames(complete_organ)) {
      if (col %in% covar) {
        imputed_value <-
          median(complete_organ[[col]], na.rm = TRUE)  # 使用中值进行插补
        complete_organ[[col]] <-
          ifelse(is.na(complete_organ[[col]]), imputed_value, complete_organ[[col]])
      }
    }
    # organ_PROCESS <-
    #   PROCESS(
    #     complete_organ,
    #     y = yvar,
    #     x = xvar,
    #     meds = mvar,
    #     covs =covar,
    #     ci = "boot",
    #     nsim = 1000,
    #     seed = 1
    #   )##
    # print(organ_PROCESS)


    outcome_formula <- as.formula(
      paste(yvar, "~", xvar, "+", mvar, "+", paste(covar, collapse = "+")))
    outcome_model <-
      lm(formula = outcome_formula,
        data = complete_organ
      )
    summary(outcome_model) # 因变量和自变量关系
    #confint(outcome_model, 'organ', level = 0.95)

    mediator_formula <- as.formula(
      paste(mvar, "~", xvar, "+", paste(covar, collapse = "+")))
    mediator_model <-
      lm(mediator_formula,
         data = complete_organ)

    summary(mediator_model) # 中介变量和自变量关系
   # confint(organ_mediator, 'smoke', level = 0.95)

    set.seed(123)
    organ_mediation_model <-
      mediation::mediate(
        mediator_model,
        outcome_model,
        treat = xvar,
        mediator = mvar,
        covariates = covar,
        boot = T
      )
   organ_mediation_model
  }
s_kidney<-organ_mediation_model
s_brain<-organ_mediation_model
trace(mediation:::print.summary.mediate,
      at = 11,
      tracer = quote({
        printCoefmat <- function(x, digits) {
          p <- x[, 4] #p-values seem to be stored rounded
          x[, 1:3] <- sprintf("%.4f", x[, 1:3])
          x[, 4] <-  sprintf("%.4f", p)
          print(x, quote = FALSE, right = TRUE)
        }
      }),
      print = FALSE)
mediation:::print.summary.mediate(summary(organ_mediation_model))
mediator_formula <- as.formula(
  paste('Kidney', "~", "smoke", "+", paste(c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink',
                                             'family_num','sleep_cate'), collapse = "+")))
set.seed(123)
s_kidney<-life_med_func(data = med_dat,
              xvar = "smoke",
              mvar = 'Kidney',
              yvar = 'Identification',
              covar = c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink',
                        'family_num','sleep_cate'))#'total_income',
summary(s_kidney)
mediation:::print.summary.mediate(summary(s_kidney))

mediator_formula <- as.formula(
  paste('Brain', "~", "smoke", "+", paste(c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink',
                                            'family_num'), collapse = "+")))
s_brain<-life_med_func(data = med_dat,
                       xvar = "smoke",
                       mvar = 'Brain',
                       yvar = 'Identification',
                       covar = c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink',
                                 'family_num'))#'total_income',
summary(s_brain)
mediation:::print.summary.mediate(summary(s_brain))
s_body<-life_med_func(data = med_dat,
              xvar = "smoke",
              mvar = 'Body',
              yvar = 'Identification',
              covar = c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink',
                        'family_num','income_average'))
summary(s_body)
mediation:::print.summary.mediate(summary(s_body))
untrace(mediation:::print.summary.mediate)

mediation_summary<-list(brain= capture.output(print(summary(s_brain))),
                        kidney= capture.output(print(summary(s_kidney))),
                        #body= capture.output(print(summary(s_body))),
                        brain1= capture.output(print(mediation:::print.summary.mediate(s_brain))),
                        kidney1= capture.output(print(mediation:::print.summary.mediate(s_kidney)))
                        #body1= capture.output(print(mediation:::print.summary.mediate(s_body)))
                        )


save(s_kidney,s_brain,s_body,file='06newoutput/11lifestyle/organ_mediation_model_summary240318.Rdata')
save(s_kidney,s_brain,file='06newoutput/11lifestyle/organ_mediation_model_summary240323.Rdata')


openxlsx::write.xlsx(mediation_summary, "06newoutput/16table/f3_mediation_summary_240323.xlsx")


life_med_func(data = med_dat,
              xvar = "uPDI",
              mvar = 'Cognition',
              yvar = 'TUG',
              covar = c('age', 'sex', 'BMI', 'MetabolicHealthsum', 'drink'))




##7.4特征描述====
med_dat_changename<-med_dat_origin
colnames(med_dat_changename)<-c('ID','Current Smoke','Olfactory Identification','Sex',#'Education Year',
                                'Hypertension','Diabetes','Current Alcohol Drinking',#, 'Hyperlipidemia',
                                'Current Tea Drinking','BMI','CVD')#,'Total Income Level',
complete_brain$`Current Smoke`<-factor(complete_brain$`Current Smoke`,levels = c(0,1),labels  = c('Smoker','Nonsmoker'))
complete_brain$Sex<-factor(complete_brain$Sex,levels = c(1,2),labels  = c('Men','Women'))
#complete_brain$Hyperlipidemia<-factor(complete_brain$Hyperlipidemia,levels = c(0,1),labels  = c('NonHyperlipidemia','Hyperlipidemia'))
complete_brain$Hypertension<-factor(complete_brain$Hypertension,levels = c(0,1),labels  = c('NonHypertension','Hypertension'))
complete_brain$Diabetes<-factor(complete_brain$Diabetes,levels = c(0,1),labels  = c('NonDiabetes','Diabetes'))
complete_brain$CVD<-factor(complete_brain$CVD,levels = c(0,1),labels  = c('NonCVD','CVD'))
complete_brain$`Current Alcohol Drinking`<-factor(complete_brain$`Current Alcohol Drinking`,levels = c(0,1),labels  = c('Yes','No'))
complete_brain$`Current Tea Drinking`<-factor(complete_brain$`Current Tea Drinking`,levels = c(0,1),labels  = c('Yes','No'))

descrTable(`Current Smoke`~.,data = complete_brain,hide = T)


#7.2肾脏====



#8.预测慢性病====
source("D:/work/01衰老/02标志物/20230417/42predict_HR_func.R")

base_med<-read.csv('./02derived_data/base_med.csv')
#预测CVE，需要矫正的协变量:年龄、性别、吸烟、体重指数、2型糖尿病、高血压治疗和他汀类药物的使用。
dat2022_1017 <- read.csv("01input_data/dat2022_1017.csv")
dat2022_1017$BMI<-dat2022_1017$Weight_now/(dat2022_1017$Height_now*0.01)^2
cve_covar<-dat2022_1017 %>% dplyr::select(ID,smoke,BMI)
cve_covar<-full_join(cve_covar,base_med)

##8.1 加上其他协变量====
#年龄、性别、吸烟、体重指数、2型糖尿病、高血压治疗和他汀类药物的使用
BaseFollowDiseaseBA <- read.csv("02derived_data/BaseFollowDiseaseBA_0722.csv")
BaseFollowDiseaseBA$date2 %>% table
length(which(!is.na(BaseFollowDiseaseBA$date2)))
median(BaseFollowDiseaseBA$date_diff2,na.rm = T)

f20_all_completeBA<-BaseFollowDiseaseBA#f20_all

bagaplm<-svm_res_0120[[1]][[4]]
bagap<-svm_res_0120[[1]][[6]]


f20_all_disease<-f20_all_completeBA[,-c(1,25:32)]
f20_all_covar<-left_join(round_df_func(bagap),round_df_func(f20_all_disease))
f20_all_covar<-left_join(round_df_func(bagaplm),round_df_func(f20_all_disease))

f20_all_covar<-left_join(f20_all_covar,cve_covar)
colnames(f20_all_covar)<-multireplace(colnames(f20_all_covar),organ_simple,organ_complete)

cal_NA(f20_all_covar)
covar<-c('sex','smoke','BMI','Diabetes','Hypertension','Hyperlipidemia','Hypotensive_Druguse','Hypoglycemic_Druguse')

correlation(f20_all_covar,select = 'cve2',select2 = 'Cardiovascular')

##有协变量，CA和BA+gap
f20_model<-CoxModelFunc_gapcovar(data = f20_all_covar,disease = 'cve',
                            organ_order = organ_complete,covar=covar,complete_organ=F)
CoxOutputFunc(f20_model,pathname_append = 'cagap',output = 'default',n=T)

##没有协变量，CA和BA+gap
f20_model<-CoxModelFunc_gap(data = f20_all_covar,disease = 'cve',organ_order = organ_complete,complete_organ=F)
CoxOutputFunc(f20_model,pathname_append = 'cagap',output = 'default',n=T)


##以下是BA和CA的对比
##交集，有协变量====
f20_model<-CoxModelFunc_bacovar(data = f20_all_covar,disease = 'cve',covar=covar,
                              organ_order = organ_complete,complete_organ=T)
CoxOutputFunc(f20_model,pathname_append = 'cagap',output = 'default',n=T)

##交集，无协变量====
f20_model<-CoxModelFunc_ba(data = f20_all_covar,disease = 'cve',
                         organ_order = organ_complete,complete_organ=T)
CoxOutputFunc(f20_model,pathname_append = 'cagap',output = 'default',n=T)

##全集，有协变量====
f20_model<-CoxModelFunc_bacovar(data = f20_all_covar,disease = 'cve',covar=covar,
                              organ_order = organ_complete,complete_organ=F)
CoxOutputFunc(f20_model,pathname_append = 'cagap',n=T,output = 'default')#

##全集，无协变量====
options(modelsummary_format_numeric_latex = "plain")

f20_model<-CoxModelFunc_ba(data = f20_all_covar,disease = 'cve',
                         organ_order = organ_complete,complete_organ=F)
CoxOutputFunc(f20_model,pathname_append = 'cagap',n=T,output ='default' )#'default'


##8.1.2森林图====
###1.先生成table====
bacacoxrestable<-CoxOutputFunc(f20_model,estimate =c('{estimate}'),
                               statistic = c('{conf.low}','{conf.high}','{p.value}'),
                               n=T,output ='data.frame' )#'default
#有协变量的
bacacoxrestable<-bacacoxrestable[,c(1:2,31:42)]
###2.画图====
bacacox_forest_p<-bacacox_forest_plot(bacacoxrestable)

pdf('06newoutput/17cox/bacacoxlm_forest_0120.pdf',width =6.6 ,height = 7)#_covar
print(bacacox_forest_p)
dev.off()

pdf('06newoutput/17cox/bacacoxlm_covar_forest_0120.pdf',width =6.6 ,height = 7)#_covar
print(bacacox_forest_p)
dev.off()

##8.2中介分析====
##8.2.1综合评分====
#a:  unhealthy lifestyles → aging; b: aging → CVE; c: unhealthy lifestyles → CVE
load('./03output/08lifestyle/lifestyle1111.Rdata')
dat2023_0905<-read.csv('./02derived_data/dat2023_0905.csv')
dput(colnames(lifestyle))
bagap_bodysepe
load("D:/work/01衰老/02标志物/20230417/01input_data/TIS_basedata20200806.Rdata")
cdat$exercise_frequency

unhealthy_lifestyle_df<-lifestyle %>%
  dplyr::select(ID,smoke,drink,hPDI_cate,uPDI_cate,PDI_cate)#,exercise_frequency
          #%>% rowSums()

unhealthy_lifestyle_df<-full_join(unhealthy_lifestyle_df,
                                  dat2023_0905[,c('ID','BMI')])
unhealthy_lifestyle_df<-full_join(unhealthy_lifestyle_df,
                                  cdat[,c('ID','exercise_frequency')],by='ID')
unhealthy_lifestyle_df$exercise_frequency<-as.numeric(unhealthy_lifestyle_df$exercise_frequency)
unhealthy_lifestyle_df$rev_exercise<-multireplace(unhealthy_lifestyle_df$exercise_frequency,c(1,2,3,4,5),c(5,4,3,2,1))

unhealthy_lifestyle_df$bmi_cate<-classify_bmi(unhealthy_lifestyle_df$BMI)


unhealthy_lifestyle_df<-unhealthy_lifestyle_df %>% mutate(
  hPDI_cate_st=standard_func(hPDI_cate),
  uPDI_cate_st=standard_func(uPDI_cate),
  PDI_cate_st=standard_func(PDI_cate),
  bmi_cate_st=standard_func(bmi_cate),
  rev_exercise_st=standard_func(rev_exercise)
)

unhealthy_lifestyle_df<-unhealthy_lifestyle_df %>% mutate(
  uPDI_cate_binary=ifelse(uPDI_cate<4,0,1),
  rev_exercise_binary=ifelse(rev_exercise<5,0,1)
)

##单项二值分类后，加权求和
unhealthy_lifestyle_df$score1<-weight_func(unhealthy_lifestyle_df,
            select = c('smoke','drink','uPDI_cate_binary','rev_exercise_binary'),
            weight=c(4,1,1,1))
unhealthy_lifestyle_df$score_cate<-my_cut(unhealthy_lifestyle_df$score)

##单项标准化后求和
unhealthy_lifestyle_df$score2<-unhealthy_lifestyle_df %>% select(smoke,drink,
                 uPDI_cate_st,rev_exercise_st)%>%#bmi_cate_st,
  rowSums()
unhealthy_lifestyle_df$score_cate<-my_cut(unhealthy_lifestyle_df$score2)

unhealthy_lifestyle_ba<-full_join(unhealthy_lifestyle_df,bagaplm)
correlation::correlation(data=unhealthy_lifestyle_ba,select = c('smoke','drink',
 'uPDI_cate_st','rev_exercise_st','score2','score_cate'),select2 = organ_simple,
                         p_adjust ='none')
correlation::correlation(data=unhealthy_lifestyle_ba,select = c('smoke'),
                         select2 = organ_simple,
                         p_adjust ='none')
correlation::correlation(data=unhealthy_lifestyle_ba,select = c('BMI'),select2 = 'smoke',
                         p_adjust ='none')
##LE8====
LE8<- read_excel("01input_data/LE8汇总.xlsx")
le_bagap<-full_join(LE8,bagap)
my_cut <- function(x) {
  q <- quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  cut(x, breaks = q, include.lowest = TRUE, labels = FALSE)
}
my_cut <- function(x) {
  q <- quantile(x, probs = c(0,1/2, 1), na.rm = TRUE)
  cut(x, breaks = q, include.lowest = TRUE, labels = FALSE)
}
le_bagap$le8_3<-my_cut(le_bagap$LE8)
correlation::correlation(data=le_bagap,select = c('LE8'),
                         select2 = organ_simple,
                         p_adjust ='none')
descrTable(le8_3~Cardiac+Bone+Gut+Kidney+Metab.+Brain,data=le_bagap)


#8.2.2吸烟中介

f20_all_covar$Bodyage<-f20_all_covar$Body+f20_all_covar$age
f20_all_covar_trans<-f20_all_covar
f20_all_covar_trans<-full_join(f20_all_covar_trans,LE8[,c('ID','LE8')])
f20_all_covar_trans<-f20_all_covar_trans %>% mutate(cve_surv := Surv(date_diff2, cve2))

#f20_all_covar_trans<-full_join(f20_all_covar_trans,unhealthy_lifestyle_df)
f20_all_covar_trans<-f20_all_covar_trans[complete.cases(f20_all_covar_trans[,
  c('cve_surv' , 'Cardiac','age','sex','LE8')]),]
c('cve_surv' , 'Cardiac','age','sex','LE8') %in% colnames(f20_all_covar_trans)
# 对其他存在缺失的变量进行中值插补
for (col in c('Diabetes','Hypertension',
              'Hypotensive_Druguse','Hypoglycemic_Druguse')) {
    imputed_value <- median(f20_all_covar_trans[[col]], na.rm = TRUE)  # 使用中值进行插补
    f20_all_covar_trans[[col]] <- ifelse(is.na(f20_all_covar_trans[[col]]),
                                         imputed_value, f20_all_covar_trans[[col]])

}


cve_Cardiac_model <- coxph( Surv(date_diff2, cve2)~ Cardiac+sex+smoke+Diabetes+Hypertension+#+Cardiac
                           Hypotensive_Druguse+Hypoglycemic_Druguse, data=f20_all_covar_trans)
cve_Cardiac_model <- survreg( cve_surv ~ LE8+age+smoke+sex+Diabetes+Hypertension+#+Cardiac
                             Hypotensive_Druguse+Hypoglycemic_Druguse, data=f20_all_covar_trans)
cve_Cardiac_model <- coxph( cve_surv ~ LE8+age+sex+Diabetes+Hypertension+#+Cardiac
                              Hypotensive_Druguse+Hypoglycemic_Druguse, data=f20_all_covar_trans)
cve_Cardiac_model <- coxph( cve_surv ~ Cardiac+LE8+age+sex+Diabetes+Hypertension+#+Cardiac
                              Hypotensive_Druguse+Hypoglycemic_Druguse, data=f20_all_covar_trans)
cve_Cardiac_model <- coxph( cve_surv ~ Cardiac+Bone+Brain+Kidney+Metab.+Body+age+sex+smoke+Diabetes+Hypertension+#+Cardiac
                              Hypotensive_Druguse+Hypoglycemic_Druguse, data=f20_all_covar_trans)
median(f20_all_covar_trans$date_diff2,na.rm = T)
summary(cve_Cardiac_model) # 因变量和自变量关系
modelsummary(cve_Cardiac_model,
             estimate = '{estimate}{stars}\n({conf.low}, {conf.high})',
             statistic = NULL,
             exponentiate = TRUE,
             coef_omit = "Intercept")

confint(cve_Cardiac_model, 'Cardiac', level=0.95)

cve_Cardiac_mediator <- lm(Cardiac~ sex+LE8+Diabetes+Hypertension+#+Cardiac
                          Hypotensive_Druguse+Hypoglycemic_Druguse, data=f20_all_covar_trans)
 # 因变量和自变量关系
summary(cve_Cardiac_mediator) # 中介变量和自变量关系
confint(cve_Cardiac_mediator, 'smoke', level=0.95)
class(cve_Cardiac_model)

Cardiac_mediation_model <-mediation::mediate(cve_Cardiac_mediator,cve_Cardiac_model,
                 treat = "LE8", mediator = "Cardiac")#,boot = T
summary(Cardiac_mediation_model)
Kidney_mediation_model$n0.ci#Prop. Mediated CI



table(unhealthy_lifestyle_df$smoke)#0   1   777 123
table(unhealthy_lifestyle_df$uPDI_cate)# 1   2   3   4   5  198 194 186 161 165
table(unhealthy_lifestyle_df$drink)# 0   1 640 258
table(unhealthy_lifestyle_df$exercise_frequency)# 1   2   3   4   5  772  11  25  36  60
length(which(is.na(unhealthy_lifestyle_df$BMI)))#缺失34
table(unhealthy_lifestyle_df$bmi_cate)#1   2   3  13 438 419

## 修正的心血管年龄====
f20_all_completeBA_modifyBP<-f20_all_completeBA
colnames(f20_all_completeBA_modifyBP)[which(colnames(f20_all_completeBA_modifyBP)=='Cardiovascular')]<-'Cardiovascular_origin'
f20_all_completeBA_modifyBP<-left_join(f20_all_completeBA_modifyBP,Cardiovascular_ba[[4]][,c(1,4)])


##函数====
##1.bmi划分
classify_bmi <- function(x) {
ifelse(x<18.5,1,ifelse(x<24,2,ifelse(28,3,4)))
}
standard_func<-function(x){
  min_x<-min(x,na.rm = T)
  max_x=max(x,na.rm = T)
  (x-min_x)/(max_x-min_x)
}
##2.四分位数划分
my_cut<-function(x){cut(x,quantile(x,na.rm = T),include.lowest=TRUE,labels=FALSE)}
#三分位数
my_cut <- function(x) {
  q <- quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  cut(x, breaks = q, include.lowest = TRUE, labels = FALSE)
}
##3.加权求和
weight_func<-function(data,select,weight){
  t<-as.matrix(data[,select])%*%weight/sum(weight)
}
#4.

#9.0====
#查看aa在cvd组间是否有差异

dat_long_disease2 = f20_all_completeBA %>% pivot_longer(all_of(organ_complete), names_to = "method",
                                       values_to = "measure") %>%
  mutate(method = factor(method,levels = organ_complete)) %>%
  filter(., CVD == 0 & !is.na(CVD2))

dat_long_disease2$CVD2<-factor(dat_long_disease2$CVD2,levels = c(0,1),labels = c('NonCVD2','CVD2'))
normal_test_df_func(dat_long_disease2,disease = 'CVD2',organs = organ_complete)

##9.前驱症状组是否年龄差更大====
dat_long_disease2 = f20_all_completeBA %>% pivot_longer(all_of(organ_complete), names_to = "method",
                                                        values_to = "measure") %>%
  mutate(method = factor(method,levels = organ_complete))
dat_long_disease3<-dat_long_disease2 %>% mutate(proHypertension=ifelse(Hypertension==0&Hypertension2==1,1,
                                                ifelse(Hypertension==0&Hypertension2==0,0,NA))) %>%
                                         mutate(proHyperlipidemia=ifelse(Hyperlipidemia==0&Hyperlipidemia2==1,1,
                                                ifelse(Hyperlipidemia==0&Hyperlipidemia2==0,0,NA))) %>%
                                        mutate(proDiabetes=ifelse(Diabetes==0&Diabetes2==1,1,
                                                ifelse(Diabetes==0&Diabetes2==0,0,NA))) %>%
                                        mutate(proCVD=ifelse(CVD==0&CVD2==1,1,
                                                ifelse(CVD==0&CVD2==0,0,NA)))
dat_long_disease3<- dat_long_disease3 %>% mutate(prothreehigh=ifelse((proHypertension+proHyperlipidemia+proDiabetes)>0,1,0))
                                                   #)
pro_gap_test<-normal_test_df_func(dat_long_disease3,disease =c('proHypertension','proHyperlipidemia',
                                                                'proDiabetes','proCVD','prothreehigh')
                                   ,organs = organ_complete,fac = F)

#10.基线药物数据====
base<-read_excel_allsheets("01input_data/Taizhou Imaging Study-Baseline-1230.xlsx")
base1<-base[[3]]
base2<-base[[4]]
base1_med<-base1 %>% dplyr::select(ID,Hypotensive_Druguse,Hypoglycemic_Druguse)#,Insulin_Druguse
base2_med<-base2 %>% dplyr::select(ID,HTD_medication,DM_medication)
colnames(base2_med)<-c('ID','Hypotensive_Druguse','Hypoglycemic_Druguse')
base_med<-rbind(base1_med,base2_med)
base_med[is.na(base_med)]<-0
write.csv('./02derived_data/base_med.csv',row.names = F)
#预测CVE，需要矫正的协变量:年龄、性别、吸烟、体重指数、2型糖尿病、高血压治疗和他汀类药物的使用。
dat2022_1017 <- read.csv("01input_data/dat2022_1017.csv")
dat2022_1017$BMI<-dat2022_1017$Weight_now/(dat2022_1017$Height_now*0.01)^2
cve_covar<-dat2022_1017 %>% dplyr::select(ID,smoke,BMI)
cve_covar<-full_join(cve_covar,base_med)
cal_NA(cve_covar)#smoke缺4个，bmi缺34个

dput(covar)
write.csv(base_med,file = './02derived_data/base_med.csv',row.names = F)

dat2023_0905<-full_join(dat2022_1017,base_med)
write.csv(dat2023_0905,file = './02derived_data/dat2023_0905.csv',row.names = F)
#附加函数====
#1.#识别df中的小数变量，保留小数====
round_df_func <- function(df, round_digit = 4) {
  computable_id <- sapply(df, class) %in% c("numeric") %>% which
  round_id <- c()
  for (i in computable_id) {
    if (any(df[, i] %% 1 != 0, na.rm = TRUE)) {
      round_id <- c(round_id, i)
    }
  }
  if (length(round_id) > 1) {
    df[, round_id] <-
      apply(df[, round_id], 2, function(x)
        round(x, round_digit))
  } else if (length(round_id) == 1) {
    df[, round_id] <- round(df[, round_id], round_digit)
  }
  return(df)
}
