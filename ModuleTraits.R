library(WGCNA)
library(pracma)
library(ggplot2)
#1.modules~traits关联函数====
# 假设你的表型数据集是一个矩阵 named 'phenotype_data'
# 假设你的模块数据是一个列表 named 'modules'，其中每个元素是一个模块中的基因的向量
# 假设你的特征数据是一个矩阵 named 'feature_data'，其中每一列是一个特征

# 计算模块特征基因
moduleEigengenes <- function(data, modules) {
  MEList <- list()
  for (module in names(modules)) {
    moduleData <- data[, modules[[module]]]
    ME <- principalComponent(moduleData)
    MEList[[module]] <- ME
  }
  return(MEList)
}

principalComponent <- function(data) {
  # 使用PCA的第一主成分作为模块特征基因
  complete_id<-which(complete.cases(data))
  data<-na.omit(data)

  pca <- prcomp(data, na.action = na.exclude)
  module_pca<-numeric(904)
  module_pca[setdiff(1:904,complete_id)]<-NA
  module_pca[complete_id]<-pca$x[,1]
print(pca[["rotation"]][,1])
  return(module_pca)
}
data=dat_phenotype;modules=com_louvain
MEs <- moduleEigengenes(dat_all, com_louvain_test)
MEs <- moduleEigengenes(dat_phenotype, com_infomap)
unlist(com_louvain)[which(unlist(com_louvain) %notin% colnames(dat_phenotype))]

# 计算模块特征基因与特征之间的关联
correlations <- matrix(0, nrow=length(MEs), ncol=ncol(feature_data))
pvalues <- matrix(0, nrow=length(MEs), ncol=ncol(feature_data))


for (i in 1:length(MEs)) {
  for (j in 1:ncol(feature_data)) {
    result <- cor.test(MEs[[i]], feature_data[,j], use="complete.obs")

    correlations[i,j] <- result$estimate
    pvalues[i,j] <- result$p.value
  }
}

rownames(correlations)<-paste0('M',names(MEs))
colnames(correlations)<-colnames(feature_data)
rownames(pvalues)<-paste0('M',names(MEs))
colnames(pvalues)<-colnames(feature_data)
correlations_abs<-abs(correlations)

#2. modules~age====

MEs <- moduleEigengenes(dat_all, com_louvain)
correlations <- matrix(0, nrow=length(MEs), ncol=1)
pvalues <- matrix(0, nrow=length(MEs), ncol=1)
n_objs <- matrix(0, nrow=length(MEs), ncol=1)

for (i in 1:length(MEs)) {
  for (j in 1:1) {
    result <- cor.test(MEs[[i]], dat_all$age, use="complete.obs")
    n_objs[i,j]<-length(which(complete.cases(cbind(MEs[[i]], dat_all$age))))
    correlations[i,j] <- result$estimate
    pvalues[i,j] <- result$p.value
  }
}

rownames(correlations)<-paste0('M',names(MEs))
colnames(correlations)<-'age'
rownames(pvalues)<-paste0('M',names(MEs))
colnames(pvalues)<-'age'
rownames(n_objs)<-paste0('M',names(MEs))
colnames(n_objs)<-'n'

correlations_abs<-abs(correlations)

#2.可视化====
#library(reshape2)
#library(ggplot2)
cor_long <- melt(correlations)
pval_long <- melt(pvalues)
nval_long <- melt(n_objs)
heatmap_data <- cbind(cor_long, pval_long,n_objs)[,-4:-5]
colnames(heatmap_data)<-c('Module','Trait','Correlation','pvalue','n.objs')
heatmap_data$fdr<-p.adjust(heatmap_data$pvalue,method = 'fdr')
heatmap_data$label =sprintf("%.3f", heatmap_data$Correlation)
#heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.1 & heatmap_data$fdr >= 0.05] = paste0(heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.1 & heatmap_data$fdr >= 0.05], "\n+")
heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.05 & heatmap_data$fdr >= 0.01] = paste0(heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.05 & heatmap_data$fdr >= 0.01], "\n*")
heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.01 & heatmap_data$fdr >= 0.001] = paste0(heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.01 & heatmap_data$fdr >= 0.001], "\n**")
heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.001 & heatmap_data$fdr >= 0.0001] = paste0(heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.001 & heatmap_data$fdr >= 0.0001], "\n***")
heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.0001 ] = paste0(heatmap_data$label[!is.na(heatmap_data$fdr) & heatmap_data$fdr < 0.0001], "\n****")

# 生成热图
heatmap_plot <- ggplot(data = heatmap_data, aes(x = Module, y = Trait, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_text(aes(label = label), color = "black", size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text = element_text(size = 8))+
  labs(x = "", y = "")+
  labs(fill = "Correlation Coefficient")+
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 5))
heatmap_plot

#3.模块之间的相关性====

MEdf<-reduce(MEs,cbind)
colnames(MEdf)<-paste0('M',1:ncol(MEdf))

MEdf<-as.data.frame(MEdf)
#MEdf_cor<-correlation(MEdf,select = colnames(MEdf),na.rm=T,p_adjust='none',method = "spearman")

correlation_me <- cor(MEdf[,c(1,3,7,12,15)],use = "pairwise.complete.obs")#计算相关系数
p_value_matrix <- matrix(nrow = 5, ncol = 5)#创建空矩阵来储存P值
for (i in 1:5) {
  for (j in 1:5) {
    if (i == j) {
      p_value_matrix[i, j] <- NA  # 代谢物与自身的P值设为NA
    } else {
      p_value <- cor.test(MEdf[,i], MEdf[,j])$p.value
      p_value_matrix[i, j] <- p_value
    }
  }
}
heatmap(correlation_me)
p_value_matrix

MEagecor<-c()
for(i in 1:15) {
  testdf <- dat_all[, c(com_louvain[[i]], 'age')]
  testdf <- na.omit(testdf)
  cancor_res <- cancor(testdf[, com_louvain[[i]]], testdf[, 'age'])
  MEagecor<-c(MEagecor,cancor_res$cor)
}
MEagecor
#install.packages('CCP')
trans_louvain<-list(com_louvain[[1]],com_louvain[[2]],com_louvain[[11]],
                    com_louvain[[14]],com_louvain[[6]], com_louvain[[12]],
                    com_louvain[[15]],com_louvain[[5]],com_louvain[[3]],
                    com_louvain[[9]],com_louvain[[10]],com_louvain[[4]],
                    com_louvain[[8]],com_louvain[[7]],com_louvain[[13]])
library(CCP)
trans_louvain<-read_excel_allsheets("06newoutput/16table/f1_trans_louvain.xlsx")
ME_cancor<-matrix(nrow = 16, ncol = 16)
ME_cancor_pvalue<-matrix(nrow = 16, ncol = 16)
com_louvain_age<-c('age',trans_louvain)

for(i in 1:16) {
  for(j in 1:16){
    testdf<-dat_all[,c(com_louvain_age[[i]],com_louvain_age[[j]])]
    testdf<-na.omit(testdf)
    cancor_res <- cancor(testdf[,com_louvain_age[[i]]],testdf[,com_louvain_age[[j]]])$cor

    n=nrow(testdf)
    p=length(com_louvain_age[[i]])
    q=length(com_louvain_age[[j]])
    if(length(cancor_res)<min(p,q)){cancor_res=c(cancor_res,cancor_res[
      length(cancor_res)])}

    p_roy<-p.asym(cancor_res, n, p, q, tstat = "Wilks")#Wilks Hotelling
    #

    if(length(p_roy$p.value)>1){ME_cancor_pvalue[i,j]<-p_roy$p.value[2]
    }else{
      ME_cancor_pvalue[i,j]<-p_roy$p.value[length(p_roy$p.value)]
    }

    ME_cancor[i,j]<-cancor_res[which(abs(cancor_res)==min(abs(cancor_res)))[1]]
    m<-m+1
    }
}

colnames(ME_cancor)<-c('age',paste0('M',1:15))
rownames(ME_cancor)<-c('age',paste0('M',1:15))


colnames(ME_cancor_pvalue)<-c('age',paste0('M',1:15))
rownames(ME_cancor_pvalue)<-c('age',paste0('M',1:15))

ME_cancor<-ME_cancor[,-1]
ME_cancor_pvalue<-ME_cancor_pvalue[,-1]
ME_cancor_pvalue_label<-apply(ME_cancor_pvalue,2,p_to_star)
#ME_cancor_pvalue<-apply(ME_cancor_pvalue,2,p_to_star)

row_split=c('',rep('  ',15))

p.asym(0.4, 200, 3, 9, tstat = "Wilks")

##3.2模块和年龄之间的相关性====
library(CCP)
trans_louvain<-read_excel_allsheets("06newoutput/16table/f1_trans_louvain.xlsx")
Mage_cancor<-matrix(nrow = 15, ncol = 1)
Mage_cancor_pvalue<-matrix(nrow = 15, ncol = 1)
com_louvain_age<-c('age',trans_louvain)


Mage_cancor_func <- function(dat_all, trans_louvain) {
  Mage_cancor <- data.frame(cor = numeric(15), p = numeric(15))
  for (i in 1:15) {
    testdf <- dat_all[, c(t(trans_louvain[[i]]), 'age')]
    testdf <- na.omit(testdf)
    cancor_res <-
      cancor(testdf[, t(trans_louvain[[i]])], testdf[, 'age'])$cor

    n = nrow(testdf)
    p = length(t(trans_louvain[[i]]))
    q = 1
    if (length(cancor_res) < min(p, q)) {
      cancor_res = c(cancor_res, cancor_res[length(cancor_res)])
    }

    p_roy <-
      p.asym(cancor_res, n, p, q, tstat = "Wilks")#Wilks Hotelling
    Mage_cancor[i, 2] <- p_roy$p.value[length(p_roy$p.value)]
    Mage_cancor[i, 1] <- cancor_res[[1]]
  }
  return(Mage_cancor)
}
Mage_cancor_all<-Mage_cancor_func(dat_all,trans_louvain)
Mage_cancor_man<-Mage_cancor_func(dat1_all,trans_louvain)
Mage_cancor_woman<-Mage_cancor_func(dat2_all,trans_louvain)

Mage_cancor<-data.frame(all=Mage_cancor_all$cor,man=Mage_cancor_man$cor,woman=Mage_cancor_woman$cor)
Mage_p<-data.frame(all=Mage_cancor_all$p,man=Mage_cancor_man$p,woman=Mage_cancor_woman$p)
Mage_p_adj<-apply(Mage_p,2,function(x)p.adjust(x,method = 'fdr')) %>% as.data.frame()
Mage_plabel<-apply(Mage_p,2,p_to_star) %>% as.data.frame()
Mage_padj_label<-apply(Mage_p_adj,2,p_to_star) %>% as.data.frame()


#4.可视化年龄、模块之间的关系====
filename='06newoutput/02agecor/moudulescor1022.pdf'
pdf(filename,width = 5,height = 6)
col_fun = colorRamp2(c(0, 1), c( "#EEEEEE", "#E64B35B2"))
p1<-Heatmap(ME_cancor,col = col_fun,name = 'cor',
            rect_gp = gpar(col = "white"), cluster_rows = F,
            border = F,
            cluster_columns = F,
            column_names_rot=45,
            row_names_gp = gpar(fontsize = 15),
            column_names_gp = gpar(fontsize = 12),
            show_heatmap_legend = T,
            row_split=row_split,
            #row_split=lifestyle_label,
)#,row_order = names(lifestyle_label),
draw(p1)
dev.off()

p_to_star <- function(p_values) {
  ifelse(p_values=='', "",
         ifelse(p_values < 0.001, "***",
                ifelse(p_values < 0.01, "**",
                       ifelse(p_values < 0.05, "*",
                              ifelse(p_values < 0.1, "+",'')))))
}
##5.年龄相关性，稳健性得分====
cluster_similarity_df_long<-read.csv('./06newoutput/100louvainInfomapLeiden/1031cluster_similarity_df_long.csv')

##5.1只有一个总体模块和年龄相关性====
sta_cor_df<-data.frame( `Correlation\nwith age`= ME_cancor[1,],
  `Moudules\nstability`=cluster_similarity_df_long[cluster_similarity_df_long$Method=='Mean',c(3)]
)
sta_cor_df<-as.matrix(sta_cor_df)
sta_cor_df<-apply(sta_cor_df,2,as.numeric)
rownames(sta_cor_df)<-paste0('M',1:15)
colnames(sta_cor_df)<-c('','')#'Correlation#\nwith age','Stability'


sta_cor_df_label<-data.frame(
  correlation=paste0(sprintf("%.2f", ME_cancor[1,]),'\n',ME_cancor_pvalue_label[1,]),
  stability=sprintf("%.2f", sta_cor_df[,2]))

sta_cor_df_label<-data.frame(
  correlation=paste0(ME_cancor_pvalue_label[1,]),
  stability='')

##5.2稳定性，总男女模块年龄相关性====
stability=cluster_similarity_df_long[cluster_similarity_df_long$Method=='Mean',c(3)]
#5.2.1 cor matrix
sta_cor_df<-cbind(stability,Mage_cancor)
sta_cor_df<-as.matrix(sta_cor_df)
sta_cor_df<-apply(sta_cor_df,2,as.numeric)
rownames(sta_cor_df)<-paste0('M',1:15)
colnames(sta_cor_df)[1]<-c('')#'Correlation#\nwith age','Stability'
sta_cor_df<-t(sta_cor_df)
#5.2.1 label matrix
sta_cor_df_label<-cbind(stability,Mage_padj_label)
sta_cor_df_label$stability<-''
sta_cor_df_label<-as.matrix(sta_cor_df_label)
rownames(sta_cor_df_label)<-paste0('M',1:15)
colnames(sta_cor_df_label)[1]<-c('')#'Correlation#\nwith age','Stability'
sta_cor_df_label<-t(sta_cor_df_label)


module_list<-list(ModuleStability=stability,ModuleAgeCor=Mage_cancor,
                  ModuleAgeP=Mage_p,ModuleAgeFDR=Mage_p_adj)
openxlsx::write.xlsx(module_list,rownames=T,file = '06newoutput/16table/f1_module_Stability&Agecor.xlsx')

#5.3画图====
col1_test<-colorRamp2(c( 0, 1), c( "#3C5488B2", "#8491B4B2"))
col1_test(0.6)

row_split=c('Stability','Correlation','Correlation','Correlation')
col_fun = colorRamp2(c( 0, 1), c( "#EEEEEE",'#E64B35B2'))# "#F21A00"
col1_fun = colorRamp2(c( 0, 1), c( "#EEEEEE", "#6878A2FF"))
p1<-Heatmap(sta_cor_df,col = col_fun,name = 'Correlation\nwith Age',
            rect_gp = gpar(col = "white"), cluster_rows = F,border = F,
            cluster_columns = F,
            column_names_rot=45,
            row_names_gp = gpar(fontsize = 11),
            column_names_gp = gpar(fontsize = 11),
            show_heatmap_legend = T,
            #column_split = col_split,
            row_split=row_split,
            row_gap = unit(6, "mm"),
           # column_names_side = "top",
            cell_fun = function(j, i, x, y, width, height, fill) {
              if(i==1) {

                grid.rect(x = x, y = y, width = width, height = height,
                          gp = gpar(fill = col1_fun(sta_cor_df[i, j]),col = "white"))

              } else{
              grid.rect(x = x, y = y, width = width, height = height,
                        gp = gpar( fill = col_fun(sta_cor_df[i, j]),col = "white"))}#,
              grid.text( sta_cor_df_label[i, j], x, y, gp = gpar(fontsize = 12))
            } ,
            #row_split=lifestyle_label,
)#,row_order = names(lifestyle_label),
lgd = Legend(col_fun = col1_fun,at = c(0,  0.5,  1),
             title = expression(bold("Module\nStability")),
             legend_width = unit(2, "cm"),grid_height =unit(0.4, "cm"),
             labels_gp = gpar(fontsize=9),title_gp = gpar(fontsize=10))

filename='06newoutput/02agecor/StabilityCorrealtion_allmanwoman.pdf'
pdf(filename,width = 6,height = 4)#width = 2,height = 3.8
draw(p1, heatmap_legend_side = "left")
draw(lgd, x = unit(0.06, "npc"), y = unit(0.15, "npc"))
dev.off()

##6.模块关联网络====
##6.1====
mnet<-graph_from_adjacency_matrix(ME_cancor[-1,],mode="undirected",weighted=TRUE,diag=FALSE)
#E(mnet)$width<-E(mnet)$weight*10
range(E(mnet)$weight*10)
V(mnet)$frame.color<-NA
V(mnet)$label.family<-'Helvetica'
#V(mnet)$size<-strength(mnet)*5
col_fun = colorRamp2(c(0, 0.5), c( "#F3DBD5FF", "#DC0000B2"))
E(mnet)$weight %>% max
E(mnet)$color<-col_fun(E(mnet)$weight)
#V(mnet)$label<-''

##6.2关联p值构建的网络====
ME_cancor_pvalue_trans<-ME_cancor_pvalue
ME_cancor_pvalue_trans[ME_cancor_pvalue_trans==0]<-10e-5
mnet_pvalue<-graph_from_adjacency_matrix(ME_cancor_pvalue_trans[-1,],mode="undirected",weighted=TRUE,diag=FALSE)
E(mnet)$width<-ifelse(E(mnet_pvalue)$weight=='', 0,
                      ifelse(E(mnet_pvalue)$weight < 0.001, 8,
                             ifelse(E(mnet_pvalue)$weight < 0.01, 6,
                                    ifelse(E(mnet_pvalue)$weight < 0.05, 4,
                                           ifelse(E(mnet_pvalue)$weight < 0.1, 2,0.5))))) %>% as.numeric()

#edges_to_delete <- E(mnet)[E(mnet)$width == 0]
#mnet_delete <- delete_edges(mnet, edges_to_delete)

filename='06newoutput/02agecor/modulesnet_12.pdf'
pdf(filename,width = 6,height = 6)

lay_mnet<-layout.circle(mnet)
plot(mnet,layout=lay_mnet)

lgd = Legend(col_fun = col_fun,
             title = "Correlation coefficient",
             direction = "horizontal",legend_width = unit(4, "cm"),grid_height =unit(0.5, "cm"),
             labels_gp = gpar(fontsize=10),title_gp = gpar(fontsize=12))
draw(lgd, x = unit(0.55, "npc"), y = unit(0.12, "npc"))
legend("top",#xy.coords(x = unit(0.15, "npc"), y = unit(0.85, "npc")),
       legend = c("***", "**",
                  "*","+",
                  "ns"),
       lwd = c(8,6,4,2,0.5),
       cex=1.5,
       col = 'grey',
       ncol=5,
       bty='n',
       x.intersp = 0.5, seg.len = 0.8, y.intersp = 1)

dev.off()


# cell_fun = function(j, i, x, y, width, height, fill) {
#   grid.rect(x = x, y = y, width = width, height = height,
#             gp = gpar( fill = col_fun(round(ME_cancor[i, j],2)),col = "white"))#,
#   grid.text( round(ME_cancor[i, j],2), x, y, gp = gpar(fontsize = 12))
# }
# lgd = Legend(col_fun = col_fun,
#              title = "Correlation coefficient",
#              direction = "horizontal",legend_width = unit(4, "cm"),grid_height =unit(0.8, "cm"),
#              labels_gp = gpar(fontsize=15),title_gp = gpar(fontsize=15))
# draw(lgd, x = unit(0.86, "npc"), y = unit(0.14, "npc"))


#ME_cancor<-as.data.frame(ME_cancor)
ME_cancor_age_related<-ME_cancor[c(1,3,7,12,15),c(1,3,7,12,15)]
rownames(ME_cancor_age_related)<-paste0('M',c(1,3,7,12,15))
colnames(ME_cancor_age_related)<-paste0('M',c(1,3,7,12,15))
heatmap(ME_cancor)
heatmap(ME_cancor_age_related)

testdf<-dat_all[,c(com_louvain[[1]],com_louvain[[7]])]
testdf<-na.omit(testdf)
cancor_obj<-cancor(testdf[,com_louvain[[1]]],testdf[,com_louvain[[7]]])
summary(cancor_obj)$correlations$P[1]

#ME_cancor<-as.data.frame(ME_cancor)
heatmap(ME_cancor)

#为什么心血管模块单个变量和年龄显著，PCA之后不显著？
cor.test(dat_phenotype[,com_louvain[[12]]], use="complete.obs")
cor.test(ME,dat_phenotype$age, use="complete.obs")
head(dat_phenotype$age)
head(dat_all$age)
head(feature_data$age)
dat_phenotype
cor_vascular<-correlation(dat_phenotype,select = com_louvain[[1]],select2 = 'age')
cor_bone<-correlation(dat_phenotype,select = com_louvain[[12]][7:12],select2 = 'age')
ME=principalComponent(dat_phenotype[,com_louvain[[12]]][7:12])
cor_bone<-correlation(dat_all,select = com_louvain[[3]],select2 = 'Identification')

#4.测试ccp
X <- matrix(rnorm(150), 50, 3)
Y <- matrix(rnorm(250), 50, 5)
## Calculate canonical correlations:
rho <- cancor(X,Y)$cor
## Define number of observations,
## and number of dependent and independent variables:
N = dim(X)[1]
p = dim(X)[2]
q = dim(Y)[2]

## Calculate p-values using F-approximations of some test statistics:
p.asym(rho, N, p, q, tstat = "Wilks")

res1 <- p.asym(rho, N, p, q)
plt.asym(res1,rhostart=1)
