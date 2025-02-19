#1.加载R包====
library(igraph)

#2.基准网络====
load('./02derived_data/Rdata/res_0424.Rdata')
res_ref<-res
sig_ref<-res_ref %>% filter(p.bonferroni<0.05)

adj_ref<-cor2adj_func(sig_ref[,1:3])
View(adj_ref)
net_ref<-graph_from_adjacency_matrix(adj_ref,mode="undirected",weighted=TRUE,diag=FALSE)


lc_0424<-read_excel_allsheets("02derived_data/lc_0424.xlsx")
lapply(lc_0424, function(x)nrow(x)) %>% unlist %>% sum

mod_ref<-lapply(lc_0424, function(x)x[['statname']])
names(mod_ref)<-1:length(mod_ref)


mem_df<-data.frame()
for (i in 1:length(mod_ref)) {
  statname <-mod_ref[[i]]
  clu <- i
  mem_df<-rbind(mem_df,data.frame(statname=statname,clu=clu))
}
mem<-multireplace(V(net_ref)$name,mem_df$statname,mem_df$clu) %>% as.numeric()



# 3.4计算模块间的连接强度====
#max(module_membership)

# 遍历网络的边
module_strength_func<-function(g){
  for (e in E(g)) {
    from_module <- V(g)[.from(e)]$clu_id
    to_module <- V(g)[.to(e)]$clu_id

    from_module_size <- sum(V(g)$clu_id == from_module)
    to_module_size <- sum(V(g)$clu_id == to_module)

    weight_correction <- sqrt(from_module_size * to_module_size)

    module_strength[from_module, to_module] <- module_strength[from_module, to_module] + E(g)$weight[e] / weight_correction
    module_strength[to_module, from_module] <- module_strength[to_module, from_module] + E(g)$weight[e] / weight_correction
    # 增加边的权重到连接强度矩阵中

    #module_strength[from_module, to_module] <- module_strength[from_module, to_module] + 1#无权，有方向
    #module_strength[to_module, from_module] <- module_strength[to_module, from_module] + 1#

  }
  return(module_strength)
}
module_strength <- matrix(0, nrow =7 , ncol = 7)
module_strength<-module_strength_func(organ_net)
colnames(module_strength)<-c('Cardiovascular','Bone','Brain','Gut','Kidney','Metabolic','Cognition')
rownames(module_strength)<-c('Cardiovascular','Bone','Brain','Gut','Kidney','Metabolic','Cognition')
#4.直接构建网络，重新矫正====
organ_ref<-res

organ_sig_ref<-organ_ref %>%
  filter(v1 %in% organ_statname_df$statname & v2 %in% organ_statname_df$statname) %>%
  mutate(p.fdr=p.adjust(p.value,method='bonferroni')) %>%#fdr
  filter(p.fdr<0.05)

organ_adj_ref<-cor2adj_func(organ_sig_ref[,1:3])
#View(organ_adj_ref)#287个节点
organ2_net<-graph_from_adjacency_matrix(organ_adj_ref,mode="undirected",weighted=TRUE,diag=FALSE)
V(organ2_net)$clu_id<-multireplace(V(organ2_net)$name,organ_statname_df$statname,
                                organ_statname_df$clu_id)%>% as.numeric()
#E(organ2_net)
module_strength2 <- matrix(0, nrow =7 , ncol = 7)
module_strength2<-module_strength_func(organ2_net)
colnames(module_strength2)<-c('Cardiovascular','Bone','Brain','Gut','Kidney','Metabolic','Cognition')
rownames(module_strength2)<-c('Cardiovascular','Bone','Brain','Gut','Kidney','Metabolic','Cognition')


V(organ2_net)$clu<-multireplace(V(organ2_net)$name,organ_statname_df$statname,organ_statname_df$Organ)
V(organ2_net)$color<-multireplace(V(organ2_net)$clu,col_df$organ,col_df$col)
createNetworkFromIgraph(organ2_net,"organ2Net")
plot(organ2_net)

#5.器官BA的偏相关性
load("D:/work/01衰老/02标志物/20230417/03output/svm_res_0602.Rdata")
ba<-svm_res[[4]]
ba<-ba[,-(1:3)]
colnames(ba)<-gsub('\nGap','',colnames(ba), fixed = TRUE)
ba<-na.omit(ba)
#library(ppcor)
ba_pcor<-cor(ba)

#矫正年龄性别====
load("D:/work/01衰老/02标志物/20230417/03output/svm_res_0602.Rdata")
ba<-svm_res[[4]]
colnames(ba)<-gsub('\nGap','',colnames(ba), fixed = TRUE)
ba_pcor<-ba_pcor_func(ba)
ba_pcor<-ba_pcor %>% filter(v1 %in% c('Cardiac','Bone','Brain','Gut','Kidney','Metab.')&
                              v2 %in% c('Cardiac','Bone','Brain','Gut','Kidney','Metab.'))
organBA_net<-graph_from_edgelist(as.matrix(ba_pcor[,1:2]),directed = F)
E(organBA_net)$weight<-ba_pcor$estimate
col<-c('#E69F00','#56B4E9','#009E73','#0072B2','#D55E00','#CC79A7')#,'#F0E442'
col_df<-data.frame(organ=c('Cardiac','Bone','Brain','Gut','Kidney','Metab.'),
                   col=col)#
V(organBA_net)$color<-multireplace(V(organBA_net)$name,col_df$organ,col_df$col)
E(organBA_net)$width<-abs(E(organBA_net)$weight*50)
(E(organBA_net)$weight) %>% range
# 设置颜色映射范围

color_range <- c(-0.2, 0.2)  # 连续型变量的最小值和最大值
colors <-  colorRampPalette(colors = c("blue","white","red"))(100)  # 颜色范围从蓝色到白色再到红色
values=E(organBA_net)$weight
# 将连续型变量映射到颜色

# Find color indices for values
color_indices <- findInterval(values, seq(fro=-0.2,to=0.2,length.out =100))

# Map color indices to colors
mapped_colors <- colors[color_indices]
E(organBA_net)$color<-mapped_colors
V(organBA_net)$frame.color<-NA

plot(organBA_net)
createNetworkFromIgraph(organBA_net,"organBANet")
