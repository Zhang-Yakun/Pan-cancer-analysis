
setwd("f:/workplace/结果1/")

library(dplyr)
library(data.table)

pss<-read.csv("PSS.csv")
pathway16<-pss$description%>%unique()

s_df<-data.frame()

for(i in 1:length(pathway16)){
  s<-subset(pss,pss$description==pathway16[i])
  s2<-subset(s,abs(s$S)>6)
  n<-length(s2$pathwayID)
  df<-data.frame(pathway=pathway16[i],cancer_number=n)
  s_df<-rbind(s_df,df)
}

pas<-read.csv("PAS.csv")
pathway16<-pas$description%>%unique()

a_df<-data.frame()

for(i in 1:length(pathway16)){
  a<-subset(pas,pas$description==pathway16[i])
  a2<-subset(a,abs(a$PAS)>6)
  n<-length(a2$pathwayID)
  df<-data.frame(pathway=pathway16[i],cancer_number=n)
  a_df<-rbind(a_df,df)
}


# 通路变异数 -------------------------------------------------------------------

setwd("f:/workplace/结果1/")

gene<-fread("GeneSymbol_90.txt")%>%as.data.frame()%>%distinct()
g<-gene$ENTREZID


kegg<-read.csv("HSA_KEGG.csv")
p7<-c("hsa04660","hsa04650","hsa05340","hsa04064","hsa04150","hsa04012","hsa04310")

for(i in 1:7){
  path<-subset(kegg,kegg$KEGGID==p7[i])
  b<-intersect(path$ENTREZID,g)
  path2<-subset(path,path$ENTREZID%in%b)
  path3<-subset(gene,gene$ENTREZID%in%path2$ENTREZID)
  write.csv(path3,file=paste(p7[i],".csv",sep = ""),row.names = F)
}

########通路突变统计
setwd("F:/immuno2.0/result2/")

cnv_matrix<-read.csv("泛癌-47基因拷贝数频数.csv")

mut_matrix<-read.csv("泛癌-47基因突变频数.csv")

###############拷贝数总数统计
row.names(cnv_matrix)<-cnv_matrix$cancer
cnv_matrix<-subset(cnv_matrix , select = - cancer)
cnv_matrix<-subset(cnv_matrix , select = - case)
cnv_matrix<-as.matrix(cnv_matrix)
cnv_result<-apply(cnv_matrix,1,sum)
cnv_result<-sort(cnv_result,decreasing = T)
write.csv(cnv_result,"泛癌-基因拷贝总数.csv")

colnames(cnv_matrix)<-stringr::str_replace(colnames(cnv_matrix),"[.]","-")

###############突变总数统计
row.names(mut_matrix)<-mut_matrix$cancer
mut_matrix<-subset(mut_matrix , select = - cancer)
mut_matrix<-subset(mut_matrix , select = - cases)
mut_matrix<-as.matrix(mut_matrix)
mut_result<-apply(mut_matrix,1,sum)
mut_result<-sort(mut_result,decreasing = T)
write.csv(mut_result,"泛癌-基因突变总数.csv")

library(stringr)
colnames(mut_matrix)<-stringr::str_replace(colnames(mut_matrix),"[.]","-")

#################通路中基因突变总数

library(dplyr)

file_name<-list.files("f:/immuno2.0/result2/pathway9/")

setwd("f:/immuno2.0/result2/pathway9/")

for(i in 11:16){
  input<-read.csv(file_name[i])
  pathway_gene<-input$SYMBOL%>%unique()
  mut_matrix2<-mut_matrix[,pathway_gene]
  
  mut_result2<-apply(mut_matrix2,1,sum)
  mut_result2<-sort(mut_result2,decreasing = T)
  mut_result3<-as.data.frame(mut_result2)
  mut_result3$CancerType<-rownames(mut_result3)
  colnames(mut_result3)<-c("MutNumber","CancerType")
  write.csv(mut_result3,paste("泛癌基因突变总数-",file_name[i],sep = ""),row.names = F)
}
#########

mut_df<-data.frame()

for(i in 1:15){
  input<-read.csv(file_name[i])
  pathway_gene<-input$SYMBOL%>%unique()
  mut_matrix2<-mut_matrix[,pathway_gene]
  
  mut_result2<-apply(mut_matrix2,1,sum)
  mut_result2<-sort(mut_result2,decreasing = T)
  mut_result3<-as.data.frame(mut_result2)
  mut_result3$CancerType<-rownames(mut_result3)
  colnames(mut_result3)<-c("MutNumber","CancerType")
  mut_sum<-sum(mut_result3$MutNumber)
  
  mut_res4<-data.frame(pathway=file_name[i],mutSum=mut_sum)
  
  mut_df<-rbind(mut_df,mut_res4)
}

mut_df<-mut_df[order(mut_df$mutSum),]
write.csv(mut_df,"通路免疫衰老基因的总突变数.csv",row.names = F)

#################通路中基因CNV总数

file_name<-list.files("f:/immuno2.0/result2/pathway9/")

setwd("f:/immuno2.0/result2/pathway9/")

cnv_matrix<-cnv_matrix[-34,]

for(i in 1:15){
  input<-read.csv(file_name[i])
  pathway_gene<-input$SYMBOL%>%unique()
  path_gene<-intersect(colnames(cnv_matrix),pathway_gene)
  cnv_matrix2<-cnv_matrix[,path_gene]
  
  cnv_result2<-apply(cnv_matrix2,1,sum)
  cnv_result2<-sort(cnv_result2,decreasing = T)
  cnv_result3<-as.data.frame(cnv_result2)
  cnv_result3$CancerType<-rownames(cnv_result3)
  colnames(cnv_result3)<-c("MutNumber","CancerType")
  write.csv(cnv_result3,paste("泛癌基因CNV总数-",file_name[i],sep = ""),row.names = F)
}

#########

cnv_df<-data.frame()

for(i in 1:15){
  input<-read.csv(file_name[i])
  pathway_gene<-input$SYMBOL%>%unique()
  path_gene<-intersect(colnames(cnv_matrix),pathway_gene)
  cnv_matrix2<-cnv_matrix[,path_gene]
  
  cnv_result2<-apply(cnv_matrix2,1,sum)
  cnv_result2<-sort(cnv_result2,decreasing = T)
  cnv_result3<-as.data.frame(cnv_result2)
  cnv_result3$CancerType<-rownames(cnv_result3)
  colnames(cnv_result3)<-c("cnvNumber","CancerType")
  
  cnv_sum<-sum(cnv_result3$cnvNumber)
  
  cnv_res4<-data.frame(pathway=file_name[i],cnvSum=cnv_sum)
  
  cnv_df<-rbind(cnv_df,cnv_res4)
}

cnv_df<-cnv_df[order(cnv_df$cnvSum),]
write.csv(cnv_df,"通路免疫衰老基因的总拷贝数.csv",row.names = F)

#################不同癌症中47基因突变总数/拷贝总数排序柱状图

cnv_all<-read.csv("F:/immuno2.0/result2/泛癌-基因拷贝总数.csv")

mut_all<-read.csv("F:/immuno2.0/result2/泛癌-基因突变总数.csv")


############拷贝数柱状图

cnv_all$CancerType <- factor(cnv_all$CancerType, levels = cnv_all$CancerType)
head(cnv_all)

library(ggplot2)

ggplot(cnv_all, aes(CancerType, CNVnumber,fill="cancer")) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  
  #写label
  geom_text(data=cnv_all,
            mapping=aes(x=1:dim(cnv_all)[1], y= 0, label= paste0(" ", CancerType)),#bar跟坐标轴间留出间隙
            size = 2, #字的大小
            hjust = "outward" ) +  #字的对齐方式
  
  xlab("") +ylab("Number of immunosenescence-related CNV in Cancer")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴


############突变数柱状图

mut_all$CancerType <- factor(mut_all$CancerType, levels = mut_all$CancerType)
head(mut_all)

library(ggplot2)

ggplot(mut_all, aes(CancerType, MutNumber,fill="cancer")) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  
  #写label
  geom_text(data=mut_all,
            mapping=aes(x=1:dim(mut_all)[1], y= 0, label= paste0(" ", CancerType)),#bar跟坐标轴间留出间隙
            size = 2, #字的大小
            hjust = "outward" ) +  #字的对齐方式
  
  xlab("") +ylab("Number of immunosenescence-related mutation in Cancer")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 0.6)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴


#####Gse41613 -------------------------------------------------------------------------


rm(list=ls())

setwd("F:/immuno2.0/revise/GEO/GSE41613/")
library(dplyr)
library(data.table)

sur<-fread("survival.txt")%>%as.data.frame()
exp<-fread("GSE41613_series_matrix.txt")%>%as.data.frame()

library(stringr)

sur<-sur[,-1]

sample=sur[1,]
sex=sur[12,]
status=sur[13,]
time=sur[14,]

sur.df<-rbind(sample,sex,status,time)
colnames(sur.df)<-sur.df[1,]

sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","sex","status","time")
sur.df2<-as.data.frame(sur.df2)

sur.df2$sex<-sur.df2$sex %>%str_remove("Sex: ")%>%as.factor()
sur.df2$status<-ifelse(sur.df2$status=="vital: Alive",0,1)%>%as.numeric()
sur.df2$time<-sur.df2$time%>%str_remove("fu time: ")%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

sur.df3_f<-subset(sur.df3,sur.df3$sex=="F")
sur.df3_m<-subset(sur.df3,sur.df3$sex=="M")

sur.df4<-rbind(sur.df3_f,sur.df3_m[1:31,])

###性别差异生存分析

library(survminer)
library(ggpubr)
library(survival)
library(ggplot2)

fit<- survfit(Surv(time, status) ~ sex, data = sur.df3)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = c("#F8766D","#00BFC4"))

######GSE42743 -------------------------------------------------------------------------


rm(list=ls())

setwd("F:/immuno2.0/revise/GEO/GSE42743/")
library(dplyr)
library(data.table)

sur<-fread("survival.txt")%>%as.data.frame()

library(stringr)

sur<-sur[,-1]

sample=sur[1,]
sex=sur[11,]

sample<-c()
status<-c()
time<-c()

s<-c("survivallastfollowup: Dead with Disease",
     "survivallastfollowup: Living NED",
     "survivallastfollowup: Dead, unknown disease status",
     "survivallastfollowup: Died, due to treatment of disease",
     "survivallastfollowup: Dead, due to disease",
     "survivallastfollowup: Dead, other cause")

for(i in 1:dim(sur)[2]){
  
  x=sur[,i]
  sample[i]<-x[1]
  
  for(j in 1:length(x)){
    if(x[j]%in%s){
      status[i]<-x[j]
      time[i]<-x[j+3]
    }
  }
}

sur.df<-rbind(sample,sex,status,time)
colnames(sur.df)<-sur.df[1,]

sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","sex","status","time")
sur.df2<-as.data.frame(sur.df2)

sur.df2$sex<-sur.df2$sex %>%str_remove("gender: ")%>%as.factor()
sur.df2$status<-ifelse(sur.df2$status=="survivallastfollowup: Living NED",0,1)%>%as.numeric()
sur.df2$time<-sur.df2$time%>%str_remove("futime: ")%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

sur.df3_f<-subset(sur.df3,sur.df3$sex=="Female")
sur.df3_m<-subset(sur.df3,sur.df3$sex=="Male")

sur.df4<-rbind(sur.df3_f,sur.df3_m[1:20,])

###性别差异生存分析

library(survminer)
library(ggpubr)
library(survival)
library(ggplot2)

fit<- survfit(Surv(time, status) ~ sex, data = sur.df3)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = c("#F8766D","#00BFC4"))

# GSE65858 ----------------------------------------------------------
rm(list=ls())

setwd("F:/immuno2.0/revise/GEO/GSE65858/")
library(dplyr)
library(data.table)

sur<-fread("survival.txt")%>%as.data.frame()

library(stringr)

sur<-sur[,-1]

sample=sur[1,]
sex=sur[10,]
status=sur[25,]
time=sur[24,]

sur.df<-rbind(sample,sex,status,time)
colnames(sur.df)<-sur.df[1,]

sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","sex","status","time")
sur.df2<-as.data.frame(sur.df2)

sur.df2$sex<-sur.df2$sex %>%str_remove("gender: ")%>%as.factor()
sur.df2$status<-ifelse(sur.df2$status=="os_event: FALSE",0,1)%>%as.numeric()
sur.df2$time<-sur.df2$time%>%str_remove("os: ")%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

sur.df3_f<-subset(sur.df3,sur.df3$sex=="F")
sur.df3_m<-subset(sur.df3,sur.df3$sex=="M")

sur.df4<-rbind(sur.df3_f,sur.df3_m[1:47,])

###性别差异生存分析

library(survminer)
library(ggpubr)
library(survival)
library(ggplot2)

fit<- survfit(Surv(time, status) ~ sex, data = sur.df4)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = c("#F8766D","#00BFC4"))

#####表达-性别差异分析

exp<-fread("GSE65858_series_matrix.txt")%>%as.data.frame()
rownames(exp)<-exp$ID_REF
exp<-exp[,-1]

exp2<-exp[,sur.df4$sample]
#表达谱和样本性别对应n=74

probe<-fread("GPL10558-50081.txt")%>%as.data.frame()
probe2<-probe[,c("ID","ILMN_Gene")]
probe3<-subset(probe2,probe2$ILMN_Gene!="")

#探针和表达谱的基因对应
id<-intersect(rownames(exp2),probe3$ID)
exp3<-exp2[id,]
exp3$ID<-rownames(exp3)

exp_prob<-inner_join(probe3,exp3) #探针和表达谱的基因对应

###多个探针均值作为基因表达值
exp5<-aggregate(exp_prob,by=list(exp_prob$ILMN_Gene),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

##通路基因
library(clusterProfiler)
hsa16<-read.csv("f:/workplace/enrichKEGG16.csv",header=T) #16条免疫衰老相关通路
pathway<-hsa16$ID%>%as.character()  #16条通路的ID向量

hsa_gene<-read.csv("f:/workplace/HSA_KEGG.csv",header=T)   #KEGG通路和基因EntrzID对应列表

pathway16<-dplyr::filter(hsa_gene,hsa_gene$KEGGID%in%pathway)  #通路的基因
pathway_gene<-bitr(geneID=pathway16$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(pathway_gene)<-c("ENTREZID","gene") %>%unique() #通路的基因和ENTRZID

pathway_gene_result<-pathway_gene$gene%>%unique()

path_exp<-exp5[pathway_gene_result,]
path_exp<-distinct(path_exp)
path_exp2<-path_exp[rownames(path_exp)!="NA",]
#去掉NA值

###limma差异分析-sex group

sex<-data.frame(sample=sur.df4$sample,group=sur.df4$sex)

limma_input<-path_exp2[,sex$sample] 
#把limma的输入表达谱按性别排序【F,M】
identical(colnames(limma_input),sex$sample)

## 强制限定顺序
group<-ifelse(sex$group=="F","Female","Male")
group_list <- factor(group,levels = c("Female","Male"),ordered = F) 

library(limma)
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
fit <- lmFit(limma_input , design)
contrast.matrix <- makeContrasts(Female - Male,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG <- na.omit(tempOutput)  #免疫衰老基因的P值和logFC

dif<-subset(nrDEG,nrDEG$P.Value<0.5)#显著差异的结果

mygene<-c("NGF","CXCL9","IGF1R","PRKAA2","BCL2","BCL2L11","TAB1","SIRT1","CX3CL1","TP53","AIRE")

dif2<-dif[mygene,]%>%distinct()
dif3<-dif2[rownames(dif2)!="NA",]

####差异表达热图

library(pheatmap)

heat_input<-limma_input[rownames(dif3),] #性别差异显著的通路基因的表达
heat_input<-heat_input[c("CXCL9","CX3CL1","AIRE"),]

cormat<-round(cor(heat_input , method = "pearson"),2)

annotation_col<-data.frame(gender = group_list) #按正常疾病分组
rownames(annotation_col) = colnames(heat_input)

ann_colors = list(
  gender = c(Female = "#E83140", Male = "#3C7DAF")
)


heat_input[which(is.na(heat_input))]=0

n <- t(scale(t(heat_input)))
n[n > 2] <- 2
n[n < -2] <- -2
df <- n
rownames(df)<-rownames(heat_input)


library(RColorBrewer)  #配色包

pheatmap(df,cellwidth =2, cellheight = 20, fontsize = 14,fontsize_row=14,
         method="pearson", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=F,#为基因做聚类
         cluster_cols=F,#为sample做聚类
         color = colorRampPalette(c("#3C7DAF", "white", "#E83140"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         treeheight_col = "0",#不画树
         border_color = "NA")





# ######GSE117973 ---------------------------------------------------------


rm(list = ls())
setwd("F:/immuno2.0/revise/GEO/GSE117973/")
library(dplyr)
library(data.table)

sur<-fread("survival.txt")%>%as.data.frame()

sur<-sur[,-1]

sample=sur[1,]
sex=sur[10,]
status=sur[20,]
time=sur[21,]

sur.df<-rbind(sample,sex,status,time)
colnames(sur.df)<-sur.df[1,]

sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","sex","status","time")
sur.df2<-as.data.frame(sur.df2)

sur.df2$sex<-sur.df2$sex %>%str_remove("Sex: ")%>%as.factor()
sur.df2$status<-ifelse(sur.df2$status=="pfs event: 0",0,1)%>%as.numeric()
sur.df2$time<-sur.df2$time%>%str_remove("pfs months: ")%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

###性别差异生存分析

library(survminer)
library(ggpubr)
library(survival)
library(ggplot2)

fit<- survfit(Surv(time, status) ~ sex, data = sur.df3)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = c("#F8766D","#00BFC4"))


# ####E-MTAB-8588 ---------------------------------------------------------



setwd("F:/immuno2.0/revise/GEO/E-MTAB-8588/")

sur<-fread("1.txt")%>%as.data.frame()

sur.df<-data.frame(sample=sur$individual,
                   sex=sur$sex,
                   status=sur$`overall survival status`,
                   time=sur$`overall survival`)
sur.df2<-distinct(sur.df)

sur.df2$sex<-sur.df2$sex %>%as.factor()
sur.df2$status<-ifelse(sur.df2$status=="alive",0,1)%>%as.numeric()
sur.df2$time<-sur.df2$time%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]


fit<- survfit(Surv(time, status) ~ sex, data = sur.df3)

ggsurvplot(fit,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette = c("#F8766D","#00BFC4"))


# CCLE药物反应验证分析 -----------------------------------------------------------
rm(list = ls())

setwd("F:/workplace/结果5/ccle/")

drug<-fread("data_drug_treatment_IC50.txt")%>%as.data.frame()
cell<-fread("Cell_lines_annotations_20181226.txt")%>%as.data.frame()

HNSC<-subset(cell,cell$tcga_code=="HNSC")
HNSC2<-HNSC[,c("CCLE_ID","Gender")]
HNSC3<-subset(HNSC2,HNSC2$Gender%in%c("female","male"))

drug2<-drug[,intersect(colnames(drug),HNSC3$CCLE_ID)]
drug3<-cbind(drug[,c(1,2,4)],drug2)

osi027<-subset(drug3,drug3$NAME=="OSI-027")

osi027<-t(osi027)

osi_df<-data.frame(CCLE_ID=rownames(osi027),OSI027=osi027)
osi_df2<-osi_df[4:20,]

table(HNSC3$Gender) #F=4,M=29

input<-inner_join(osi_df2,HNSC3)
input$Gender<-as.factor(input$Gender)

input_m<-subset(input,input$Gender=="male")
input_f<-subset(input,input$Gender=="female")

aver_m<-as.numeric(input_m$X162)%>%mean()
aver_f<-as.numeric(input_f$X162)%>%mean()

bar_df<-data.frame(group=c("male","female")%>%as.factor(),IC50=c(aver_m,aver_f))


#IC50性别差异-HNSC-OSI-027 ----------------------------------------------



library(ggplot2)
library(ggpubr)

wilcox.test(as.numeric(input_m$X162),as.numeric(input_f$X162))
#p-value = 0.4706

ggplot(data=bar_df, aes(x=group, y=IC50,fill=group)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=IC50), vjust=1.6, color="white", size=3.5)+
  theme_minimal()

# pRRophetic --------------------------------------------------------------
library(data.table)
library(dplyr)

setwd("F:/immuno2.0/revise/drug/")

mydrug<-fread("通路-癌症-药物对应关联表.csv")%>%as.data.frame()

d<-mydrug$DRUG_NAME%>%unique()
####20种免疫衰老潜在药物


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.14')
BiocManager::install("sva",force = TRUE)

library(pRRophetic)
library(ggplot2)

##可以看到 Cancer Genome Project (CGP)里面的也是一千多个细胞系的2万个基因的表达量矩阵，关键是它药物超过了200种
data(cgp2016ExprRma) 
dim(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
length(unique(drugData2016$Drug.name))
unique(drugData2016$Drug.name)

###验证7个免疫衰老潜在药物在HNSC中的抗肿瘤效果

d7<-intersect(d,unique(drugData2016$Drug.name))

# 输入自己的表达谱矩阵 --------------------------------------------------------------

HNSC<-exp
rownames(HNSC)<-HNSC$V1
HNSC<-HNSC[,-1]

HNSC<-HNSC%>%as.matrix()

######预测7个药物的治疗反应

for(i in 1:7){
  
  predictedPtype <- pRRopheticPredict(testMatrix=HNSC, 
                                      drug=d7[i],
                                      tissueType = "all", 
                                      batchCorrect = "eb",
                                      selection=1,
                                      dataset = "cgp2016")
  
  
  df <- stack(list(NR = predictedPtype[((studyResponse == 'PGx_Responder = NR') & bortIndex)],
                   R = predictedPtype[((studyResponse == 'PGx_Responder = R') & bortIndex)]))
  head(df)
  
  write.csv(df,paste(d7[i],"group.csv"))
  
  pvalue<-t.test(predictedPtype[((studyResponse == "PGx_Responder = NR") & bortIndex)], 
                 predictedPtype[((studyResponse == "PGx_Responder = R") & bortIndex)], 
                 alternative="greater")
  
  ggplot(data = df,
         aes(y = values,
             x = ind))+
    geom_boxplot(alpha = 0.3,
                 fill = c('#e94753','#47a4e9'))+
    theme_bw()+
    ylab(paste('Predicted', d7[i], 'Sensitivity',sep=" ")) +
    xlab(paste('Clinical Response',round(pvalue$p.value,3),sep = ":"))
  
  ggsave(paste( i, ".pdf", sep=""))
}
