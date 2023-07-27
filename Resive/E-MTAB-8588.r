
setwd("D:/workplace/mywork/revise/GEO/E-MTAB-8588/")
rm(list=ls())

library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

# 通路基因的表达谱 ----------------------------------------------------------------

exp<-fread("normalized_matrix.txt")%>%as.data.frame()
rownames(exp)<-exp$V1
exp2<-exp[,-1]

probe<-fread("GPL6947-13512.txt")%>%as.data.frame()
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
expMatrix_tpm<-exp5[,-1]

# 通路基因 --------------------------------------------------------------------

pathway_gene<-read.table("D:/workplace/mywork/revise/GEO/pathway_gene.txt")

###########通路基因的表达数据

pathway_exp<-expMatrix_tpm[intersect(pathway_gene$x,rownames(expMatrix_tpm)),] #1216个基因，576个样本


# 临床数据 --------------------------------------------------------------------

sur<-fread("1.txt")%>%as.data.frame()

sur.df<-data.frame(sample=sur$`Source Name`,
                   sex=sur$sex,
                   status=sur$`overall survival status`,
                   time=sur$`overall survival`,
                   cancer=sur$`sampling site`,
                   age=sur$age,
                   stage=sur$`disease staging`,
                   grade=sur$`tumor grading`,
                   smoking=sur$`smoking behavior`,
                   alcohol=sur$`alcohol drinking`)
sur.df<-distinct(sur.df)

sur.df2<-subset(sur.df,sur.df$cancer=="tumour")
sur.df2<-sur.df2[,-5]

sur.df2$sex<-sur.df2$sex %>%as.factor()
sur.df2$status<-ifelse(sur.df2$status=="alive",0,1)%>%as.numeric()
sur.df2$time<-sur.df2$time%>%as.numeric()

sur.df2$stage<-ifelse(sur.df2$stage=="I",1,
                      ifelse(sur.df2$stage=="II",2,
                      ifelse(sur.df2$stage=="III",3,
                      ifelse(sur.df2$stage=="IV",4,sur.df2$stage))))%>%as.numeric()

sur.df2$grade<-sur.df2$grade%>%str_remove("G")%>%as.numeric()
sur.df2$smoking<-ifelse(sur.df2$smoking=="no",0,
                        ifelse(sur.df2$smoking=="yes",1,sur.df2$smoking))%>%as.numeric()
sur.df2$alcohol<-ifelse(sur.df2$alcohol=="no",0,
                        ifelse(sur.df2$alcohol=="yes",1,sur.df2$alcohol))%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]


# PSM的数据准备—结合表达和临床信息 ------------------------------------------------------

library(dplyr)
library(ggplot2)
library(MatchIt)
library(stringr)

exp<-t(pathway_exp)
exp2<-as.data.frame(exp)

rownames(exp2)<-rownames(exp)
colnames(exp2)<-colnames(exp)

exp3<-exp2%>%mutate(sample=rownames(exp2))%>%dplyr::select(sample,everything())

expression_clinical<-inner_join(exp3,sur.df3)  #将表达谱和临床表型信息整合在一个数据框中[484,89]

expression_clinical %>% group_by(sex) %>% summarise(sample_number = n(),
                                                    age = mean(age),
                                                    survival_time = mean(time),
                                                    smoking=median(smoking),
                                                    alcohol=median(alcohol))
#校正之前先大致看男女样本之间的协变量差异

expression_clinical<-expression_clinical %>% mutate(sex=ifelse(sex=="male",1,0)) 
#把性别变成1,0格式,男性[60]为1，女性[17]为0
mc_ps <- glm(sex ~ age +stage +grade + smoking +alcohol,
             family = binomial(), data = expression_clinical) #广义线性模型
summary(mc_ps)

pre_mc_ps<- data.frame(pr_score = predict(mc_ps, type = "response"),
                       sex = mc_ps$model$sex) #用模型预测得分作为倾向性得分
head(pre_mc_ps)

library(ggplot2)
pre_mc_ps %>%
  ggplot(aes(x = pr_score)) +
  geom_histogram(color = "white") +
  facet_wrap(~sex) +
  xlab("Probability") +
  theme_bw() #柱状图展示男女样本PS得分的分布

# 对倾向性得分进行匹配算法 ------------------------------------------------------------

clinical_cov<-c("age","stage","grade","smoking","alcohol") #挑选协变量
gene<-colnames(exp3)[-1] #列出在HNSC中表达的免疫衰老基因名
expression_clinical_nomiss <- expression_clinical %>%  # MatchIt包不允许缺失值，先去除NA值
  dplyr::select(sample,gene, sex, one_of(clinical_cov)) %>%
  na.omit() 

rownames(expression_clinical_nomiss)<-expression_clinical_nomiss$sample
expression_clinical_nomiss<-expression_clinical_nomiss[,-1]  # 缺失值处理之后

mod_match <- matchit(sex ~ age +stage +grade + smoking +alcohol,
                     method = "nearest", data = expression_clinical_nomiss) #最近邻法的倾向性得分匹配过程
summary(mod_match)
plot(mod_match,type="hist")

dta_m <- match.data(mod_match) #倾向性得分匹配后的结果,有distance和weights
dim(dta_m)  #匹配后剩余14个样本

write.csv(dta_m,"E-MTAB-8588_PSM.csv")

#检查匹配样本中协变量的平衡

gender=as.factor(ifelse(dta_m$sex==0,"Female","Male"))

fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$sex <- as.factor(dta$sex)
  support <- c(min(dta$variable), max(dta$variable))
  ggplot(dta, aes(x = distance, y = variable, color = gender)) +
    geom_point(alpha = 0.2, size = 1.3) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw() +
    ylim(support)
} #自己写的函数，用来画协变量匹配后的倾向得分（distance）分布

dev.off()
library(gridExtra) #用来拼图，类似plot_grid()

grid.arrange(
  fn_bal(dta_m, "age") + theme(legend.position = "none"),
  fn_bal(dta_m, "stage")+ theme(legend.position = "none"),
  fn_bal(dta_m, "grade")+ theme(legend.position = "none"),
  fn_bal(dta_m, "smoking")+ theme(legend.position = "none"),
  fn_bal(dta_m, "alcohol")+ theme(legend.position = "none"),
  nrow = 3, widths = c(1, 0.8)
) #展示4种协变量匹配后的两组倾向得分分布

dta_m %>%
  group_by(sex) %>%
  dplyr::select(one_of(clinical_cov)) %>%
  summarise_all(funs(mean))  #psm之后两组样本协变量的均值差异

lapply(clinical_cov, function(v) {
  t.test(dta_m[, v] ~ dta_m$sex)
})  #T检验检查两组样本协变量的差异，
####βx代表匹配后样本中治疗组和对照组协变量均值之间的差异。
####平均绝对标准化差值最好接近0，因为这表明在匹配的样本中，对照组和治疗组之间存在小的差异。

# 协变量平衡之后就可以用limma评估两组样本表达值的差异了 ----------------------------------------

gene90<-read.table("D:/workplace/mywork/revise/GEO/gene90.txt")

g<-intersect(unique(gene90$V1),colnames(dta_m))

dta_m2<-dta_m[,g] #94gene expression

limma_input<-t(dta_m2)
limma_input<-as.data.frame(limma_input)  #行是基因（1216），列是样本（414）的表达谱，之后作为limma输入做差异

sex<-dta_m$sex
limma_input2<-limma_input[,order(sex)] #把limma的输入表达谱按性别排序【0,1】


group_list=c(rep('Female',dim(limma_input2)[2]/2),rep('Male',dim(limma_input2)[2]/2))
## 强制限定顺序
group_list <- factor(group_list,levels = c("Female","Male"),ordered = F) 

library(limma)
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
fit <- lmFit(limma_input2 , design)
contrast.matrix <- makeContrasts(Female - Male,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG <- na.omit(tempOutput)  #免疫衰老基因的P值和logFC，作为GSEA的输入

dif<-subset(nrDEG,nrDEG$P.Value<0.5)#limma找到的性别差异免疫衰老基因409个

dif<-dif[order(dif$logFC),]
diff_gene<-rownames(dif)

# 差异基因的表达热图 ---------------------------------------------------------------

library(pheatmap)

heat_input<-limma_input2[diff_gene,] #性别差异显著的通路基因的表达（414样本）

cormat<-round(cor(heat_input , method = "pearson"),2)

heat_sample<-data.frame(sample=colnames(heat_input))
heat_cov<-left_join(heat_sample,sur.df3[,c(1:2,5:9)])  #热图中的临床信息
heat_cov$age=ifelse(heat_cov$age>=50,">=50","<50")

annotation_col<-data.frame(gender = factor(rep(c("Female","Male"),c(dim(limma_input2)[2]/2,dim(limma_input2)[2]/2))),
                           age = as.factor(heat_cov$age),
                           stage = as.factor(heat_cov$stage),
                           grade = as.factor(heat_cov$grade),
                           smoking = as.factor(heat_cov$smoking),
                           alcohol=as.factor(heat_cov$alcohol)) #按正常疾病分组
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

pheatmap(df,cellwidth =10, cellheight = 6, fontsize = 6,fontsize_row=4,
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


# 性别差异分析 ------------------------------------------------------------------

###性别差异生存分析

library(survminer)
library(ggpubr)
library(survival)
library(ggplot2)

rownames(sur.df3)<-sur.df3$sample
sur.df4<-sur.df3[colnames(heat_input),]

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


# oncoPredict -------------------------------------------------------------
#比较正常和癌旁组织对药物敏感性差异

library(oncoPredict)
library(gtools)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggpubr)

###表达谱和tumour-normal分组
expMatrix_tpm[1:5,1:5]

sample_df<-data.frame(sur.df$sample,sur.df$cancer)
sample_df<-sample_df[order(sample_df$sur.df.cancer),]

exp_sam<-expMatrix_tpm[,sample_df$sur.df.sample]
identical(colnames(exp_sam),sample_df$sur.df.sample)
onco_input<-exp_sam%>%as.matrix()

###预测药物敏感性
#以GDSC2-训练集的表达量矩阵和药物处理信息
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

GDSC2_Expr = readRDS(file="GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res = readRDS(file = "GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 

#预测药物反应

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = onco_input,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )


###分组比较药反差异
drug<-read.csv("D:/workplace/mywork/revise/GEO/E-MTAB-8588/calcPhenotype_Output/DrugPredictions.csv")
library(dplyr)
library(data.table)
rownames(drug)<-drug$X
drug<-drug[,-1]

drug$group<-sample_df$sur.df.cancer%>%as.factor()

drug2<-drug%>%dplyr::select(group,everything())

library(stringr)
colnames(drug2)<-str_remove(colnames(drug2),"_[0-9]*")
colnames(drug2)<-str_replace(colnames(drug2),"[._]","-")
###

pt<-c()

for(i in 2:dim(drug2)[2]){
  
  f<-drug2[,i]
  
  t<-t.test(f ~ drug2$group)
  pt[i-1]=t$p.value
  names(pt)[i-1]<-colnames(drug2)[i]
  
}

pt_sig<-pt[pt<0.05]

names(pt_sig)

mydrug<-fread("D:/workplace/mywork/revise/通路-癌症-药物对应关联表.csv")%>%as.data.frame()

d<-mydrug$DRUG_NAME%>%unique()

drug6<-intersect(names(pt_sig),d)

####20种免疫衰老潜在药物

pt_sig<-pt_sig[drug6]

######显著差异的药物画IC50-boxplot

setwd("D:/workplace/mywork/revise/")

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)

plot_list<-list()
for(i in 1:length(pt_sig)){
  
  df<-drug2[,names(pt_sig)[i]]%>%as.data.frame()
  colnames(df)<-names(pt_sig)[i]
  
  df2<-cbind(drug2[,1],df)
  colnames(df2)[1]<-"group"
  colnames(df2)[2]<-"drug"
  plot_list[[i]]<-ggboxplot(df2, x = "group", y = "drug" , color = "group",add.params = list(color = "group",size=1),add = "jitter",
                            palette = c("#1F78B4","#A65628"),title=names(pt_sig)[i])+stat_compare_means(method = "t.test")+ylab("Estimated IC50")+xlab(names(pt_sig)[i])+theme(legend.position = "none")
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[5]],plot_list[[6]])
