setwd("f:/workplace/HNSC/HNSC_multidata_sex_diff/")

library(data.table)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 通路基因的表达谱 ----------------------------------------------------------------


HNSC_tpm<-fread("f:/workplace/HNSC/HNSC_readcount.genes.tpm.txt")%>%as.data.frame()
rownames(HNSC_tpm)<-HNSC_tpm$V1
HNSC_tpm<-HNSC_tpm[,-1]

expMatrix_tpm<-(log2(HNSC_tpm+1))  #所有基因在HNSC的表达log2+1(TPM值)

# 通路基因 --------------------------------------------------------------------

hsa16<-read.csv("f:/workplace/enrichKEGG16.csv",header=T) #16条免疫衰老相关通路
pathway<-hsa16$ID%>%as.character()  #16条通路的ID向量

hsa_gene<-read.csv("f:/workplace/HSA_KEGG.csv",header=T)   #KEGG通路和基因EntrzID对应列表

pathway16<-dplyr::filter(hsa_gene,hsa_gene$KEGGID%in%pathway)  #通路的基因
pathway_gene<-bitr(geneID=pathway16$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(pathway_gene)<-c("ENTREZID","gene") %>%unique() #通路的基因和ENTRZID

pathway_gene_result<-pathway_gene$gene%>%unique()

#write.table(pathway_gene_result,"f:/workplace/pathway_gene.txt",row.names = F,col.names = F,quote = F)

###########通路基因的表达数据

pathway_exp<-expMatrix_tpm[intersect(pathway_gene$gene,rownames(expMatrix_tpm)),] #1216个基因，576个样本


#############HNSC样本的临床信息

sample_cli<-fread("f:/workplace/HNSC/HNSC_clinical_survival.csv")%>%as.data.frame()
#HNSC样本临床和生存信息

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

expression_clinical<-inner_join(exp3,sample_cli)  #将表达谱和临床表型信息整合在一个数据框中[484,89]

expression_clinical$group<-ifelse(as.numeric(str_sub(expression_clinical$sample,14,15))<10,"Tumor","Normal")

expression_clinical<-subset(expression_clinical ,expression_clinical$group=="Tumor") #剔除正常样本

expression_clinical$histology_subtype<-as.factor(expression_clinical$histology_subtype)
expression_clinical$histology_subtype=ifelse(expression_clinical$histology_subtype== "Head & Neck Squamous Cell Carcinoma",1,expression_clinical$histology_subtype)
#把组织学亚型的分类变成数字编号

expression_clinical$stage<-as.factor(expression_clinical$stage)
expression_clinical$stage=ifelse(expression_clinical$stage== "[Discrepancy]",1,expression_clinical$stage)
#把stage的分类变成数字编号

expression_clinical$Grade<-as.factor(expression_clinical$Grade)
expression_clinical$Grade=ifelse(expression_clinical$Grade== "G1",1,expression_clinical$Grade)

expression_clinical$smoking<-as.factor(expression_clinical$smoking)
expression_clinical$smoking<-ifelse(expression_clinical$smoking=="NA",0,expression_clinical$smoking)

expression_clinical %>% group_by(sex) %>% summarise(sample_number = n(),
                                                    age = mean(age_at_diagnosis),
                                                    histology_subtype = median(histology_subtype),
                                                    survival_time = mean(OS.time),
                                                    stage=median(stage),
                                                    Grade=median(Grade),
                                                    smoking = median(smoking)
                                                    #purity = median(purity)
                                                    )#校正之前先大致看男女样本之间的协变量差异

expression_clinical<-expression_clinical %>% mutate(sex=ifelse(sex=="MALE",1,0)) 
#把性别变成1,0格式,男性[222]为1，女性[262]为0
mc_ps <- glm(sex ~ age_at_diagnosis  + histology_subtype + stage +Grade +smoking ,
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

clinical_cov<-c("age_at_diagnosis","histology_subtype","stage","Grade","smoking") #挑选协变量
gene<-colnames(exp3)[-1] #列出在HNSC中表达的免疫衰老基因名【87个】
expression_clinical_nomiss <- expression_clinical %>%  # MatchIt包不允许缺失值，先去除NA值
  dplyr::select(sample,gene, sex, one_of(clinical_cov)) %>%
  na.omit() 

rownames(expression_clinical_nomiss)<-expression_clinical_nomiss$sample
expression_clinical_nomiss<-expression_clinical_nomiss[,-1]  # 缺失值处理之后剩537个样本

mod_match <- matchit(sex ~ age_at_diagnosis + histology_subtype +stage+Grade+smoking ,
                     method = "nearest", data = expression_clinical_nomiss) #最近邻法的倾向性得分匹配过程
summary(mod_match)
plot(mod_match,type="hist")

dta_m <- match.data(mod_match) #倾向性得分匹配后的结果,有distance和weights
dim(dta_m)  #414*87 匹配后剩余414个样本

write.csv(dta_m,"f:/workplace/HNSC/HNSC_multidata_sex_diff/HNSC_表达匹配后的样本.csv")
# 检查匹配样本中协变量的平衡 -----------------------------------------------------------

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

library(gridExtra) #用来拼图，类似plot_grid()


grid.arrange(
  fn_bal(dta_m, "age_at_diagnosis") + theme(legend.position = "none"),
  fn_bal(dta_m, "histology_subtype")+ theme(legend.position = "none"),
  fn_bal(dta_m, "Grade")+ theme(legend.position = "none"),
  fn_bal(dta_m, "smoking")+ theme(legend.position = "none"),
  fn_bal(dta_m, "stage")+ theme(legend.position = "none"),
  fn_bal(dta_m, "stage"),
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

dta_m2<-dta_m[,1:length(gene)] #87gene expression

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
nrDEG87 <- na.omit(tempOutput)  #87个免疫衰老基因的P值和logFC，作为GSEA的输入

write.csv(nrDEG87,"F:/workplace/pathway multi data/HNSC_limma_pathway_expression.csv",row.names = T)

dif<-subset(nrDEG87,nrDEG87$P.Value<0.05)#limma找到的性别差异免疫衰老基因409个
write.csv(dif,"f:/workplace/HNSC/HNSC_multidata_sex_diff/HNSC_diff_exp.csv")

diff_gene<-subset(nrDEG87,nrDEG87$P.Value<0.0005)%>%rownames()

# 差异基因的表达热图 ---------------------------------------------------------------

library(pheatmap)

heat_input<-limma_input2[diff_gene,] #性别差异显著的通路基因的表达（414样本）

cormat<-round(cor(heat_input , method = "pearson"),2)

heat_sample<-data.frame(sample=colnames(heat_input))
heat_cov<-left_join(heat_sample,sample_cli[,1:6])  #热图中的临床信息
heat_cov$age=ifelse(heat_cov$age_at_diagnosis>=50,">=50","<50")

annotation_col<-data.frame(gender = factor(rep(c("Female","Male"),c(dim(limma_input2)[2]/2,dim(limma_input2)[2]/2))),
                           age = as.factor(heat_cov$age),
                           #histology_subtype = as.factor(heat_cov$histology_subtype),
                           stage = as.factor(heat_cov$stage),
                           smoking = as.factor(heat_cov$smoking)) #按正常疾病分组
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

pheatmap(df,cellwidth =0.5, cellheight = 4, fontsize = 6,fontsize_row=4,
         method="pearson", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类
         cluster_cols=F,#为sample做聚类
         color = colorRampPalette(c("#3C7DAF", "white", "#E83140"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         treeheight_col = "0",#不画树
         border_color = "NA")


