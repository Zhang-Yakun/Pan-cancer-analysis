setwd("f:/workplace/LUAD/LUAD_multidata_sex_diff/")

library(data.table)
library(dplyr)


# 蛋白质数据 -------------------------------------------------------------------

protein<-fread("f:/workplace/LUAD/LUAD_RPPA_RBN")%>%as.data.frame()
rownames(protein)<-protein$Sample_description
protein<-protein[,-1]


# 输入通路基因 ------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

hsa16<-read.csv("f:/workplace/enrichKEGG16.csv",header=T) #16条免疫衰老相关通路
pathway<-hsa16$ID%>%as.character()  #16条通路的ID向量

hsa_gene<-read.csv("f:/workplace/HSA_KEGG.csv",header=T)   #KEGG通路和基因EntrzID对应列表

pathway16<-dplyr::filter(hsa_gene,hsa_gene$KEGGID%in%pathway)  #通路的基因
pathway_gene<-bitr(geneID=pathway16$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(pathway_gene)<-c("ENTREZID","gene") %>%unique() #通路的基因和ENTRZID

pathway_gene_result<-pathway_gene$gene%>%unique()


# 免疫衰老蛋白质表达 ---------------------------------------------------------------

prot<-intersect(rownames(protein),pathway_gene_result)
im_prot<-protein[prot,]

#############LUAD样本的临床信息

sample_cli<-fread("f:/workplace/LUAD/LUAD_clinical_survival.csv")%>%as.data.frame()
#THCA样本临床和生存信息

# PSM的数据准备—结合表达和临床信息 ------------------------------------------------------

pathway_prot<-t(im_prot)%>%as.data.frame()
pathway_prot$sample<-rownames(pathway_prot)

path_prot<-inner_join(pathway_prot,sample_cli)
table(path_prot$sex)
rownames(path_prot)<-path_prot$sample

library(dplyr)
library(ggplot2)
library(MatchIt)

path_prot$histology_subtype<-as.factor(path_prot$histology_subtype)
path_prot$histology_subtype=ifelse(path_prot$histology_subtype== "Lung Acinar Adenocarcinoma",1,path_prot$histology_subtype)
#把组织学亚型的分类变成数字编号

path_prot$stage<-as.factor(path_prot$stage)
path_prot$stage=ifelse(path_prot$stage== "[Discrepancy]",1,path_prot$stage)
#把stage的分类变成数字编号

path_prot$Grade<-as.factor(path_prot$Grade)
path_prot$Grade=ifelse(path_prot$Grade== "G1",1,path_prot$Grade)

path_prot$smoking<-as.factor(path_prot$smoking)
path_prot$smoking<-ifelse(path_prot$smoking=="[Discrepancy]",0,path_prot$smoking)

path_prot %>% group_by(sex) %>% summarise(sample_number = n(),
                                                    age = mean(age_at_diagnosis),
                                                    histology_subtype = median(histology_subtype),
                                                    survival_time = mean(os.time),
                                                    stage=median(stage),
                                                    #Grade=median(Grade),
                                                    smoking = median(smoking)
                                                    #purity = median(purity)
                                          )#校正之前先大致看男女样本之间的协变量差异

path_prot<-path_prot %>% mutate(sex=ifelse(sex=="MALE",1,0)) 
#把性别变成1,0格式,男性[222]为1，女性[262]为0
mc_ps <- glm(sex ~ age_at_diagnosis   + stage +smoking+histology_subtype ,
             family = binomial(), data = path_prot) #广义线性模型
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

clinical_cov<-c("age_at_diagnosis","stage","smoking","histology_subtype") #挑选协变量
gene<-colnames(pathway_prot)[-dim(pathway_prot)[2]] #列出在THCA中表达的免疫衰老基因名【87个】
path_prot_nomiss <- path_prot %>%  # MatchIt包不允许缺失值，先去除NA值
  dplyr::select(sample,gene, sex, one_of(clinical_cov)) %>%
  na.omit() 

rownames(path_prot_nomiss)<-path_prot_nomiss$sample
path_prot_nomiss<-path_prot_nomiss[,-1]  # 缺失值处理之后剩537个样本

mod_match <- matchit(sex ~ age_at_diagnosis  +stage +smoking+histology_subtype,
                     method = "nearest", data = path_prot_nomiss) #最近邻法的倾向性得分匹配过程
dta_m <- match.data(mod_match) #倾向性得分匹配后的结果,有distance和weights
dim(dta_m)  #414*87 匹配后剩余414个样本

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
  fn_bal(dta_m, "histology_subtype"),
  fn_bal(dta_m, "smoking"),
  fn_bal(dta_m, "stage"),
  nrow = 2, widths = c(1, 0.8)
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

write.csv(nrDEG87,"F:/workplace/pathway multi data/LUAD_limma_pathway_protein.csv",row.names = T)

dif<-subset(nrDEG87,nrDEG87$P.Value<0.05)#limma找到的性别差异免疫衰老基因409个
write.csv(dif,"f:/workplace/LUAD/LUAD_multidata_sex_diff/LUAD_diff_prot.csv")
