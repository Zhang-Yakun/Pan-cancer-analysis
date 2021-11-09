
library(dplyr)
library(tibble)
library(data.table)
library(stringr)

# 服务器上：免疫衰老通路基因的甲基化谱 -----------------------------------------------------------

methy<-read.table("/pub2/zhoudianshuang/Zhangyakun/methyl/LUAD_HumanMethylation450.txt",header=T)
colnames(methy)[1]="probe"

pathway_gene<-read.table("/pub2/zhoudianshuang/Zhangyakun/methyl/pathway_gene.txt")

probe <- read.table(file="/pub2/zhoudianshuang/Zhangyakun/methyl/illuminaMethyl450_hg19_GPL16304_TCGAlegacy.txt",sep = "\t")

in_gene<-intersect(probe$V2,pathway_gene$V1)
pathway_probe<-subset(probe,probe$V2%in%in_gene)

pb<-data.frame(probe=pathway_probe$V1,gene=pathway_probe$V2 )

gene_methy<-inner_join(pb,methy)

write.csv(gene_methy,"/pub2/zhoudianshuang/Zhangyakun/methyl/LUAD_pathway_methylation.csv",row.names=F)


# 结合临床信息，统计男女癌症样本量 ---------------------------------------------------

setwd("f:/workplace/pathway multi data/")

methyl<-read.csv("LUAD_pathway_methylation.csv")
colnames(methyl)<-str_replace_all(colnames(methyl),pattern = "[.*]",replacement = "-")

rownames(methyl)<-methyl$probe
methyl<-methyl[,-1]

sample_cli<-fread("f:/workplace/LUAD/LUAD_clinical_survival.csv")%>%as.data.frame()

in_sam<-intersect(colnames(methyl)[-1],sample_cli$sample)
methy_cli<-subset(sample_cli,sample_cli$sample%in%in_sam)
methy_cli$group<-ifelse(as.numeric(str_sub(methy_cli$sample,14,15))<10,"Tumor","Normal")

methy_cli<-subset(methy_cli,methy_cli$group=="Tumor")

table(methy_cli$sex)



# 结合基因表达谱，将甲基化探针和基因唯一对应 ---------------------------------------------------

######表达谱

LUAD_tpm<-fread("f:/workplace/LUAD/LUAD_readcount.genes.tpm.txt")%>%as.data.frame()
rownames(LUAD_tpm)<-LUAD_tpm$V1
LUAD_tpm<-LUAD_tpm[,-1]

expMatrix_tpm<-(log2(LUAD_tpm+1))  #所有基因在LUAD的表达log2+1(TPM值)

inter_gene<-intersect(rownames(expMatrix_tpm),methyl$gene)
expr<-expMatrix_tpm[inter_gene,]

inter_sample<-intersect(colnames(expr),colnames(methyl))

expr_match<-expr[,inter_sample]
methy_match<-methyl[,inter_sample]

identical(colnames(expr_match),colnames(methy_match))  #表达谱和甲基化谱的基因和样本相对应
methy_match$gene<-methyl$gene

######结合基因表达，将甲基化基因和探针一一对应，取最负相关的探针

df=data.frame()
pg_match<-data.frame()  #甲基化探针和基因名一一对应

for(i in 1:dim(expr_match)[1]){
  pro<-subset(methy_match,methy_match$gene==rownames(expr_match)[i])
  pr<-pro[,-dim(pro)[2]]
  
  for(j in 1:length(rownames(pr))){
    r=cor.test(as.numeric(pr[j,]),as.numeric(expr_match[i,]) )
    df2=data.frame(r=r$estimate,p=rownames(pr)[j])
    df=rbind(df2,df)
  }
  min_r<-min(df$r)
  if(!is.na(min_r)){
    final_probe<-subset(df,df$r==min_r)$p%>%as.character()
    df3<-data.frame(probe=final_probe,gene=rownames(expr_match)[i])
    pg_match<-rbind(pg_match,df3)
  }
}

########甲基化谱

setwd("f:/workplace/pathway multi data/")

methyl<-read.csv("LUAD_pathway_methylation.csv")
colnames(methyl)<-str_replace_all(colnames(methyl),pattern = "[.*]",replacement = "-")

rownames(methyl)<-methyl$probe
methyl<-methyl[,-1]

sample_cli<-fread("f:/workplace/LUAD/LUAD_clinical_survival.csv")%>%as.data.frame()

in_sam<-intersect(colnames(methyl)[-1],sample_cli$sample)
methy2<-methyl[,in_sam]

methy2$probe<-rownames(methy2)
methy3<-inner_join(pg_match,methy2)
LUAD_methy<-methy3[,-1]
rownames(LUAD_methy)<-LUAD_methy$gene

LUAD_methy<-LUAD_methy[,-1]
#探针转化为基因的一一对应甲基化谱
write.csv(LUAD_methy,"f:/workplace/pathway multi data/LUAD_unique_methyl_matrix.csv")

# PSM的数据准备—结合表达和临床信息 ------------------------------------------------------

library(dplyr)
library(ggplot2)
library(MatchIt)

met<-t(LUAD_methy)
met2<-as.data.frame(met)

rownames(met2)<-rownames(met)
colnames(met2)<-colnames(met)

met3<-met2%>%mutate(sample=rownames(met2))%>%dplyr::select(sample,everything())

methylation_clinical<-inner_join(met3,sample_cli)  #将表达谱和临床表型信息整合在一个数据框中

methylation_clinical$group<-ifelse(as.numeric(str_sub(methylation_clinical$sample,14,15))<10,"Tumor","Normal")

methylation_clinical<-subset(methylation_clinical ,methylation_clinical$group=="Tumor") #剔除正常样本

methylation_clinical$histology_subtype<-as.factor(methylation_clinical$histology_subtype)
methylation_clinical$histology_subtype=ifelse(methylation_clinical$histology_subtype== "Lung Adenocarcinoma- Not Otherwise Specified (NOS)",1,methylation_clinical$histology_subtype)
#把组织学亚型的分类变成数字编号

methylation_clinical$stage<-as.factor(methylation_clinical$stage)
methylation_clinical$stage=ifelse(methylation_clinical$stage== "[Discrepancy]",1,methylation_clinical$stage)
#把stage的分类变成数字编号

methylation_clinical$Grade<-as.factor(methylation_clinical$Grade)
methylation_clinical$Grade=ifelse(methylation_clinical$Grade== "G1",1,methylation_clinical$Grade)

methylation_clinical$smoking<-as.factor(methylation_clinical$smoking)
methylation_clinical$smoking<-ifelse(methylation_clinical$smoking=="NA",0,methylation_clinical$smoking)

methylation_clinical %>% group_by(sex) %>% summarise(sample_nuOSmber = n(),
                                                     age = mean(age_at_diagnosis),
                                                     histology_subtype = median(histology_subtype),
                                                     survival_time = mean(OS.time),
                                                     stage=median(stage),
                                                     #Grade=median(Grade),
                                                     smoking = median(smoking),
                                                     purity = median(purity)
)#校正之前先大致看男女样本之间的协变量差异

methylation_clinical<-methylation_clinical %>% mutate(sex=ifelse(sex=="MALE",1,0)) 
#把性别变成1,0格式,男性[222]为1，女性[262]为0
mc_ps <- glm(sex ~ age_at_diagnosis  + histology_subtype + stage + smoking+purity,
             family = binomial(), data = methylation_clinical) #广义线性模型
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

clinical_cov<-c("age_at_diagnosis","histology_subtype","stage","smoking","purity") #挑选协变量
gene<-colnames(met3)[-1] #列出在LUAD中表达的通路基因名
methylation_clinical_nomiss <- methylation_clinical %>%  # MatchIt包不允许缺失值，先去除NA值
  dplyr::select(sample,gene, sex, one_of(clinical_cov)) %>%
  na.omit() 

rownames(methylation_clinical_nomiss)<-methylation_clinical_nomiss$sample
methylation_clinical_nomiss<-methylation_clinical_nomiss[,-1]  # 缺失值处理之后剩的样本

mod_match <- matchit(sex ~ age_at_diagnosis + histology_subtype +stage +smoking+purity,
                     method = "nearest", data = methylation_clinical_nomiss) #最近邻法的倾向性得分匹配过程
dta_m <- match.data(mod_match) #倾向性得分匹配后的结果,有distance和weights
dim(dta_m)  #匹配后剩余的样本

write.csv(dta_m,"f:/workplace/LUAD/LUAD_multidata_sex_diff/LUAD_甲基化匹配后的样本.csv")

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

dta_m2<-dta_m[,1:length(gene)] # pathway gene expression

limma_input<-t(dta_m2)
limma_input<-as.data.frame(limma_input)  #行是基因，列是样本）的表达谱，之后作为limma输入做差异

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
nrDEG87 <- na.omit(tempOutput)  #通路基因的P值和logFC，作为GSEA的输入

write.csv(nrDEG87,"F:/workplace/pathway multi data/LUAD_limma_pathway_methylation.csv",row.names = T)

dif<-subset(nrDEG87,nrDEG87$P.Value<0.05)#limma找到的性别差异通路基因
write.csv(dif,"f:/workplace/LUAD/LUAD_multidata_sex_diff/LUAD_diff_methyl.csv")


