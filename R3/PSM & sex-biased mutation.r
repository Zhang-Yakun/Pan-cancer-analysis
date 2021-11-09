setwd("f:/workplace/KIRP/KIRP_multidata_sex_diff/")

library(dplyr)
library(data.table)
library(tibble)
library(clusterProfiler) #转换基因ID的包
library(org.Hs.eg.db)
library(stringr)


# 挑出通路基因的突变数据 -------------------------------------------------------------

sample_cli<-fread("f:/workplace/KIRP/KIRP_clinical_survival.csv")%>%as.data.frame()
#KIRP样本临床和生存信息

mutation_01<-fread("f:/workplace/KIRP/KIRP_mutation_preprocess_01.csv")%>%as.data.frame()
#预处理后的突变0-1谱
rownames(mutation_01)<-mutation_01$V1
mutation_01<-mutation_01[,-1]

mutation<-fread("f:/workplace/KIRP/KIRP_mutation_preprocess.csv")%>%as.data.frame()
#预处理后的突变数据

# 通路基因 --------------------------------------------------------------------

hsa16<-read.csv("f:/workplace/enrichKEGG16.csv",header=T) #16条免疫衰老相关通路
pathway<-hsa16$ID%>%as.character()  #16条通路的ID向量

hsa_gene<-read.csv("f:/workplace/HSA_KEGG.csv",header=T)   #KEGG通路和基因EntrzID对应列表

pathway16<-dplyr::filter(hsa_gene,hsa_gene$KEGGID==pathway)  #16条通路的基因
pathway_gene16<-bitr(geneID=pathway16$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(pathway_gene16)<-c("ENTREZID","gene") %>%unique() #通路的基因和ENTRZID


m_g<-colnames(mutation_01)%>%intersect(pathway_gene16$gene) #通路基因在KIRP中发生突变0-1谱
path_mut_01<-mutation_01[,m_g]
rownames(path_mut_01)<-rownames(mutation_01)
#60基因，488样本突变0-1谱


# PS方法找突变在男女之间有差异的基因 ----------------------------------------------------------

library(ggplot2)
library(MatchIt)

m_488<-as.matrix.data.frame(path_mut_01)
rownames(m_488)<-rownames(path_mut_01)
colnames(m_488)<-colnames(path_mut_01)

mutation_488<-as.data.frame(m_488) #把table格式的突变普变成数据框格式
m488<-mutation_488%>%mutate(sample=rownames(mutation_488))%>%dplyr::select(sample,everything())
mutation_clinical<-inner_join(m488,sample_cli)  #将突变谱和临床表型信息整合在一个数据框中[488,19]

mutation_clinical$histology_subtype<-as.factor(mutation_clinical$histology_subtype)
mutation_clinical$histology_subtype=ifelse(mutation_clinical$histology_subtype== "Kidney Papillary Renal Cell Carcinoma",1,mutation_clinical$histology_subtype)
#把组织学亚型的分类变成数字编号

mutation_clinical$stage<-as.factor(mutation_clinical$stage)
mutation_clinical$stage=ifelse(mutation_clinical$stage== "[Discrepancy]",1,mutation_clinical$stage)
#把stage的分类变成数字编号

mutation_clinical$Grade<-as.factor(mutation_clinical$Grade)
mutation_clinical$Grade=ifelse(mutation_clinical$Grade== "G1",1,mutation_clinical$Grade)


mutation_clinical %>% group_by(sex) %>% summarise(sample_number = n(),
                                                  age = mean(age_at_diagnosis),
                                                  #histology_subtype = median(histology_subtype),
                                                  survival_time = mean(OS.time),
                                                  stage=median(stage),
                                                  #Grade=median(Grade),
                                                  smoking = median(smoking)
                                                  #purity = median(purity)
                                                  )#校正之前先大致看男女样本之间的协变量差异

mutation_clinical<-mutation_clinical %>% mutate(sex=ifelse(sex=="MALE",1,0)) 
#把性别变成1,0格式,男性[225]为1，女性[263]为0
mc_ps <- glm(sex ~ age_at_diagnosis +  stage+smoking,
             family = binomial(), data = mutation_clinical) #广义线性模型
summary(mc_ps)

pre_mc_ps<- data.frame(pr_score = predict(mc_ps, type = "response"),
                       sex = mc_ps$model$sex) #用模型预测得分作为倾向性得分
head(pre_mc_ps)
pre_mc_ps %>%
  ggplot(aes(x = pr_score)) +
  geom_histogram(color = "white") +
  facet_wrap(~sex) +
  xlab("Probability") +
  theme_bw() #柱状图展示男女样本PS得分的分布

# 对倾向性得分进行匹配算法 ------------------------------------------------------------

clinical_cov<-c("age_at_diagnosis","stage","smoking") #挑选协变量
gene<-colnames(m488)[-1] #列出基因名【19个】
mutation_clinical_nomiss <- mutation_clinical %>%  # MatchIt包不允许缺失值，先去除NA值
  dplyr::select(gene, sex,sample, one_of(clinical_cov)) %>%
  na.omit() 

mod_match <- matchit(sex ~ age_at_diagnosis+stage+smoking,
                     method = "nearest", data = mutation_clinical_nomiss) #最近邻法的倾向性得分匹配过程
dta_m <- match.data(mod_match) #倾向性得分匹配后的结果,有distance和weights
dim(dta_m)  #420*68 匹配后剩余420个样本

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
  fn_bal(dta_m, "Grade"),
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

# 协变量平衡之后就可以评估两组样本表达值的差异了 ----------------------------------------

diff<-lapply(gene, function(v) {
  wilcox.test(dta_m[, v] ~ dta_m$sex)
})  #秩和检验检查两组样本基因表达值的差异

diff_p<-c()
for(i in 1:length(gene)){
  
  diff_p<-c(diff_p,diff[[i]][3]%>%as.numeric())
  
}

diff_gene<-gene[diff_p<0.05] #男女样本间突变差异的4个基因【"ADCY1" "ADCY5" 】

result<-data.frame(Gene=gene,wilcox.Pvalue=diff_p)

write.csv(result,"KIRP_diff_mut.csv")


# 把差异基因的突变信息输出，用mutationmapper展示 ------------------------------------------


mm<-subset(mutation,mutation$gene%in%diff_gene)
mm2<-mm[,1:7]
cancer_type<-rep("Liver Cancer",dim(mm2)[1])
mm3<-cbind(mm2$sample,cancer_type,mm2[,2:7])
colnames(mm3)<-c("Sample_ID","Cancer_Type","Chromosome","Start_Position","End_Position","Reference_Allele","Variant_Allele","Hugo_Symbol")
write.table(mm3,"f:/workplace/KIRP/KIRP_multidata_sex_diff/KIRP_MutationMaper.txt",row.names = F,col.names = T,quote = F,sep = "\t") #用来MutationMaper的输入突变文件


# 瀑布图 ---------------------------------------------------------------------

library(pheatmap)
dta_m2<-dta_m[,c("CACNA1D","PTEN")]
p<-t(dta_m2)
colnames(p)<-dta_m$sample
#rownames(p)<-"TP53"
#p<-p[-3,] #4个基因的突变01谱

max(p)

p[p == 1] <- "mut"
p[p == 0] <- ""
p[p == 2] <- "mut"
p[p == 3] <- "mut"

library(ComplexHeatmap)

mytype<-data.frame(sample=dta_m$sample,sex=ifelse(dta_m$sex==0,"Female","Male")) #样本的性别信息

#用fill = 设置mutation以及背景用什么颜色
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = NA, col = NA)) #不要背景色
  },
  mut = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#008000", col = NA)) #mut是绿色
  })

#mutation bar plot的颜色，跟瀑布图一致
col = c("mut" = "#008000")

#定义足够多的颜色
mycol <- c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

#提取瀑布图中画出的sample
p<-as.data.frame(p)
p1 <- oncoPrint(p[1:(ncol(p))], get_type = function(x) x,
                alter_fun = alter_fun, col = col,
                remove_empty_columns = TRUE #删除没有突变的sample
                #column_order = NULL, #不按突变频率给sample排序
                #row_order = NULL, #不按突变频率给基因排序
)


rownames(mytype) <- mytype$sample

#现在，这个mytype就跟瀑布图中的sample一致了。

#画临床数据的heatmap
my_annotation = HeatmapAnnotation(df = data.frame(mytype[2]),
                                  col = list(sex = c("Female" = mycol[2], "Male" = mycol[1])))

#画瀑布图
p2 <- oncoPrint(p[1:(ncol(p))], get_type = function(x) x,
                alter_fun = alter_fun, col = col,
                remove_empty_columns = TRUE,#删除没有突变的sample
                #column_order = NULL, #不按突变频率给sample排序
                #row_order = NULL, #不按突变频率给基因排序
                #show_pct = FALSE, #左侧不显示百分比
                bottom_annotation = my_annotation,#把临床信息画在下面
                show_heatmap_legend = FALSE) #不显示突变的图例
p2

