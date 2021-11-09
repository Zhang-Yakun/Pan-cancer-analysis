library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)
library("RColorBrewer")
display.brewer.all(type = "qual")

setwd("f:/workplace/结果2/")

# 47个通路免疫衰老基因 -------------------------------------------------------------
gene47<-read.table("f:/workplace/结果2/pathway_gene_47.txt")
path_immu_gene<-gene47$V1%>%as.character()

# 47免疫衰老基因的mRNA表达谱及其突变分组 ----------------------------------------------------------
############mRNA表达谱
BLCA_tpm<-read.table("F:/workplace/BLCA/BLCA_readcount.genes.tpm.txt")

BLCA_tpm<-(log2(BLCA_tpm+1))  #表达谱变为log2(TPM+1)
colnames(BLCA_tpm)<-str_replace_all(colnames(BLCA_tpm),"[.*]",replacement = "-")

BLCA_gene<-intersect(rownames(BLCA_tpm),path_immu_gene) 
BLCA_tpm2<-BLCA_tpm[BLCA_gene,]  #免疫衰老基因的表达谱

########突变分组
BLCA<-fread("f:/workplace/BLCA/BLCA_mc3.txt")%>%as.data.frame()%>%distinct()
BLCA_mut<-subset(BLCA,effect!="Silent") #去除沉默突变,剩16种突变类型

BLCA_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

BLCA_mut_47<-subset(BLCA_mut,BLCA_mut$gene %in% path_immu_gene)
#47基因在BLCA中的突变数据
BLCA_m<-table(BLCA_mut_47$sample,BLCA_mut_47$gene)
BLCA_mutat<-t(BLCA_m)%>%as.matrix.data.frame()
rownames(BLCA_mutat)<-colnames(BLCA_m)
colnames(BLCA_mutat)<-rownames(BLCA_m)

BLCA_l<-list()
for(i in 1:length(path_immu_gene)){
  
  x<-t(BLCA_tpm2[path_immu_gene[i],])
  BLCA_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  BLCA_exp_gene$sample<-as.character(BLCA_exp_gene$sample)
  
  x2<-t(BLCA_mutat[path_immu_gene[i],])
  BLCA_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  BLCA_mut_gene$sample<-as.character(BLCA_mut_gene$sample)
  
  BLCA_union<-inner_join(BLCA_exp_gene,BLCA_mut_gene)
  BLCA_union$cancer<-rep("BLCA",dim(BLCA_union)[1])
  BLCA_union$gene<-rep(path_immu_gene[i],dim(BLCA_union)[1])
  BLCA_union$mut_group<-ifelse(BLCA_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(BLCA_union)[2]<-"mRNA"
  BLCA_l<-rbind(BLCA_union,BLCA_l)
}
##############
LUAD_tpm<-read.table("F:/workplace/LUAD/LUAD_readcount.genes.tpm.txt")

LUAD_tpm<-(log2(LUAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LUAD_tpm)<-str_replace_all(colnames(LUAD_tpm),"[.*]",replacement = "-")

LUAD_gene<-intersect(rownames(LUAD_tpm),path_immu_gene) 
LUAD_tpm2<-LUAD_tpm[LUAD_gene,]  #免疫衰老基因的表达谱

########突变分组
LUAD<-fread("f:/workplace/LUAD/LUAD_mc3.txt")%>%as.data.frame()%>%distinct()
LUAD_mut<-subset(LUAD,effect!="Silent") #去除沉默突变,剩16种突变类型

LUAD_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

LUAD_mut_47<-subset(LUAD_mut,LUAD_mut$gene %in% path_immu_gene)
#47基因在LUAD中的突变数据
LUAD_m<-table(LUAD_mut_47$sample,LUAD_mut_47$gene)
LUAD_mutat<-t(LUAD_m)%>%as.matrix.data.frame()
rownames(LUAD_mutat)<-colnames(LUAD_m)
colnames(LUAD_mutat)<-rownames(LUAD_m)

LUAD_l<-list()
for(i in 1:length(path_immu_gene)){
  
  x<-t(LUAD_tpm2[path_immu_gene[i],])
  LUAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LUAD_exp_gene$sample<-as.character(LUAD_exp_gene$sample)
  
  x2<-t(LUAD_mutat[path_immu_gene[i],])
  LUAD_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  LUAD_mut_gene$sample<-as.character(LUAD_mut_gene$sample)
  
  LUAD_union<-inner_join(LUAD_exp_gene,LUAD_mut_gene)
  LUAD_union$cancer<-rep("LUAD",dim(LUAD_union)[1])
  LUAD_union$gene<-rep(path_immu_gene[i],dim(LUAD_union)[1])
  LUAD_union$mut_group<-ifelse(LUAD_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(LUAD_union)[2]<-"mRNA"
  LUAD_l<-rbind(LUAD_union,LUAD_l)
}
###########
HNSC_tpm<-read.table("F:/workplace/HNSC/HNSC_readcount.genes.tpm.txt")

HNSC_tpm<-(log2(HNSC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(HNSC_tpm)<-str_replace_all(colnames(HNSC_tpm),"[.*]",replacement = "-")

HNSC_gene<-intersect(rownames(HNSC_tpm),path_immu_gene) 
HNSC_tpm2<-HNSC_tpm[HNSC_gene,]  #免疫衰老基因的表达谱

########突变分组
HNSC<-fread("f:/workplace/HNSC/HNSC_mc3.txt")%>%as.data.frame()%>%distinct()
HNSC_mut<-subset(HNSC,effect!="Silent") #去除沉默突变,剩16种突变类型

HNSC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

HNSC_mut_47<-subset(HNSC_mut,HNSC_mut$gene %in% path_immu_gene)
#47基因在HNSC中的突变数据
HNSC_m<-table(HNSC_mut_47$sample,HNSC_mut_47$gene)
HNSC_mutat<-t(HNSC_m)%>%as.matrix.data.frame()
rownames(HNSC_mutat)<-colnames(HNSC_m)
colnames(HNSC_mutat)<-rownames(HNSC_m)

HNSC_l<-list()
for(i in 1:dim(HNSC_m)[2]){
  
  x<-t(HNSC_tpm2[colnames(HNSC_m)[i],])
  HNSC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  HNSC_exp_gene$sample<-as.character(HNSC_exp_gene$sample)
  
  x2<-t(HNSC_mutat[colnames(HNSC_m)[i],])
  HNSC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  HNSC_mut_gene$sample<-as.character(HNSC_mut_gene$sample)
  
  HNSC_union<-inner_join(HNSC_exp_gene,HNSC_mut_gene)
  HNSC_union$cancer<-rep("HNSC",dim(HNSC_union)[1])
  HNSC_union$gene<-rep(colnames(HNSC_m)[i],dim(HNSC_union)[1])
  HNSC_union$mut_group<-ifelse(HNSC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(HNSC_union)[2]<-"mRNA"
  HNSC_l<-rbind(HNSC_union,HNSC_l)
}

##########
KIRC_tpm<-read.table("F:/workplace/KIRC/KIRC_readcount.genes.tpm.txt")

KIRC_tpm<-(log2(KIRC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(KIRC_tpm)<-str_replace_all(colnames(KIRC_tpm),"[.*]",replacement = "-")

KIRC_gene<-intersect(rownames(KIRC_tpm),path_immu_gene) 
KIRC_tpm2<-KIRC_tpm[KIRC_gene,]  #免疫衰老基因的表达谱

########突变分组
KIRC<-fread("f:/workplace/KIRC/KIRC_mc3.txt")%>%as.data.frame()%>%distinct()
KIRC_mut<-subset(KIRC,effect!="Silent") #去除沉默突变,剩16种突变类型

KIRC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

KIRC_mut_47<-subset(KIRC_mut,KIRC_mut$gene %in% path_immu_gene)
#47基因在KIRC中的突变数据
KIRC_m<-table(KIRC_mut_47$sample,KIRC_mut_47$gene)
KIRC_mutat<-t(KIRC_m)%>%as.matrix.data.frame()
rownames(KIRC_mutat)<-colnames(KIRC_m)
colnames(KIRC_mutat)<-rownames(KIRC_m)

KIRC_l<-list()
for(i in 1:dim(KIRC_m)[2]){
  
  x<-t(KIRC_tpm2[colnames(KIRC_m)[i],])
  KIRC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  KIRC_exp_gene$sample<-as.character(KIRC_exp_gene$sample)
  
  x2<-t(KIRC_mutat[colnames(KIRC_m)[i],])
  KIRC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  KIRC_mut_gene$sample<-as.character(KIRC_mut_gene$sample)
  
  KIRC_union<-inner_join(KIRC_exp_gene,KIRC_mut_gene)
  KIRC_union$cancer<-rep("KIRC",dim(KIRC_union)[1])
  KIRC_union$gene<-rep(colnames(KIRC_m)[i],dim(KIRC_union)[1])
  KIRC_union$mut_group<-ifelse(KIRC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(KIRC_union)[2]<-"mRNA"
  KIRC_l<-rbind(KIRC_union,KIRC_l)
}

############
KIRP_tpm<-read.table("F:/workplace/KIRP/KIRP_readcount.genes.tpm.txt")

KIRP_tpm<-(log2(KIRP_tpm+1))  #表达谱变为log2(TPM+1)
colnames(KIRP_tpm)<-str_replace_all(colnames(KIRP_tpm),"[.*]",replacement = "-")

KIRP_gene<-intersect(rownames(KIRP_tpm),path_immu_gene) 
KIRP_tpm2<-KIRP_tpm[KIRP_gene,]  #免疫衰老基因的表达谱

########突变分组
KIRP<-fread("f:/workplace/KIRP/KIRP_mc3.txt")%>%as.data.frame()%>%distinct()
KIRP_mut<-subset(KIRP,effect!="Silent") #去除沉默突变,剩16种突变类型

KIRP_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

KIRP_mut_47<-subset(KIRP_mut,KIRP_mut$gene %in% path_immu_gene)
#47基因在KIRP中的突变数据
KIRP_m<-table(KIRP_mut_47$sample,KIRP_mut_47$gene)
KIRP_mutat<-t(KIRP_m)%>%as.matrix.data.frame()
rownames(KIRP_mutat)<-colnames(KIRP_m)
colnames(KIRP_mutat)<-rownames(KIRP_m)

KIRP_l<-list()
for(i in 1:dim(KIRP_m)[2]){
  
  x<-t(KIRP_tpm2[colnames(KIRP_m)[i],])
  KIRP_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  KIRP_exp_gene$sample<-as.character(KIRP_exp_gene$sample)
  
  x2<-t(KIRP_mutat[colnames(KIRP_m)[i],])
  KIRP_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  KIRP_mut_gene$sample<-as.character(KIRP_mut_gene$sample)
  
  KIRP_union<-inner_join(KIRP_exp_gene,KIRP_mut_gene)
  KIRP_union$cancer<-rep("KIRP",dim(KIRP_union)[1])
  KIRP_union$gene<-rep(colnames(KIRP_m)[i],dim(KIRP_union)[1])
  KIRP_union$mut_group<-ifelse(KIRP_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(KIRP_union)[2]<-"mRNA"
  KIRP_l<-rbind(KIRP_union,KIRP_l)
}

#############

LIHC_tpm<-read.table("F:/workplace/LIHC/LIHC_readcount.genes.tpm.txt")

LIHC_tpm<-(log2(LIHC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LIHC_tpm)<-str_replace_all(colnames(LIHC_tpm),"[.*]",replacement = "-")

LIHC_gene<-intersect(rownames(LIHC_tpm),path_immu_gene) 
LIHC_tpm2<-LIHC_tpm[LIHC_gene,]  #免疫衰老基因的表达谱

########突变分组
LIHC<-fread("f:/workplace/LIHC/LIHC_mc3.txt")%>%as.data.frame()%>%distinct()
LIHC_mut<-subset(LIHC,effect!="Silent") #去除沉默突变,剩16种突变类型

LIHC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

LIHC_mut_47<-subset(LIHC_mut,LIHC_mut$gene %in% path_immu_gene)
#47基因在LIHC中的突变数据
LIHC_m<-table(LIHC_mut_47$sample,LIHC_mut_47$gene)
LIHC_mutat<-t(LIHC_m)%>%as.matrix.data.frame()
rownames(LIHC_mutat)<-colnames(LIHC_m)
colnames(LIHC_mutat)<-rownames(LIHC_m)

LIHC_l<-list()
for(i in 1:dim(LIHC_m)[2]){
  
  x<-t(LIHC_tpm2[colnames(LIHC_m)[i],])
  LIHC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LIHC_exp_gene$sample<-as.character(LIHC_exp_gene$sample)
  
  x2<-t(LIHC_mutat[colnames(LIHC_m)[i],])
  LIHC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  LIHC_mut_gene$sample<-as.character(LIHC_mut_gene$sample)
  
  LIHC_union<-inner_join(LIHC_exp_gene,LIHC_mut_gene)
  LIHC_union$cancer<-rep("LIHC",dim(LIHC_union)[1])
  LIHC_union$gene<-rep(colnames(LIHC_m)[i],dim(LIHC_union)[1])
  LIHC_union$mut_group<-ifelse(LIHC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(LIHC_union)[2]<-"mRNA"
  LIHC_l<-rbind(LIHC_union,LIHC_l)
}

#############
LUSC_tpm<-read.table("F:/workplace/LUSC/LUSC_readcount.genes.tpm.txt")

LUSC_tpm<-(log2(LUSC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LUSC_tpm)<-str_replace_all(colnames(LUSC_tpm),"[.*]",replacement = "-")

LUSC_gene<-intersect(rownames(LUSC_tpm),path_immu_gene) 
LUSC_tpm2<-LUSC_tpm[LUSC_gene,]  #免疫衰老基因的表达谱

########突变分组
LUSC<-fread("f:/workplace/LUSC/LUSC_mc3.txt")%>%as.data.frame()%>%distinct()
LUSC_mut<-subset(LUSC,effect!="Silent") #去除沉默突变,剩16种突变类型

LUSC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

LUSC_mut_47<-subset(LUSC_mut,LUSC_mut$gene %in% path_immu_gene)
#47基因在LUSC中的突变数据
LUSC_m<-table(LUSC_mut_47$sample,LUSC_mut_47$gene)
LUSC_mutat<-t(LUSC_m)%>%as.matrix.data.frame()
rownames(LUSC_mutat)<-colnames(LUSC_m)
colnames(LUSC_mutat)<-rownames(LUSC_m)

LUSC_l<-list()
for(i in 1:dim(LUSC_m)[2]){
  
  x<-t(LUSC_tpm2[colnames(LUSC_m)[i],])
  LUSC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LUSC_exp_gene$sample<-as.character(LUSC_exp_gene$sample)
  
  x2<-t(LUSC_mutat[colnames(LUSC_m)[i],])
  LUSC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  LUSC_mut_gene$sample<-as.character(LUSC_mut_gene$sample)
  
  LUSC_union<-inner_join(LUSC_exp_gene,LUSC_mut_gene)
  LUSC_union$cancer<-rep("LUSC",dim(LUSC_union)[1])
  LUSC_union$gene<-rep(colnames(LUSC_m)[i],dim(LUSC_union)[1])
  LUSC_union$mut_group<-ifelse(LUSC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(LUSC_union)[2]<-"mRNA"
  LUSC_l<-rbind(LUSC_union,LUSC_l)
}
##############
THCA_tpm<-read.table("F:/workplace/THCA/THCA_readcount.genes.tpm.txt")

THCA_tpm<-(log2(THCA_tpm+1))  #表达谱变为log2(TPM+1)
colnames(THCA_tpm)<-str_replace_all(colnames(THCA_tpm),"[.*]",replacement = "-")

THCA_gene<-intersect(rownames(THCA_tpm),path_immu_gene) 
THCA_tpm2<-THCA_tpm[THCA_gene,]  #免疫衰老基因的表达谱

########突变分组
THCA<-fread("f:/workplace/THCA/THCA_mc3.txt")%>%as.data.frame()%>%distinct()
THCA_mut<-subset(THCA,effect!="Silent") #去除沉默突变,剩16种突变类型

THCA_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

THCA_mut_47<-subset(THCA_mut,THCA_mut$gene %in% path_immu_gene)
#47基因在THCA中的突变数据
THCA_m<-table(THCA_mut_47$sample,THCA_mut_47$gene)
THCA_mutat<-t(THCA_m)%>%as.matrix.data.frame()
rownames(THCA_mutat)<-colnames(THCA_m)
colnames(THCA_mutat)<-rownames(THCA_m)

THCA_l<-list()
for(i in 1:dim(THCA_m)[2]){
  
  x<-t(THCA_tpm2[colnames(THCA_m)[i],])
  THCA_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  THCA_exp_gene$sample<-as.character(THCA_exp_gene$sample)
  
  x2<-t(THCA_mutat[colnames(THCA_m)[i],])
  THCA_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  THCA_mut_gene$sample<-as.character(THCA_mut_gene$sample)
  
  THCA_union<-inner_join(THCA_exp_gene,THCA_mut_gene)
  THCA_union$cancer<-rep("THCA",dim(THCA_union)[1])
  THCA_union$gene<-rep(colnames(THCA_m)[i],dim(THCA_union)[1])
  THCA_union$mut_group<-ifelse(THCA_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(THCA_union)[2]<-"mRNA"
  THCA_l<-rbind(THCA_union,THCA_l)
}
##########
COAD_tpm<-read.table("F:/workplace/COAD/COAD_readcount.genes.tpm.txt")

COAD_tpm<-(log2(COAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(COAD_tpm)<-str_replace_all(colnames(COAD_tpm),"[.*]",replacement = "-")

COAD_gene<-intersect(rownames(COAD_tpm),path_immu_gene) 
COAD_tpm2<-COAD_tpm[COAD_gene,]  #免疫衰老基因的表达谱

########突变分组
COAD<-fread("f:/workplace/COAD/COAD_mc3.txt")%>%as.data.frame()%>%distinct()
COAD_mut<-subset(COAD,effect!="Silent") #去除沉默突变,剩16种突变类型

COAD_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

COAD_mut_47<-subset(COAD_mut,COAD_mut$gene %in% path_immu_gene)
#47基因在COAD中的突变数据
COAD_m<-table(COAD_mut_47$sample,COAD_mut_47$gene)
COAD_mutat<-t(COAD_m)%>%as.matrix.data.frame()
rownames(COAD_mutat)<-colnames(COAD_m)
colnames(COAD_mutat)<-rownames(COAD_m)

COAD_l<-list()
for(i in 1:dim(COAD_m)[2]){
  
  x<-t(COAD_tpm2[colnames(COAD_m)[i],])
  COAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  COAD_exp_gene$sample<-as.character(COAD_exp_gene$sample)
  
  x2<-t(COAD_mutat[colnames(COAD_m)[i],])
  COAD_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  COAD_mut_gene$sample<-as.character(COAD_mut_gene$sample)
  
  COAD_union<-inner_join(COAD_exp_gene,COAD_mut_gene)
  COAD_union$cancer<-rep("COAD",dim(COAD_union)[1])
  COAD_union$gene<-rep(colnames(COAD_m)[i],dim(COAD_union)[1])
  COAD_union$mut_group<-ifelse(COAD_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(COAD_union)[2]<-"mRNA"
  COAD_l<-rbind(COAD_union,COAD_l)
}
##############
GBM_tpm<-read.table("F:/workplace/GBM/GBM_readcount.genes.tpm.txt")

GBM_tpm<-(log2(GBM_tpm+1))  #表达谱变为log2(TPM+1)
colnames(GBM_tpm)<-str_replace_all(colnames(GBM_tpm),"[.*]",replacement = "-")

GBM_gene<-intersect(rownames(GBM_tpm),path_immu_gene) 
GBM_tpm2<-GBM_tpm[GBM_gene,]  #免疫衰老基因的表达谱

########突变分组
GBM<-fread("f:/workplace/GBM/GBM_mc3.txt")%>%as.data.frame()%>%distinct()
GBM_mut<-subset(GBM,effect!="Silent") #去除沉默突变,剩16种突变类型

GBM_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

GBM_mut_47<-subset(GBM_mut,GBM_mut$gene %in% path_immu_gene)
#47基因在GBM中的突变数据
GBM_m<-table(GBM_mut_47$sample,GBM_mut_47$gene)
GBM_mutat<-t(GBM_m)%>%as.matrix.data.frame()
rownames(GBM_mutat)<-colnames(GBM_m)
colnames(GBM_mutat)<-rownames(GBM_m)

GBM_l<-list()
for(i in 1:dim(GBM_m)[2]){
  
  x<-t(GBM_tpm2[colnames(GBM_m)[i],])
  GBM_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  GBM_exp_gene$sample<-as.character(GBM_exp_gene$sample)
  
  x2<-t(GBM_mutat[colnames(GBM_m)[i],])
  GBM_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  GBM_mut_gene$sample<-as.character(GBM_mut_gene$sample)
  
  GBM_union<-inner_join(GBM_exp_gene,GBM_mut_gene)
  GBM_union$cancer<-rep("GBM",dim(GBM_union)[1])
  GBM_union$gene<-rep(colnames(GBM_m)[i],dim(GBM_union)[1])
  GBM_union$mut_group<-ifelse(GBM_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(GBM_union)[2]<-"mRNA"
  GBM_l<-rbind(GBM_union,GBM_l)
}
##########3
KICH_tpm<-read.table("F:/workplace/KICH/KICH_readcount.genes.tpm.txt")

KICH_tpm<-(log2(KICH_tpm+1))  #表达谱变为log2(TPM+1)
colnames(KICH_tpm)<-str_replace_all(colnames(KICH_tpm),"[.*]",replacement = "-")

KICH_gene<-intersect(rownames(KICH_tpm),path_immu_gene) 
KICH_tpm2<-KICH_tpm[KICH_gene,]  #免疫衰老基因的表达谱

########突变分组
KICH<-fread("f:/workplace/KICH/KICH_mc3.txt")%>%as.data.frame()%>%distinct()
KICH_mut<-subset(KICH,effect!="Silent") #去除沉默突变,剩16种突变类型

KICH_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

KICH_mut_47<-subset(KICH_mut,KICH_mut$gene %in% path_immu_gene)
#47基因在KICH中的突变数据
KICH_m<-table(KICH_mut_47$sample,KICH_mut_47$gene)
KICH_mutat<-t(KICH_m)%>%as.matrix.data.frame()
rownames(KICH_mutat)<-colnames(KICH_m)
colnames(KICH_mutat)<-rownames(KICH_m)

KICH_l<-list()
for(i in 1:dim(KICH_m)[2]){
  
  x<-t(KICH_tpm2[colnames(KICH_m)[i],])
  KICH_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  KICH_exp_gene$sample<-as.character(KICH_exp_gene$sample)
  
  x2<-t(KICH_mutat[colnames(KICH_m)[i],])
  KICH_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  KICH_mut_gene$sample<-as.character(KICH_mut_gene$sample)
  
  KICH_union<-inner_join(KICH_exp_gene,KICH_mut_gene)
  KICH_union$cancer<-rep("KICH",dim(KICH_union)[1])
  KICH_union$gene<-rep(colnames(KICH_m)[i],dim(KICH_union)[1])
  KICH_union$mut_group<-ifelse(KICH_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(KICH_union)[2]<-"mRNA"
  KICH_l<-rbind(KICH_union,KICH_l)
}
############
LGG_tpm<-read.table("F:/workplace/LGG/LGG_readcount.genes.tpm.txt")

LGG_tpm<-(log2(LGG_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LGG_tpm)<-str_replace_all(colnames(LGG_tpm),"[.*]",replacement = "-")

LGG_gene<-intersect(rownames(LGG_tpm),path_immu_gene) 
LGG_tpm2<-LGG_tpm[LGG_gene,]  #免疫衰老基因的表达谱

########突变分组
LGG<-fread("f:/workplace/LGG/LGG_mc3.txt")%>%as.data.frame()%>%distinct()
LGG_mut<-subset(LGG,effect!="Silent") #去除沉默突变,剩16种突变类型

LGG_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

LGG_mut_47<-subset(LGG_mut,LGG_mut$gene %in% path_immu_gene)
#47基因在LGG中的突变数据
LGG_m<-table(LGG_mut_47$sample,LGG_mut_47$gene)
LGG_mutat<-t(LGG_m)%>%as.matrix.data.frame()
rownames(LGG_mutat)<-colnames(LGG_m)
colnames(LGG_mutat)<-rownames(LGG_m)

LGG_l<-list()
for(i in 1:dim(LGG_m)[2]){
  
  x<-t(LGG_tpm2[colnames(LGG_m)[i],])
  LGG_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LGG_exp_gene$sample<-as.character(LGG_exp_gene$sample)
  
  x2<-t(LGG_mutat[colnames(LGG_m)[i],])
  LGG_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  LGG_mut_gene$sample<-as.character(LGG_mut_gene$sample)
  
  LGG_union<-inner_join(LGG_exp_gene,LGG_mut_gene)
  LGG_union$cancer<-rep("LGG",dim(LGG_union)[1])
  LGG_union$gene<-rep(colnames(LGG_m)[i],dim(LGG_union)[1])
  LGG_union$mut_group<-ifelse(LGG_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(LGG_union)[2]<-"mRNA"
  LGG_l<-rbind(LGG_union,LGG_l)
}
#############
BRCA_tpm<-read.table("F:/workplace/BRCA/BRCA_readcount.genes.tpm.txt")

BRCA_tpm<-(log2(BRCA_tpm+1))  #表达谱变为log2(TPM+1)
colnames(BRCA_tpm)<-str_replace_all(colnames(BRCA_tpm),"[.*]",replacement = "-")

BRCA_gene<-intersect(rownames(BRCA_tpm),path_immu_gene) 
BRCA_tpm2<-BRCA_tpm[BRCA_gene,]  #免疫衰老基因的表达谱

########突变分组
BRCA<-fread("f:/workplace/BRCA/BRCA_mc3.txt")%>%as.data.frame()%>%distinct()
BRCA_mut<-subset(BRCA,effect!="Silent") #去除沉默突变,剩16种突变类型

BRCA_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

BRCA_mut_47<-subset(BRCA_mut,BRCA_mut$gene %in% path_immu_gene)
#47基因在BRCA中的突变数据
BRCA_m<-table(BRCA_mut_47$sample,BRCA_mut_47$gene)
BRCA_mutat<-t(BRCA_m)%>%as.matrix.data.frame()
rownames(BRCA_mutat)<-colnames(BRCA_m)
colnames(BRCA_mutat)<-rownames(BRCA_m)

BRCA_l<-list()
for(i in 1:dim(BRCA_m)[2]){
  
  x<-t(BRCA_tpm2[colnames(BRCA_m)[i],])
  BRCA_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  BRCA_exp_gene$sample<-as.character(BRCA_exp_gene$sample)
  
  x2<-t(BRCA_mutat[colnames(BRCA_m)[i],])
  BRCA_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  BRCA_mut_gene$sample<-as.character(BRCA_mut_gene$sample)
  
  BRCA_union<-inner_join(BRCA_exp_gene,BRCA_mut_gene)
  BRCA_union$cancer<-rep("BRCA",dim(BRCA_union)[1])
  BRCA_union$gene<-rep(colnames(BRCA_m)[i],dim(BRCA_union)[1])
  BRCA_union$mut_group<-ifelse(BRCA_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(BRCA_union)[2]<-"mRNA"
  BRCA_l<-rbind(BRCA_union,BRCA_l)
}
############
UCS_tpm<-read.table("F:/workplace/UCS/UCS_readcount.genes.tpm.txt")

UCS_tpm<-(log2(UCS_tpm+1))  #表达谱变为log2(TPM+1)
colnames(UCS_tpm)<-str_replace_all(colnames(UCS_tpm),"[.*]",replacement = "-")

UCS_gene<-intersect(rownames(UCS_tpm),path_immu_gene) 
UCS_tpm2<-UCS_tpm[UCS_gene,]  #免疫衰老基因的表达谱

########突变分组
UCS<-fread("f:/workplace/UCS/UCS_mc3.txt")%>%as.data.frame()%>%distinct()
UCS_mut<-subset(UCS,effect!="Silent") #去除沉默突变,剩16种突变类型

UCS_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

UCS_mut_47<-subset(UCS_mut,UCS_mut$gene %in% path_immu_gene)
#47基因在UCS中的突变数据
UCS_m<-table(UCS_mut_47$sample,UCS_mut_47$gene)
UCS_mutat<-t(UCS_m)%>%as.matrix.data.frame()
rownames(UCS_mutat)<-colnames(UCS_m)
colnames(UCS_mutat)<-rownames(UCS_m)

UCS_l<-list()
for(i in 1:dim(UCS_m)[2]){
  
  x<-t(UCS_tpm2[colnames(UCS_m)[i],])
  UCS_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  UCS_exp_gene$sample<-as.character(UCS_exp_gene$sample)
  
  x2<-t(UCS_mutat[colnames(UCS_m)[i],])
  UCS_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  UCS_mut_gene$sample<-as.character(UCS_mut_gene$sample)
  
  UCS_union<-inner_join(UCS_exp_gene,UCS_mut_gene)
  UCS_union$cancer<-rep("UCS",dim(UCS_union)[1])
  UCS_union$gene<-rep(colnames(UCS_m)[i],dim(UCS_union)[1])
  UCS_union$mut_group<-ifelse(UCS_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(UCS_union)[2]<-"mRNA"
  UCS_l<-rbind(UCS_union,UCS_l)
}
###########
UCEC_tpm<-read.table("F:/workplace/UCEC/UCEC_readcount.genes.tpm.txt")

UCEC_tpm<-(log2(UCEC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(UCEC_tpm)<-str_replace_all(colnames(UCEC_tpm),"[.*]",replacement = "-")

UCEC_gene<-intersect(rownames(UCEC_tpm),path_immu_gene) 
UCEC_tpm2<-UCEC_tpm[UCEC_gene,]  #免疫衰老基因的表达谱

########突变分组
UCEC<-fread("f:/workplace/UCEC/UCEC_mc3.txt")%>%as.data.frame()%>%distinct()
UCEC_mut<-subset(UCEC,effect!="Silent") #去除沉默突变,剩16种突变类型

UCEC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

UCEC_mut_47<-subset(UCEC_mut,UCEC_mut$gene %in% path_immu_gene)
#47基因在UCEC中的突变数据
UCEC_m<-table(UCEC_mut_47$sample,UCEC_mut_47$gene)
UCEC_mutat<-t(UCEC_m)%>%as.matrix.data.frame()
rownames(UCEC_mutat)<-colnames(UCEC_m)
colnames(UCEC_mutat)<-rownames(UCEC_m)

UCEC_l<-list()
for(i in 1:dim(UCEC_m)[2]){
  
  x<-t(UCEC_tpm2[colnames(UCEC_m)[i],])
  UCEC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  UCEC_exp_gene$sample<-as.character(UCEC_exp_gene$sample)
  
  x2<-t(UCEC_mutat[colnames(UCEC_m)[i],])
  UCEC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  UCEC_mut_gene$sample<-as.character(UCEC_mut_gene$sample)
  
  UCEC_union<-inner_join(UCEC_exp_gene,UCEC_mut_gene)
  UCEC_union$cancer<-rep("UCEC",dim(UCEC_union)[1])
  UCEC_union$gene<-rep(colnames(UCEC_m)[i],dim(UCEC_union)[1])
  UCEC_union$mut_group<-ifelse(UCEC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(UCEC_union)[2]<-"mRNA"
  UCEC_l<-rbind(UCEC_union,UCEC_l)
}
#############
OV_tpm<-read.table("F:/workplace/OV/OV_readcount.genes.tpm.txt")

OV_tpm<-(log2(OV_tpm+1))  #表达谱变为log2(TPM+1)
colnames(OV_tpm)<-str_replace_all(colnames(OV_tpm),"[.*]",replacement = "-")

OV_gene<-intersect(rownames(OV_tpm),path_immu_gene) 
OV_tpm2<-OV_tpm[OV_gene,]  #免疫衰老基因的表达谱
########突变分组
OV<-fread("f:/workplace/OV/OV_mc3.txt")%>%as.data.frame()%>%distinct()
OV_mut<-subset(OV,effect!="Silent") #去除沉默突变,剩16种突变类型

OV_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

OV_mut_47<-subset(OV_mut,OV_mut$gene %in% path_immu_gene)
#47基因在OV中的突变数据
OV_m<-table(OV_mut_47$sample,OV_mut_47$gene)
OV_mutat<-t(OV_m)%>%as.matrix.data.frame()
rownames(OV_mutat)<-colnames(OV_m)
colnames(OV_mutat)<-rownames(OV_m)

OV_l<-list()
for(i in 1:dim(OV_m)[2]){
  
  x<-t(OV_tpm2[colnames(OV_m)[i],])
  OV_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  OV_exp_gene$sample<-as.character(OV_exp_gene$sample)
  
  x2<-t(OV_mutat[colnames(OV_m)[i],])
  OV_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  OV_mut_gene$sample<-as.character(OV_mut_gene$sample)
  
  OV_union<-inner_join(OV_exp_gene,OV_mut_gene)
  OV_union$cancer<-rep("OV",dim(OV_union)[1])
  OV_union$gene<-rep(colnames(OV_m)[i],dim(OV_union)[1])
  OV_union$mut_group<-ifelse(OV_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(OV_union)[2]<-"mRNA"
  OV_l<-rbind(OV_union,OV_l)
}
############
PRAD_tpm<-read.table("F:/workplace/PRAD/PRAD_readcount.genes.tpm.txt")

PRAD_tpm<-(log2(PRAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(PRAD_tpm)<-str_replace_all(colnames(PRAD_tpm),"[.*]",replacement = "-")

PRAD_gene<-intersect(rownames(PRAD_tpm),path_immu_gene) 
PRAD_tpm2<-PRAD_tpm[PRAD_gene,]  #免疫衰老基因的表达谱

########突变分组
PRAD<-fread("f:/workplace/PRAD/PRAD_mc3.txt")%>%as.data.frame()%>%distinct()
PRAD_mut<-subset(PRAD,effect!="Silent") #去除沉默突变,剩16种突变类型

PRAD_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

PRAD_mut_47<-subset(PRAD_mut,PRAD_mut$gene %in% path_immu_gene)
#47基因在PRAD中的突变数据
PRAD_m<-table(PRAD_mut_47$sample,PRAD_mut_47$gene)
PRAD_mutat<-t(PRAD_m)%>%as.matrix.data.frame()
rownames(PRAD_mutat)<-colnames(PRAD_m)
colnames(PRAD_mutat)<-rownames(PRAD_m)

PRAD_l<-list()
for(i in 1:dim(PRAD_m)[2]){
  
  x<-t(PRAD_tpm2[colnames(PRAD_m)[i],])
  PRAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  PRAD_exp_gene$sample<-as.character(PRAD_exp_gene$sample)
  
  x2<-t(PRAD_mutat[colnames(PRAD_m)[i],])
  PRAD_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  PRAD_mut_gene$sample<-as.character(PRAD_mut_gene$sample)
  
  PRAD_union<-inner_join(PRAD_exp_gene,PRAD_mut_gene)
  PRAD_union$cancer<-rep("PRAD",dim(PRAD_union)[1])
  PRAD_union$gene<-rep(colnames(PRAD_m)[i],dim(PRAD_union)[1])
  PRAD_union$mut_group<-ifelse(PRAD_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(PRAD_union)[2]<-"mRNA"
  PRAD_l<-rbind(PRAD_union,PRAD_l)
}
#########
STAD_tpm<-read.table("F:/workplace/STAD/STAD_readcount.genes.tpm.txt")

STAD_tpm<-(log2(STAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(STAD_tpm)<-str_replace_all(colnames(STAD_tpm),"[.*]",replacement = "-")

STAD_gene<-intersect(rownames(STAD_tpm),path_immu_gene) 
STAD_tpm2<-STAD_tpm[STAD_gene,]  #免疫衰老基因的表达谱

########突变分组
STAD<-fread("f:/workplace/STAD/STAD_mc3.txt")%>%as.data.frame()%>%distinct()
STAD_mut<-subset(STAD,effect!="Silent") #去除沉默突变,剩16种突变类型

STAD_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

STAD_mut_47<-subset(STAD_mut,STAD_mut$gene %in% path_immu_gene)
#47基因在STAD中的突变数据
STAD_m<-table(STAD_mut_47$sample,STAD_mut_47$gene)
STAD_mutat<-t(STAD_m)%>%as.matrix.data.frame()
rownames(STAD_mutat)<-colnames(STAD_m)
colnames(STAD_mutat)<-rownames(STAD_m)

STAD_l<-list()
for(i in 1:dim(STAD_m)[2]){
  
  x<-t(STAD_tpm2[colnames(STAD_m)[i],])
  STAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  STAD_exp_gene$sample<-as.character(STAD_exp_gene$sample)
  
  x2<-t(STAD_mutat[colnames(STAD_m)[i],])
  STAD_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  STAD_mut_gene$sample<-as.character(STAD_mut_gene$sample)
  
  STAD_union<-inner_join(STAD_exp_gene,STAD_mut_gene)
  STAD_union$cancer<-rep("STAD",dim(STAD_union)[1])
  STAD_union$gene<-rep(colnames(STAD_m)[i],dim(STAD_union)[1])
  STAD_union$mut_group<-ifelse(STAD_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(STAD_union)[2]<-"mRNA"
  STAD_l<-rbind(STAD_union,STAD_l)
}
############
SKCM_tpm<-read.table("F:/workplace/SKCM/SKCM_readcount.genes.tpm.txt")

SKCM_tpm<-(log2(SKCM_tpm+1))  #表达谱变为log2(TPM+1)
colnames(SKCM_tpm)<-str_replace_all(colnames(SKCM_tpm),"[.*]",replacement = "-")

SKCM_gene<-intersect(rownames(SKCM_tpm),path_immu_gene) 
SKCM_tpm2<-SKCM_tpm[SKCM_gene,]  #免疫衰老基因的表达谱

########突变分组
SKCM<-fread("f:/workplace/SKCM/SKCM_mc3.txt")%>%as.data.frame()%>%distinct()
SKCM_mut<-subset(SKCM,effect!="Silent") #去除沉默突变,剩16种突变类型

SKCM_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

SKCM_mut_47<-subset(SKCM_mut,SKCM_mut$gene %in% path_immu_gene)
#47基因在SKCM中的突变数据
SKCM_m<-table(SKCM_mut_47$sample,SKCM_mut_47$gene)
SKCM_mutat<-t(SKCM_m)%>%as.matrix.data.frame()
rownames(SKCM_mutat)<-colnames(SKCM_m)
colnames(SKCM_mutat)<-rownames(SKCM_m)

SKCM_l<-list()
for(i in 1:dim(SKCM_m)[2]){
  
  x<-t(SKCM_tpm2[colnames(SKCM_m)[i],])
  SKCM_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  SKCM_exp_gene$sample<-as.character(SKCM_exp_gene$sample)
  
  x2<-t(SKCM_mutat[colnames(SKCM_m)[i],])
  SKCM_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  SKCM_mut_gene$sample<-as.character(SKCM_mut_gene$sample)
  
  SKCM_union<-inner_join(SKCM_exp_gene,SKCM_mut_gene)
  SKCM_union$cancer<-rep("SKCM",dim(SKCM_union)[1])
  SKCM_union$gene<-rep(colnames(SKCM_m)[i],dim(SKCM_union)[1])
  SKCM_union$mut_group<-ifelse(SKCM_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(SKCM_union)[2]<-"mRNA"
  SKCM_l<-rbind(SKCM_union,SKCM_l)
}
#########
CESC_tpm<-read.table("F:/workplace/CESC/CESC_readcount.genes.tpm.txt")

CESC_tpm<-(log2(CESC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(CESC_tpm)<-str_replace_all(colnames(CESC_tpm),"[.*]",replacement = "-")

CESC_gene<-intersect(rownames(CESC_tpm),path_immu_gene) 
CESC_tpm2<-CESC_tpm[CESC_gene,]  #免疫衰老基因的表达谱

########突变分组
CESC<-fread("f:/workplace/CESC/CESC_mc3.txt")%>%as.data.frame()%>%distinct()
CESC_mut<-subset(CESC,effect!="Silent") #去除沉默突变,剩16种突变类型

CESC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

CESC_mut_47<-subset(CESC_mut,CESC_mut$gene %in% path_immu_gene)
#47基因在CESC中的突变数据
CESC_m<-table(CESC_mut_47$sample,CESC_mut_47$gene)
CESC_mutat<-t(CESC_m)%>%as.matrix.data.frame()
rownames(CESC_mutat)<-colnames(CESC_m)
colnames(CESC_mutat)<-rownames(CESC_m)

CESC_l<-list()
for(i in 1:dim(CESC_m)[2]){
  
  x<-t(CESC_tpm2[colnames(CESC_m)[i],])
  CESC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  CESC_exp_gene$sample<-as.character(CESC_exp_gene$sample)
  
  x2<-t(CESC_mutat[colnames(CESC_m)[i],])
  CESC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  CESC_mut_gene$sample<-as.character(CESC_mut_gene$sample)
  
  CESC_union<-inner_join(CESC_exp_gene,CESC_mut_gene)
  CESC_union$cancer<-rep("CESC",dim(CESC_union)[1])
  CESC_union$gene<-rep(colnames(CESC_m)[i],dim(CESC_union)[1])
  CESC_union$mut_group<-ifelse(CESC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(CESC_union)[2]<-"mRNA"
  CESC_l<-rbind(CESC_union,CESC_l)
}
##########
ACC_tpm<-read.table("F:/workplace/ACC/ACC_readcount.genes.tpm.txt")

ACC_tpm<-(log2(ACC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(ACC_tpm)<-str_replace_all(colnames(ACC_tpm),"[.*]",replacement = "-")

ACC_gene<-intersect(rownames(ACC_tpm),path_immu_gene) 
ACC_tpm2<-ACC_tpm[ACC_gene,]  #免疫衰老基因的表达谱

########突变分组
ACC<-fread("f:/workplace/ACC/ACC_mc3.txt")%>%as.data.frame()%>%distinct()
ACC_mut<-subset(ACC,effect!="Silent") #去除沉默突变,剩16种突变类型

ACC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

ACC_mut_47<-subset(ACC_mut,ACC_mut$gene %in% path_immu_gene)
#47基因在ACC中的突变数据
ACC_m<-table(ACC_mut_47$sample,ACC_mut_47$gene)
ACC_mutat<-t(ACC_m)%>%as.matrix.data.frame()
rownames(ACC_mutat)<-colnames(ACC_m)
colnames(ACC_mutat)<-rownames(ACC_m)

ACC_l<-list()
for(i in 1:dim(ACC_m)[2]){
  
  x<-t(ACC_tpm2[colnames(ACC_m)[i],])
  ACC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  ACC_exp_gene$sample<-as.character(ACC_exp_gene$sample)
  
  x2<-t(ACC_mutat[colnames(ACC_m)[i],])
  ACC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  ACC_mut_gene$sample<-as.character(ACC_mut_gene$sample)
  
  ACC_union<-inner_join(ACC_exp_gene,ACC_mut_gene)
  ACC_union$cancer<-rep("ACC",dim(ACC_union)[1])
  ACC_union$gene<-rep(colnames(ACC_m)[i],dim(ACC_union)[1])
  ACC_union$mut_group<-ifelse(ACC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(ACC_union)[2]<-"mRNA"
  ACC_l<-rbind(ACC_union,ACC_l)
}
###########
PCPG_tpm<-read.table("F:/workplace/PCPG/PCPG_readcount.genes.tpm.txt")

PCPG_tpm<-(log2(PCPG_tpm+1))  #表达谱变为log2(TPM+1)
colnames(PCPG_tpm)<-str_replace_all(colnames(PCPG_tpm),"[.*]",replacement = "-")

PCPG_gene<-intersect(rownames(PCPG_tpm),path_immu_gene) 
PCPG_tpm2<-PCPG_tpm[PCPG_gene,]  #免疫衰老基因的表达谱

########突变分组
PCPG<-fread("f:/workplace/PCPG/PCPG_mc3.txt")%>%as.data.frame()%>%distinct()
PCPG_mut<-subset(PCPG,effect!="Silent") #去除沉默突变,剩16种突变类型

PCPG_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

PCPG_mut_47<-subset(PCPG_mut,PCPG_mut$gene %in% path_immu_gene)
#47基因在PCPG中的突变数据
PCPG_m<-table(PCPG_mut_47$sample,PCPG_mut_47$gene)
PCPG_mutat<-t(PCPG_m)%>%as.matrix.data.frame()
rownames(PCPG_mutat)<-colnames(PCPG_m)
colnames(PCPG_mutat)<-rownames(PCPG_m)

PCPG_l<-list()
for(i in 1:dim(PCPG_m)[2]){
  
  x<-t(PCPG_tpm2[colnames(PCPG_m)[i],])
  PCPG_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  PCPG_exp_gene$sample<-as.character(PCPG_exp_gene$sample)
  
  x2<-t(PCPG_mutat[colnames(PCPG_m)[i],])
  PCPG_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  PCPG_mut_gene$sample<-as.character(PCPG_mut_gene$sample)
  
  PCPG_union<-inner_join(PCPG_exp_gene,PCPG_mut_gene)
  PCPG_union$cancer<-rep("PCPG",dim(PCPG_union)[1])
  PCPG_union$gene<-rep(colnames(PCPG_m)[i],dim(PCPG_union)[1])
  PCPG_union$mut_group<-ifelse(PCPG_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(PCPG_union)[2]<-"mRNA"
  PCPG_l<-rbind(PCPG_union,PCPG_l)
}
###########
SARC_tpm<-read.table("F:/workplace/SARC/SARC_readcount.genes.tpm.txt")

SARC_tpm<-(log2(SARC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(SARC_tpm)<-str_replace_all(colnames(SARC_tpm),"[.*]",replacement = "-")

SARC_gene<-intersect(rownames(SARC_tpm),path_immu_gene) 
SARC_tpm2<-SARC_tpm[SARC_gene,]  #免疫衰老基因的表达谱
########突变分组
SARC<-fread("f:/workplace/SARC/SARC_mc3.txt")%>%as.data.frame()%>%distinct()
SARC_mut<-subset(SARC,effect!="Silent") #去除沉默突变,剩16种突变类型

SARC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

SARC_mut_47<-subset(SARC_mut,SARC_mut$gene %in% path_immu_gene)
#47基因在SARC中的突变数据
SARC_m<-table(SARC_mut_47$sample,SARC_mut_47$gene)
SARC_mutat<-t(SARC_m)%>%as.matrix.data.frame()
rownames(SARC_mutat)<-colnames(SARC_m)
colnames(SARC_mutat)<-rownames(SARC_m)

SARC_l<-list()
for(i in 1:dim(SARC_m)[2]){
  
  x<-t(SARC_tpm2[colnames(SARC_m)[i],])
  SARC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  SARC_exp_gene$sample<-as.character(SARC_exp_gene$sample)
  
  x2<-t(SARC_mutat[colnames(SARC_m)[i],])
  SARC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  SARC_mut_gene$sample<-as.character(SARC_mut_gene$sample)
  
  SARC_union<-inner_join(SARC_exp_gene,SARC_mut_gene)
  SARC_union$cancer<-rep("SARC",dim(SARC_union)[1])
  SARC_union$gene<-rep(colnames(SARC_m)[i],dim(SARC_union)[1])
  SARC_union$mut_group<-ifelse(SARC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(SARC_union)[2]<-"mRNA"
  SARC_l<-rbind(SARC_union,SARC_l)
}
#############
LAML_tpm<-read.table("F:/workplace/LAML/LAML_readcount.genes.tpm.txt")

LAML_tpm<-(log2(LAML_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LAML_tpm)<-str_replace_all(colnames(LAML_tpm),"[.*]",replacement = "-")

LAML_gene<-intersect(rownames(LAML_tpm),path_immu_gene) 
LAML_tpm2<-LAML_tpm[LAML_gene,]  #免疫衰老基因的表达谱

########突变分组
LAML<-fread("f:/workplace/LAML/LAML_mc3.txt")%>%as.data.frame()%>%distinct()
LAML_mut<-subset(LAML,effect!="Silent") #去除沉默突变,剩16种突变类型

LAML_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

LAML_mut_47<-subset(LAML_mut,LAML_mut$gene %in% path_immu_gene)
#47基因在LAML中的突变数据
LAML_m<-table(LAML_mut_47$sample,LAML_mut_47$gene)
LAML_mutat<-t(LAML_m)%>%as.matrix.data.frame()
rownames(LAML_mutat)<-colnames(LAML_m)
colnames(LAML_mutat)<-rownames(LAML_m)

LAML_l<-list()
for(i in 1:dim(LAML_m)[2]){
  
  x<-t(LAML_tpm2[colnames(LAML_m)[i],])
  LAML_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LAML_exp_gene$sample<-as.character(LAML_exp_gene$sample)
  
  x2<-t(LAML_mutat[colnames(LAML_m)[i],])
  LAML_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  LAML_mut_gene$sample<-as.character(LAML_mut_gene$sample)
  
  LAML_union<-inner_join(LAML_exp_gene,LAML_mut_gene)
  LAML_union$cancer<-rep("LAML",dim(LAML_union)[1])
  LAML_union$gene<-rep(colnames(LAML_m)[i],dim(LAML_union)[1])
  LAML_union$mut_group<-ifelse(LAML_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(LAML_union)[2]<-"mRNA"
  LAML_l<-rbind(LAML_union,LAML_l)
}
##########
PAAD_tpm<-read.table("F:/workplace/PAAD/PAAD_readcount.genes.tpm.txt")

PAAD_tpm<-(log2(PAAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(PAAD_tpm)<-str_replace_all(colnames(PAAD_tpm),"[.*]",replacement = "-")

PAAD_gene<-intersect(rownames(PAAD_tpm),path_immu_gene) 
PAAD_tpm2<-PAAD_tpm[PAAD_gene,]  #免疫衰老基因的表达谱
########突变分组
PAAD<-fread("f:/workplace/PAAD/PAAD_mc3.txt")%>%as.data.frame()%>%distinct()
PAAD_mut<-subset(PAAD,effect!="Silent") #去除沉默突变,剩16种突变类型

PAAD_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

PAAD_mut_47<-subset(PAAD_mut,PAAD_mut$gene %in% path_immu_gene)
#47基因在PAAD中的突变数据
PAAD_m<-table(PAAD_mut_47$sample,PAAD_mut_47$gene)
PAAD_mutat<-t(PAAD_m)%>%as.matrix.data.frame()
rownames(PAAD_mutat)<-colnames(PAAD_m)
colnames(PAAD_mutat)<-rownames(PAAD_m)

PAAD_l<-list()
for(i in 1:dim(PAAD_m)[2]){
  
  x<-t(PAAD_tpm2[colnames(PAAD_m)[i],])
  PAAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  PAAD_exp_gene$sample<-as.character(PAAD_exp_gene$sample)
  
  x2<-t(PAAD_mutat[colnames(PAAD_m)[i],])
  PAAD_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  PAAD_mut_gene$sample<-as.character(PAAD_mut_gene$sample)
  
  PAAD_union<-inner_join(PAAD_exp_gene,PAAD_mut_gene)
  PAAD_union$cancer<-rep("PAAD",dim(PAAD_union)[1])
  PAAD_union$gene<-rep(colnames(PAAD_m)[i],dim(PAAD_union)[1])
  PAAD_union$mut_group<-ifelse(PAAD_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(PAAD_union)[2]<-"mRNA"
  PAAD_l<-rbind(PAAD_union,PAAD_l)
}
##########
ESCA_tpm<-read.table("F:/workplace/ESCA/ESCA_readcount.genes.tpm.txt")

ESCA_tpm<-(log2(ESCA_tpm+1))  #表达谱变为log2(TPM+1)
colnames(ESCA_tpm)<-str_replace_all(colnames(ESCA_tpm),"[.*]",replacement = "-")

ESCA_gene<-intersect(rownames(ESCA_tpm),path_immu_gene) 
ESCA_tpm2<-ESCA_tpm[ESCA_gene,]  #免疫衰老基因的表达谱

########突变分组
ESCA<-fread("f:/workplace/ESCA/ESCA_mc3.txt")%>%as.data.frame()%>%distinct()
ESCA_mut<-subset(ESCA,effect!="Silent") #去除沉默突变,剩16种突变类型

ESCA_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

ESCA_mut_47<-subset(ESCA_mut,ESCA_mut$gene %in% path_immu_gene)
#47基因在ESCA中的突变数据
ESCA_m<-table(ESCA_mut_47$sample,ESCA_mut_47$gene)
ESCA_mutat<-t(ESCA_m)%>%as.matrix.data.frame()
rownames(ESCA_mutat)<-colnames(ESCA_m)
colnames(ESCA_mutat)<-rownames(ESCA_m)

ESCA_l<-list()
for(i in 1:dim(ESCA_m)[2]){
  
  x<-t(ESCA_tpm2[colnames(ESCA_m)[i],])
  ESCA_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  ESCA_exp_gene$sample<-as.character(ESCA_exp_gene$sample)
  
  x2<-t(ESCA_mutat[colnames(ESCA_m)[i],])
  ESCA_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  ESCA_mut_gene$sample<-as.character(ESCA_mut_gene$sample)
  
  ESCA_union<-inner_join(ESCA_exp_gene,ESCA_mut_gene)
  ESCA_union$cancer<-rep("ESCA",dim(ESCA_union)[1])
  ESCA_union$gene<-rep(colnames(ESCA_m)[i],dim(ESCA_union)[1])
  ESCA_union$mut_group<-ifelse(ESCA_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(ESCA_union)[2]<-"mRNA"
  ESCA_l<-rbind(ESCA_union,ESCA_l)
}
#########
TGCT_tpm<-read.table("F:/workplace/TGCT/TGCT_readcount.genes.tpm.txt")

TGCT_tpm<-(log2(TGCT_tpm+1))  #表达谱变为log2(TPM+1)
colnames(TGCT_tpm)<-str_replace_all(colnames(TGCT_tpm),"[.*]",replacement = "-")

TGCT_gene<-intersect(rownames(TGCT_tpm),path_immu_gene) 
TGCT_tpm2<-TGCT_tpm[TGCT_gene,]  #免疫衰老基因的表达谱

########突变分组
TGCT<-fread("f:/workplace/TGCT/TGCT_mc3.txt")%>%as.data.frame()%>%distinct()
TGCT_mut<-subset(TGCT,effect!="Silent") #去除沉默突变,剩16种突变类型

TGCT_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

TGCT_mut_47<-subset(TGCT_mut,TGCT_mut$gene %in% path_immu_gene)
#47基因在TGCT中的突变数据
TGCT_m<-table(TGCT_mut_47$sample,TGCT_mut_47$gene)
TGCT_mutat<-t(TGCT_m)%>%as.matrix.data.frame()
rownames(TGCT_mutat)<-colnames(TGCT_m)
colnames(TGCT_mutat)<-rownames(TGCT_m)

TGCT_l<-list()
for(i in 1:dim(TGCT_m)[2]){
  
  x<-t(TGCT_tpm2[colnames(TGCT_m)[i],])
  TGCT_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  TGCT_exp_gene$sample<-as.character(TGCT_exp_gene$sample)
  
  x2<-t(TGCT_mutat[colnames(TGCT_m)[i],])
  TGCT_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  TGCT_mut_gene$sample<-as.character(TGCT_mut_gene$sample)
  
  TGCT_union<-inner_join(TGCT_exp_gene,TGCT_mut_gene)
  TGCT_union$cancer<-rep("TGCT",dim(TGCT_union)[1])
  TGCT_union$gene<-rep(colnames(TGCT_m)[i],dim(TGCT_union)[1])
  TGCT_union$mut_group<-ifelse(TGCT_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(TGCT_union)[2]<-"mRNA"
  TGCT_l<-rbind(TGCT_union,TGCT_l)
}
##########
THYM_tpm<-read.table("F:/workplace/THYM/THYM_readcount.genes.tpm.txt")

THYM_tpm<-(log2(THYM_tpm+1))  #表达谱变为log2(TPM+1)
colnames(THYM_tpm)<-str_replace_all(colnames(THYM_tpm),"[.*]",replacement = "-")

THYM_gene<-intersect(rownames(THYM_tpm),path_immu_gene) 
THYM_tpm2<-THYM_tpm[THYM_gene,]  #免疫衰老基因的表达谱

########突变分组
THYM<-fread("f:/workplace/THYM/THYM_mc3.txt")%>%as.data.frame()%>%distinct()
THYM_mut<-subset(THYM,effect!="Silent") #去除沉默突变,剩16种突变类型

THYM_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

THYM_mut_47<-subset(THYM_mut,THYM_mut$gene %in% path_immu_gene)
#47基因在THYM中的突变数据
THYM_m<-table(THYM_mut_47$sample,THYM_mut_47$gene)
THYM_mutat<-t(THYM_m)%>%as.matrix.data.frame()
rownames(THYM_mutat)<-colnames(THYM_m)
colnames(THYM_mutat)<-rownames(THYM_m)

THYM_l<-list()
for(i in 1:dim(THYM_m)[2]){
  
  x<-t(THYM_tpm2[colnames(THYM_m)[i],])
  THYM_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  THYM_exp_gene$sample<-as.character(THYM_exp_gene$sample)
  
  x2<-t(THYM_mutat[colnames(THYM_m)[i],])
  THYM_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  THYM_mut_gene$sample<-as.character(THYM_mut_gene$sample)
  
  THYM_union<-inner_join(THYM_exp_gene,THYM_mut_gene)
  THYM_union$cancer<-rep("THYM",dim(THYM_union)[1])
  THYM_union$gene<-rep(colnames(THYM_m)[i],dim(THYM_union)[1])
  THYM_union$mut_group<-ifelse(THYM_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(THYM_union)[2]<-"mRNA"
  THYM_l<-rbind(THYM_union,THYM_l)
}
###########3
MESO_tpm<-read.table("F:/workplace/MESO/MESO_readcount.genes.tpm.txt")

MESO_tpm<-(log2(MESO_tpm+1))  #表达谱变为log2(TPM+1)
colnames(MESO_tpm)<-str_replace_all(colnames(MESO_tpm),"[.*]",replacement = "-")

MESO_gene<-intersect(rownames(MESO_tpm),path_immu_gene) 
MESO_tpm2<-MESO_tpm[MESO_gene,]  #免疫衰老基因的表达谱

########突变分组
MESO<-fread("f:/workplace/MESO/MESO_mc3.txt")%>%as.data.frame()%>%distinct()
MESO_mut<-subset(MESO,effect!="Silent") #去除沉默突变,剩16种突变类型

MESO_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

MESO_mut_47<-subset(MESO_mut,MESO_mut$gene %in% path_immu_gene)
#47基因在MESO中的突变数据
MESO_m<-table(MESO_mut_47$sample,MESO_mut_47$gene)
MESO_mutat<-t(MESO_m)%>%as.matrix.data.frame()
rownames(MESO_mutat)<-colnames(MESO_m)
colnames(MESO_mutat)<-rownames(MESO_m)

MESO_l<-list()
for(i in 1:dim(MESO_m)[2]){
  
  x<-t(MESO_tpm2[colnames(MESO_m)[i],])
  MESO_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  MESO_exp_gene$sample<-as.character(MESO_exp_gene$sample)
  
  x2<-t(MESO_mutat[colnames(MESO_m)[i],])
  MESO_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  MESO_mut_gene$sample<-as.character(MESO_mut_gene$sample)
  
  MESO_union<-inner_join(MESO_exp_gene,MESO_mut_gene)
  MESO_union$cancer<-rep("MESO",dim(MESO_union)[1])
  MESO_union$gene<-rep(colnames(MESO_m)[i],dim(MESO_union)[1])
  MESO_union$mut_group<-ifelse(MESO_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(MESO_union)[2]<-"mRNA"
  MESO_l<-rbind(MESO_union,MESO_l)
}
###########
UVM_tpm<-read.table("F:/workplace/UVM/UVM_readcount.genes.tpm.txt")

UVM_tpm<-(log2(UVM_tpm+1))  #表达谱变为log2(TPM+1)
colnames(UVM_tpm)<-str_replace_all(colnames(UVM_tpm),"[.*]",replacement = "-")

UVM_gene<-intersect(rownames(UVM_tpm),path_immu_gene) 
UVM_tpm2<-UVM_tpm[UVM_gene,]  #免疫衰老基因的表达谱
########突变分组
UVM<-fread("f:/workplace/UVM/UVM_mc3.txt")%>%as.data.frame()%>%distinct()
UVM_mut<-subset(UVM,effect!="Silent") #去除沉默突变,剩16种突变类型

UVM_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

UVM_mut_47<-subset(UVM_mut,UVM_mut$gene %in% path_immu_gene)
#47基因在UVM中的突变数据
UVM_m<-table(UVM_mut_47$sample,UVM_mut_47$gene)
UVM_mutat<-t(UVM_m)%>%as.matrix.data.frame()
rownames(UVM_mutat)<-colnames(UVM_m)
colnames(UVM_mutat)<-rownames(UVM_m)
#####UVM没有免疫衰老基因突变的分组
UVM_l<-data.frame()
for(i in 1:dim(UVM_m)[2]){
  
  x<-t(UVM_tpm2[colnames(UVM_m)[i],])
  UVM_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  UVM_exp_gene$sample<-as.character(UVM_exp_gene$sample)
  
  x2<-t(UVM_mutat[colnames(UVM_m)[i],])
  UVM_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  UVM_mut_gene$sample<-as.character(UVM_mut_gene$sample)
  
  UVM_union<-inner_join(UVM_exp_gene,UVM_mut_gene)
  UVM_union$cancer<-rep("UVM",dim(UVM_union)[1])
  UVM_union$gene<-rep(colnames(UVM_m)[i],dim(UVM_union)[1])
  UVM_union$mut_group<-ifelse(UVM_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(UVM_union)[2]<-"mRNA"
  UVM_l<-rbind(UVM_union,UVM_l)
}
###########
DLBC_tpm<-read.table("F:/workplace/DLBC/DLBC_readcount.genes.tpm.txt")

DLBC_tpm<-(log2(DLBC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(DLBC_tpm)<-str_replace_all(colnames(DLBC_tpm),"[.*]",replacement = "-")

DLBC_gene<-intersect(rownames(DLBC_tpm),path_immu_gene) 
DLBC_tpm2<-DLBC_tpm[DLBC_gene,]  #免疫衰老基因的表达谱

########突变分组
DLBC<-fread("f:/workplace/DLBC/DLBC_mc3.txt")%>%as.data.frame()%>%distinct()
DLBC_mut<-subset(DLBC,effect!="Silent") #去除沉默突变,剩16种突变类型

DLBC_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

DLBC_mut_47<-subset(DLBC_mut,DLBC_mut$gene %in% path_immu_gene)
#47基因在DLBC中的突变数据
DLBC_m<-table(DLBC_mut_47$sample,DLBC_mut_47$gene)
DLBC_mutat<-t(DLBC_m)%>%as.matrix.data.frame()
rownames(DLBC_mutat)<-colnames(DLBC_m)
colnames(DLBC_mutat)<-rownames(DLBC_m)

DLBC_l<-list()
for(i in 1:dim(DLBC_m)[2]){
  
  x<-t(DLBC_tpm2[colnames(DLBC_m)[i],])
  DLBC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  DLBC_exp_gene$sample<-as.character(DLBC_exp_gene$sample)
  
  x2<-t(DLBC_mutat[colnames(DLBC_m)[i],])
  DLBC_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  DLBC_mut_gene$sample<-as.character(DLBC_mut_gene$sample)
  
  DLBC_union<-inner_join(DLBC_exp_gene,DLBC_mut_gene)
  DLBC_union$cancer<-rep("DLBC",dim(DLBC_union)[1])
  DLBC_union$gene<-rep(colnames(DLBC_m)[i],dim(DLBC_union)[1])
  DLBC_union$mut_group<-ifelse(DLBC_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(DLBC_union)[2]<-"mRNA"
  DLBC_l<-rbind(DLBC_union,DLBC_l)
}
###########
CHOL_tpm<-read.table("F:/workplace/CHOL/CHOL_readcount.genes.tpm.txt")

CHOL_tpm<-(log2(CHOL_tpm+1))  #表达谱变为log2(TPM+1)
colnames(CHOL_tpm)<-str_replace_all(colnames(CHOL_tpm),"[.*]",replacement = "-")

CHOL_gene<-intersect(rownames(CHOL_tpm),path_immu_gene) 
CHOL_tpm2<-CHOL_tpm[CHOL_gene,]  #免疫衰老基因的表达谱

########突变分组
CHOL<-fread("f:/workplace/CHOL/CHOL_mc3.txt")%>%as.data.frame()%>%distinct()
CHOL_mut<-subset(CHOL,effect!="Silent") #去除沉默突变,剩16种突变类型

CHOL_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

CHOL_mut_47<-subset(CHOL_mut,CHOL_mut$gene %in% path_immu_gene)
#47基因在CHOL中的突变数据
CHOL_m<-table(CHOL_mut_47$sample,CHOL_mut_47$gene)
CHOL_mutat<-t(CHOL_m)%>%as.matrix.data.frame()
rownames(CHOL_mutat)<-colnames(CHOL_m)
colnames(CHOL_mutat)<-rownames(CHOL_m)

CHOL_l<-list()
for(i in 1:dim(CHOL_m)[2]){
  
  x<-t(CHOL_tpm2[colnames(CHOL_m)[i],])
  CHOL_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  CHOL_exp_gene$sample<-as.character(CHOL_exp_gene$sample)
  
  x2<-t(CHOL_mutat[colnames(CHOL_m)[i],])
  CHOL_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  CHOL_mut_gene$sample<-as.character(CHOL_mut_gene$sample)
  
  CHOL_union<-inner_join(CHOL_exp_gene,CHOL_mut_gene)
  CHOL_union$cancer<-rep("CHOL",dim(CHOL_union)[1])
  CHOL_union$gene<-rep(colnames(CHOL_m)[i],dim(CHOL_union)[1])
  CHOL_union$mut_group<-ifelse(CHOL_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(CHOL_union)[2]<-"mRNA"
  CHOL_l<-rbind(CHOL_union,CHOL_l)
}
############
READ_tpm<-read.table("F:/workplace/READ/READ_readcount.genes.tpm.txt")

READ_tpm<-(log2(READ_tpm+1))  #表达谱变为log2(TPM+1)
colnames(READ_tpm)<-str_replace_all(colnames(READ_tpm),"[.*]",replacement = "-")

READ_gene<-intersect(rownames(READ_tpm),path_immu_gene) 
READ_tpm2<-READ_tpm[READ_gene,]  #免疫衰老基因的表达谱

########突变分组
READ<-fread("f:/workplace/READ/READ_mc3.txt")%>%as.data.frame()%>%distinct()
READ_mut<-subset(READ,effect!="Silent") #去除沉默突变,剩16种突变类型

READ_mut$sample%>%unique()%>%length()
#统计突变数据的样本数

READ_mut_47<-subset(READ_mut,READ_mut$gene %in% path_immu_gene)
#47基因在READ中的突变数据
READ_m<-table(READ_mut_47$sample,READ_mut_47$gene)
READ_mutat<-t(READ_m)%>%as.matrix.data.frame()
rownames(READ_mutat)<-colnames(READ_m)
colnames(READ_mutat)<-rownames(READ_m)

READ_l<-data.frame()
for(i in 1:dim(READ_m)[2]){
  
  x<-t(READ_tpm2[colnames(READ_m)[i],])
  READ_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  READ_exp_gene$sample<-as.character(READ_exp_gene$sample)
  
  x2<-t(READ_mutat[colnames(READ_m)[i],])
  READ_mut_gene<-data.frame(sample=colnames(x2),mutation=x2[1,])
  READ_mut_gene$sample<-as.character(READ_mut_gene$sample)
  
  READ_union<-inner_join(READ_exp_gene,READ_mut_gene)
  READ_union$cancer<-rep("READ",dim(READ_union)[1])
  READ_union$gene<-rep(colnames(READ_m)[i],dim(READ_union)[1])
  READ_union$mut_group<-ifelse(READ_union$mutation>0,"MUT","Non-MUT")%>%as.factor()
  colnames(READ_union)[2]<-"mRNA"
  READ_l<-rbind(READ_union,READ_l)
}

# 47个基因boxplot的输入文件 -------------------------------------------------------

box_mut<-rbind(LUAD_l,BLCA_l,HNSC_l,KIRC_l,KIRP_l,LIHC_l,LUSC_l,THCA_l,COAD_l,GBM_l,KICH_l,LGG_l,BRCA_l,READ_l,UCS_l,UCEC_l,OV_l,PRAD_l,STAD_l,SKCM_l,CESC_l,ACC_l,PCPG_l,SARC_l,LAML_l,PAAD_l,ESCA_l,TGCT_l,THYM_l,MESO_l,UVM_l,DLBC_l,CHOL_l)

write.csv(box_mut,"f:/workplace/结果2/mRNA-boxplot/突变分组的显著基因-boxplot/按突变分组的mRNA表达数据.csv",row.names=F)
mRNA_mut<-c()

for(i in 1:length(path_immu_gene)){
    gene_mut<-subset(box_mut,box_mut$gene%in%path_immu_gene[[i]])
    mRNA_diff<-compare_means(mRNA ~ mut_group, data = gene_mut, method = "t.test")
    mRNA_mut[i]<-mRNA_diff$p.adj
}

names(mRNA_mut)<-path_immu_gene
mRNA_mut[which(mRNA_mut<0.05)]
#有8个基因的表达值差异显著，样本按是否mutation分组

mut_effect_mRNA<-data.frame()
for(i in 1:length(path_immu_gene)){
  gene_mut<-subset(box_mut,box_mut$gene%in%path_immu_gene[[i]])
  mut_compare<-gene_mut %>%
    group_by(mut_group) %>%
    dplyr::summarise(mean_mRNA=mean(mRNA))
  
  effect<-ifelse(mut_compare$mean_mRNA[1]>mut_compare$mean_mRNA[2],"activate","inactivate")
  df<-data.frame(gene=path_immu_gene[[i]],cnv_effect=effect)
  mut_effect_mRNA<-rbind(mut_effect_mRNA,df)
}
mRNA_result<-data.frame(gene=names(mRNA_mut),diff_pvalue=mRNA_mut)
mRNA_mut_result<-inner_join(mRNA_result,mut_effect_mRNA)
#45个基因按mut，non-mut分组后，mRNA的差异，以及mut对MRNA的激活或失活统计

write.csv(mRNA_mut_result,"f:/workplace/结果2/mRNA-boxplot/突变分组的显著基因-boxplot/免疫衰老基因激活-失活统计.csv",row.names=F)

# 差异显著的8个基因-mRNA_boxplot ----------------------------------------------------

gene8<-mRNA_mut[which(mRNA_mut<0.05)]%>%names
df8<-subset(box_mut,box_mut$gene%in%gene8)

setwd("f:/workplace/结果2/mRNA-boxplot/突变分组的显著基因-boxplot/")
for(i in 1:length(gene8)){
  
  df<-subset(df8,df8$gene==gene8[i])
  pic<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
                 palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                    "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                    "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                    "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
  ggsave(pic, file=paste(gene8[i],".pdf",sep=""), width=8, height=8)
}

###
i=1
df<-subset(df8,df8$gene==gene8[i])
p1<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
               palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                  "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                  "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                  "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
i=2
df<-subset(df8,df8$gene==gene8[i])
p2<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
              palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                 "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                 "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                 "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
i=3
df<-subset(df8,df8$gene==gene8[i])
p3<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
              palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                 "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                 "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                 "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
i=4
df<-subset(df8,df8$gene==gene8[i])
p4<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
              palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                 "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                 "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                 "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
i=5
df<-subset(df8,df8$gene==gene8[i])
p5<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
              palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                 "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                 "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                 "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
i=6
df<-subset(df8,df8$gene==gene8[i])
p6<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
              palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                 "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                 "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                 "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
i=7
df<-subset(df8,df8$gene==gene8[i])
p7<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
              palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                 "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                 "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                 "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
i=8
df<-subset(df8,df8$gene==gene8[i])
p8<-ggboxplot(df, x = "mut_group", y = "mRNA",add.params = list(color = "cancer",size=0.1),add = "jitter",
              palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                 "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                 "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                 "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene8[i])+xlab("")+theme(legend.position = "none")
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8)

# 突变激活的基因boxplot ----------------------------------------------------------

gene_ac<-subset(mRNA_mut_result,mRNA_mut_result$cnv_effect=="activate")[,1]
df_ac<-subset(box_mut,box_mut$gene%in%gene_ac)

setwd("f:/workplace/结果2/mRNA-boxplot/突变分组的显著基因-boxplot/")
plot_list<-list()
for(i in 1:length(gene_ac)){
  
  df<-subset(df_ac,df_ac$gene==gene_ac[i])
  plot_list[[i]]<-ggboxplot(df, x = "mut_group", y = "mRNA",outlier.size = 0,add.params = list(color = "cancer",size=0.001,alpha=0.4),add = "jitter",
                 palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                    "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                    "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                    "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(label = "p.signif",method = "t.test",size=2.8)+ylab(paste(gene_ac[i],"expression"," "))+xlab("")+theme(axis.title.y = element_text(size = 6),axis.text.x = element_text(size = 6,color="black"),axis.text.y = element_text(size = 6,color="black") ,legend.position = "none")
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
          plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]],
          plot_list[[19]],plot_list[[20]],plot_list[[21]],plot_list[[22]])

# 突变失活的基因boxplot ----------------------------------------------------------

gene_ia<-subset(mRNA_mut_result,mRNA_mut_result$cnv_effect=="inactivate")[,1]
df_ia<-subset(box_mut,box_mut$gene%in%gene_ia)

setwd("f:/workplace/结果2/mRNA-boxplot/突变分组的显著基因-boxplot/")
plot_list<-list()
for(i in 1:length(gene_ia)){
  
  df<-subset(df_ia,df_ia$gene==gene_ia[i])
  plot_list[[i]]<-ggboxplot(df, x = "mut_group", y = "mRNA",outlier.size = 0,add.params = list(color = "cancer",size=0.001,alpha=0.4),add = "jitter",
                            palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                               "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                               "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                               "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(label = "p.signif",method = "t.test",size=2.8)+ylab(paste(gene_ac[i],"expression"," "))+xlab("")+theme(axis.title.y = element_text(size = 6),axis.text.x = element_text(size = 6,color="black"),axis.text.y = element_text(size = 6,color="black") ,legend.position = "none")
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
          plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]],
          plot_list[[19]],plot_list[[20]],plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]])


#############
df_ac<-box_mut
gene_ac<-mRNA_mut_result[,1]
plot_list<-list()
for(i in 1:length(gene_ac)){
  
  df<-subset(df_ac,df_ac$gene==gene_ac[i])
  plot_list[[i]]<-ggboxplot(df, x = "mut_group", y = "mRNA",outlier.size = 0,add.params = list(color = "cancer",size=0.001,alpha=0.4),add = "jitter",
                            palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                               "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                               "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                               "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(label = "p.signif",method = "t.test",size=2.8)+ylab(paste(gene_ac[i],"expression"," "))+xlab("")+theme(axis.title.y = element_text(size = 6),axis.text.x = element_text(size = 6,color="black"),axis.text.y = element_text(size = 6,color="black") ,legend.position = "none")
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
          plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]],
          plot_list[[19]],plot_list[[20]],plot_list[[21]],plot_list[[22]],plot_list[[23]])

plot_grid(plot_list[[24]],plot_list[[25]],plot_list[[26]],
          plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],plot_list[[31]],plot_list[[32]],
          plot_list[[33]],plot_list[[34]],plot_list[[35]],plot_list[[36]],plot_list[[37]],plot_list[[38]],
          plot_list[[39]],plot_list[[40]],plot_list[[41]],plot_list[[42]],plot_list[[43]],plot_list[[44]],
          plot_list[[45]],plot_list[[46]],plot_list[[47]])