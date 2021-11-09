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

# 47免疫衰老基因的mRNA表达谱及其CNV分组 ----------------------------------------------------------
############
BLCA_tpm<-read.table("F:/workplace/BLCA/BLCA_readcount.genes.tpm.txt")

BLCA_tpm<-(log2(BLCA_tpm+1))  #表达谱变为log2(TPM+1)
colnames(BLCA_tpm)<-str_replace_all(colnames(BLCA_tpm),"[.*]",replacement = "-")

BLCA_gene<-intersect(rownames(BLCA_tpm),path_immu_gene) 
BLCA_tpm2<-BLCA_tpm[BLCA_gene,]  #免疫衰老基因的表达谱

BLCA<-fread("f:/workplace/BLCA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(BLCA)<-BLCA$`Gene Symbol`
BLCA<-BLCA[,-1]

both_gene<-intersect(path_immu_gene,rownames(BLCA))

BLCA_CNV<-BLCA[both_gene,]
BLCA_cnv<-ifelse(BLCA_CNV<0,1,ifelse(BLCA_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
BLCA_cnv<-as.data.frame(BLCA_cnv)

BLCA_l<-list()

for(i in 1:length(path_immu_gene)){
  
 x<-t(BLCA_tpm2[path_immu_gene[i],])
 BLCA_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
 BLCA_exp_gene$sample<-as.character(BLCA_exp_gene$sample)

 x2<-t(BLCA_cnv[path_immu_gene[i],])
 BLCA_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
 BLCA_cnv_gene$sample<-as.character(BLCA_cnv_gene$sample)

 BLCA_union<-inner_join(BLCA_exp_gene,BLCA_cnv_gene)
 BLCA_union$cancer<-rep("BLCA",dim(BLCA_union)[1])
 BLCA_union$CNV_group<-ifelse(BLCA_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
 BLCA_l[[i]]<-BLCA_union
}
###########


LUAD_tpm<-read.table("F:/workplace/LUAD/LUAD_readcount.genes.tpm.txt")

LUAD_tpm<-(log2(LUAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LUAD_tpm)<-str_replace_all(colnames(LUAD_tpm),"[.*]",replacement = "-")

LUAD_gene<-intersect(rownames(LUAD_tpm),path_immu_gene) 
LUAD_tpm2<-LUAD_tpm[LUAD_gene,]  #免疫衰老基因的表达谱

LUAD<-fread("f:/workplace/LUAD/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(LUAD)<-LUAD$`Gene Symbol`
LUAD<-LUAD[,-1]

both_gene<-intersect(path_immu_gene,rownames(LUAD))

LUAD_CNV<-LUAD[both_gene,]
LUAD_cnv<-ifelse(LUAD_CNV<0,1,ifelse(LUAD_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
LUAD_cnv<-as.data.frame(LUAD_cnv)

LUAD_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(LUAD_tpm2[path_immu_gene[i],])
  LUAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LUAD_exp_gene$sample<-as.character(LUAD_exp_gene$sample)
  
  x2<-t(LUAD_cnv[path_immu_gene[i],])
  LUAD_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  LUAD_cnv_gene$sample<-as.character(LUAD_cnv_gene$sample)
  
  LUAD_union<-inner_join(LUAD_exp_gene,LUAD_cnv_gene)
  LUAD_union$cancer<-rep("LUAD",dim(LUAD_union)[1])
  LUAD_union$CNV_group<-ifelse(LUAD_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  LUAD_l[[i]]<-LUAD_union
}
###################


HNSC_tpm<-read.table("F:/workplace/HNSC/HNSC_readcount.genes.tpm.txt")

HNSC_tpm<-(log2(HNSC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(HNSC_tpm)<-str_replace_all(colnames(HNSC_tpm),"[.*]",replacement = "-")

HNSC_gene<-intersect(rownames(HNSC_tpm),path_immu_gene) 
HNSC_tpm2<-HNSC_tpm[HNSC_gene,]  #免疫衰老基因的表达谱

HNSC<-fread("f:/workplace/HNSC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(HNSC)<-HNSC$`Gene Symbol`
HNSC<-HNSC[,-1]

both_gene<-intersect(path_immu_gene,rownames(HNSC))

HNSC_CNV<-HNSC[both_gene,]
HNSC_cnv<-ifelse(HNSC_CNV<0,1,ifelse(HNSC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
HNSC_cnv<-as.data.frame(HNSC_cnv)

HNSC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(HNSC_tpm2[path_immu_gene[i],])
  HNSC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  HNSC_exp_gene$sample<-as.character(HNSC_exp_gene$sample)
  
  x2<-t(HNSC_cnv[path_immu_gene[i],])
  HNSC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  HNSC_cnv_gene$sample<-as.character(HNSC_cnv_gene$sample)
  
  HNSC_union<-inner_join(HNSC_exp_gene,HNSC_cnv_gene)
  HNSC_union$cancer<-rep("HNSC",dim(HNSC_union)[1])
  HNSC_union$CNV_group<-ifelse(HNSC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  HNSC_l[[i]]<-HNSC_union
}
########


KIRC_tpm<-read.table("F:/workplace/KIRC/KIRC_readcount.genes.tpm.txt")

KIRC_tpm<-(log2(KIRC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(KIRC_tpm)<-str_replace_all(colnames(KIRC_tpm),"[.*]",replacement = "-")

KIRC_gene<-intersect(rownames(KIRC_tpm),path_immu_gene) 
KIRC_tpm2<-KIRC_tpm[KIRC_gene,]  #免疫衰老基因的表达谱

KIRC<-fread("f:/workplace/KIRC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(KIRC)<-KIRC$`Gene Symbol`
KIRC<-KIRC[,-1]

both_gene<-intersect(path_immu_gene,rownames(KIRC))

KIRC_CNV<-KIRC[both_gene,]
KIRC_cnv<-ifelse(KIRC_CNV<0,1,ifelse(KIRC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
KIRC_cnv<-as.data.frame(KIRC_cnv)

KIRC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(KIRC_tpm2[path_immu_gene[i],])
  KIRC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  KIRC_exp_gene$sample<-as.character(KIRC_exp_gene$sample)
  
  x2<-t(KIRC_cnv[path_immu_gene[i],])
  KIRC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  KIRC_cnv_gene$sample<-as.character(KIRC_cnv_gene$sample)
  
  KIRC_union<-inner_join(KIRC_exp_gene,KIRC_cnv_gene)
  KIRC_union$cancer<-rep("KIRC",dim(KIRC_union)[1])
  KIRC_union$CNV_group<-ifelse(KIRC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  KIRC_l[[i]]<-KIRC_union
}
##############


KIRP_tpm<-read.table("F:/workplace/KIRP/KIRP_readcount.genes.tpm.txt")

KIRP_tpm<-(log2(KIRP_tpm+1))  #表达谱变为log2(TPM+1)
colnames(KIRP_tpm)<-str_replace_all(colnames(KIRP_tpm),"[.*]",replacement = "-")

KIRP_gene<-intersect(rownames(KIRP_tpm),path_immu_gene) 
KIRP_tpm2<-KIRP_tpm[KIRP_gene,]  #免疫衰老基因的表达谱

KIRP<-fread("f:/workplace/KIRP/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(KIRP)<-KIRP$`Gene Symbol`
KIRP<-KIRP[,-1]

both_gene<-intersect(path_immu_gene,rownames(KIRP))

KIRP_CNV<-KIRP[both_gene,]
KIRP_cnv<-ifelse(KIRP_CNV<0,1,ifelse(KIRP_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
KIRP_cnv<-as.data.frame(KIRP_cnv)

KIRP_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(KIRP_tpm2[path_immu_gene[i],])
  KIRP_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  KIRP_exp_gene$sample<-as.character(KIRP_exp_gene$sample)
  
  x2<-t(KIRP_cnv[path_immu_gene[i],])
  KIRP_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  KIRP_cnv_gene$sample<-as.character(KIRP_cnv_gene$sample)
  
  KIRP_union<-inner_join(KIRP_exp_gene,KIRP_cnv_gene)
  KIRP_union$cancer<-rep("KIRP",dim(KIRP_union)[1])
  KIRP_union$CNV_group<-ifelse(KIRP_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  KIRP_l[[i]]<-KIRP_union
}
########


LIHC_tpm<-read.table("F:/workplace/LIHC/LIHC_readcount.genes.tpm.txt")

LIHC_tpm<-(log2(LIHC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LIHC_tpm)<-str_replace_all(colnames(LIHC_tpm),"[.*]",replacement = "-")

LIHC_gene<-intersect(rownames(LIHC_tpm),path_immu_gene) 
LIHC_tpm2<-LIHC_tpm[LIHC_gene,]  #免疫衰老基因的表达谱

LIHC<-fread("f:/workplace/LIHC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(LIHC)<-LIHC$`Gene Symbol`
LIHC<-LIHC[,-1]

both_gene<-intersect(path_immu_gene,rownames(LIHC))

LIHC_CNV<-LIHC[both_gene,]
LIHC_cnv<-ifelse(LIHC_CNV<0,1,ifelse(LIHC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
LIHC_cnv<-as.data.frame(LIHC_cnv)

LIHC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(LIHC_tpm2[path_immu_gene[i],])
  LIHC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LIHC_exp_gene$sample<-as.character(LIHC_exp_gene$sample)
  
  x2<-t(LIHC_cnv[path_immu_gene[i],])
  LIHC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  LIHC_cnv_gene$sample<-as.character(LIHC_cnv_gene$sample)
  
  LIHC_union<-inner_join(LIHC_exp_gene,LIHC_cnv_gene)
  LIHC_union$cancer<-rep("LIHC",dim(LIHC_union)[1])
  LIHC_union$CNV_group<-ifelse(LIHC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  LIHC_l[[i]]<-LIHC_union
}
##############


LUSC_tpm<-read.table("F:/workplace/LUSC/LUSC_readcount.genes.tpm.txt")

LUSC_tpm<-(log2(LUSC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LUSC_tpm)<-str_replace_all(colnames(LUSC_tpm),"[.*]",replacement = "-")

LUSC_gene<-intersect(rownames(LUSC_tpm),path_immu_gene) 
LUSC_tpm2<-LUSC_tpm[LUSC_gene,]  #免疫衰老基因的表达谱

LUSC<-fread("f:/workplace/LUSC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(LUSC)<-LUSC$`Gene Symbol`
LUSC<-LUSC[,-1]

both_gene<-intersect(path_immu_gene,rownames(LUSC))

LUSC_CNV<-LUSC[both_gene,]
LUSC_cnv<-ifelse(LUSC_CNV<0,1,ifelse(LUSC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
LUSC_cnv<-as.data.frame(LUSC_cnv)

LUSC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(LUSC_tpm2[path_immu_gene[i],])
  LUSC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LUSC_exp_gene$sample<-as.character(LUSC_exp_gene$sample)
  
  x2<-t(LUSC_cnv[path_immu_gene[i],])
  LUSC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  LUSC_cnv_gene$sample<-as.character(LUSC_cnv_gene$sample)
  
  LUSC_union<-inner_join(LUSC_exp_gene,LUSC_cnv_gene)
  LUSC_union$cancer<-rep("LUSC",dim(LUSC_union)[1])
  LUSC_union$CNV_group<-ifelse(LUSC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  LUSC_l[[i]]<-LUSC_union
}
##########


THCA_tpm<-read.table("F:/workplace/THCA/THCA_readcount.genes.tpm.txt")

THCA_tpm<-(log2(THCA_tpm+1))  #表达谱变为log2(TPM+1)
colnames(THCA_tpm)<-str_replace_all(colnames(THCA_tpm),"[.*]",replacement = "-")

THCA_gene<-intersect(rownames(THCA_tpm),path_immu_gene) 
THCA_tpm2<-THCA_tpm[THCA_gene,]  #免疫衰老基因的表达谱

THCA<-fread("f:/workplace/THCA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(THCA)<-THCA$`Gene Symbol`
THCA<-THCA[,-1]

both_gene<-intersect(path_immu_gene,rownames(THCA))

THCA_CNV<-THCA[both_gene,]
THCA_cnv<-ifelse(THCA_CNV<0,1,ifelse(THCA_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
THCA_cnv<-as.data.frame(THCA_cnv)

THCA_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(THCA_tpm2[path_immu_gene[i],])
  THCA_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  THCA_exp_gene$sample<-as.character(THCA_exp_gene$sample)
  
  x2<-t(THCA_cnv[path_immu_gene[i],])
  THCA_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  THCA_cnv_gene$sample<-as.character(THCA_cnv_gene$sample)
  
  THCA_union<-inner_join(THCA_exp_gene,THCA_cnv_gene)
  THCA_union$cancer<-rep("THCA",dim(THCA_union)[1])
  THCA_union$CNV_group<-ifelse(THCA_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  THCA_l[[i]]<-THCA_union
}
######


COAD_tpm<-read.table("F:/workplace/COAD/COAD_readcount.genes.tpm.txt")

COAD_tpm<-(log2(COAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(COAD_tpm)<-str_replace_all(colnames(COAD_tpm),"[.*]",replacement = "-")

COAD_gene<-intersect(rownames(COAD_tpm),path_immu_gene) 
COAD_tpm2<-COAD_tpm[COAD_gene,]  #免疫衰老基因的表达谱

COAD<-fread("f:/workplace/COAD/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(COAD)<-COAD$`Gene Symbol`
COAD<-COAD[,-1]

both_gene<-intersect(path_immu_gene,rownames(COAD))

COAD_CNV<-COAD[both_gene,]
COAD_cnv<-ifelse(COAD_CNV<0,1,ifelse(COAD_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
COAD_cnv<-as.data.frame(COAD_cnv)

COAD_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(COAD_tpm2[path_immu_gene[i],])
  COAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  COAD_exp_gene$sample<-as.character(COAD_exp_gene$sample)
  
  x2<-t(COAD_cnv[path_immu_gene[i],])
  COAD_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  COAD_cnv_gene$sample<-as.character(COAD_cnv_gene$sample)
  
  COAD_union<-inner_join(COAD_exp_gene,COAD_cnv_gene)
  COAD_union$cancer<-rep("COAD",dim(COAD_union)[1])
  COAD_union$CNV_group<-ifelse(COAD_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  COAD_l[[i]]<-COAD_union
}
###########


GBM_tpm<-read.table("F:/workplace/GBM/GBM_readcount.genes.tpm.txt")

GBM_tpm<-(log2(GBM_tpm+1))  #表达谱变为log2(TPM+1)
colnames(GBM_tpm)<-str_replace_all(colnames(GBM_tpm),"[.*]",replacement = "-")

GBM_gene<-intersect(rownames(GBM_tpm),path_immu_gene) 
GBM_tpm2<-GBM_tpm[GBM_gene,]  #免疫衰老基因的表达谱

GBM<-fread("f:/workplace/GBM/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(GBM)<-GBM$`Gene Symbol`
GBM<-GBM[,-1]

both_gene<-intersect(path_immu_gene,rownames(GBM))

GBM_CNV<-GBM[both_gene,]
GBM_cnv<-ifelse(GBM_CNV<0,1,ifelse(GBM_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
GBM_cnv<-as.data.frame(GBM_cnv)

GBM_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(GBM_tpm2[path_immu_gene[i],])
  GBM_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  GBM_exp_gene$sample<-as.character(GBM_exp_gene$sample)
  
  x2<-t(GBM_cnv[path_immu_gene[i],])
  GBM_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  GBM_cnv_gene$sample<-as.character(GBM_cnv_gene$sample)
  
  GBM_union<-inner_join(GBM_exp_gene,GBM_cnv_gene)
  GBM_union$cancer<-rep("GBM",dim(GBM_union)[1])
  GBM_union$CNV_group<-ifelse(GBM_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  GBM_l[[i]]<-GBM_union
}
#######


KICH_tpm<-read.table("F:/workplace/KICH/KICH_readcount.genes.tpm.txt")

KICH_tpm<-(log2(KICH_tpm+1))  #表达谱变为log2(TPM+1)
colnames(KICH_tpm)<-str_replace_all(colnames(KICH_tpm),"[.*]",replacement = "-")

KICH_gene<-intersect(rownames(KICH_tpm),path_immu_gene) 
KICH_tpm2<-KICH_tpm[KICH_gene,]  #免疫衰老基因的表达谱

KICH<-fread("f:/workplace/KICH/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(KICH)<-KICH$`Gene Symbol`
KICH<-KICH[,-1]

both_gene<-intersect(path_immu_gene,rownames(KICH))

KICH_CNV<-KICH[both_gene,]
KICH_cnv<-ifelse(KICH_CNV<0,1,ifelse(KICH_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
KICH_cnv<-as.data.frame(KICH_cnv)

KICH_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(KICH_tpm2[path_immu_gene[i],])
  KICH_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  KICH_exp_gene$sample<-as.character(KICH_exp_gene$sample)
  
  x2<-t(KICH_cnv[path_immu_gene[i],])
  KICH_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  KICH_cnv_gene$sample<-as.character(KICH_cnv_gene$sample)
  
  KICH_union<-inner_join(KICH_exp_gene,KICH_cnv_gene)
  KICH_union$cancer<-rep("KICH",dim(KICH_union)[1])
  KICH_union$CNV_group<-ifelse(KICH_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  KICH_l[[i]]<-KICH_union
}
##########
LGG_tpm<-read.table("F:/workplace/LGG/LGG_readcount.genes.tpm.txt")

LGG_tpm<-(log2(LGG_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LGG_tpm)<-str_replace_all(colnames(LGG_tpm),"[.*]",replacement = "-")

LGG_gene<-intersect(rownames(LGG_tpm),path_immu_gene) 
LGG_tpm2<-LGG_tpm[LGG_gene,]  #免疫衰老基因的表达谱

LGG<-fread("f:/workplace/LGG/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(LGG)<-LGG$`Gene Symbol`
LGG<-LGG[,-1]

both_gene<-intersect(path_immu_gene,rownames(LGG))

LGG_CNV<-LGG[both_gene,]
LGG_cnv<-ifelse(LGG_CNV<0,1,ifelse(LGG_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
LGG_cnv<-as.data.frame(LGG_cnv)

LGG_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(LGG_tpm2[path_immu_gene[i],])
  LGG_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LGG_exp_gene$sample<-as.character(LGG_exp_gene$sample)
  
  x2<-t(LGG_cnv[path_immu_gene[i],])
  LGG_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  LGG_cnv_gene$sample<-as.character(LGG_cnv_gene$sample)
  
  LGG_union<-inner_join(LGG_exp_gene,LGG_cnv_gene)
  LGG_union$cancer<-rep("LGG",dim(LGG_union)[1])
  LGG_union$CNV_group<-ifelse(LGG_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  LGG_l[[i]]<-LGG_union
}
#########


BRCA_tpm<-read.table("F:/workplace/BRCA/BRCA_readcount.genes.tpm.txt")

BRCA_tpm<-(log2(BRCA_tpm+1))  #表达谱变为log2(TPM+1)
colnames(BRCA_tpm)<-str_replace_all(colnames(BRCA_tpm),"[.*]",replacement = "-")

BRCA_gene<-intersect(rownames(BRCA_tpm),path_immu_gene) 
BRCA_tpm2<-BRCA_tpm[BRCA_gene,]  #免疫衰老基因的表达谱

BRCA<-fread("f:/workplace/BRCA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(BRCA)<-BRCA$`Gene Symbol`
BRCA<-BRCA[,-1]

both_gene<-intersect(path_immu_gene,rownames(BRCA))

BRCA_CNV<-BRCA[both_gene,]
BRCA_cnv<-ifelse(BRCA_CNV<0,1,ifelse(BRCA_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
BRCA_cnv<-as.data.frame(BRCA_cnv)

BRCA_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(BRCA_tpm2[path_immu_gene[i],])
  BRCA_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  BRCA_exp_gene$sample<-as.character(BRCA_exp_gene$sample)
  
  x2<-t(BRCA_cnv[path_immu_gene[i],])
  BRCA_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  BRCA_cnv_gene$sample<-as.character(BRCA_cnv_gene$sample)
  
  BRCA_union<-inner_join(BRCA_exp_gene,BRCA_cnv_gene)
  BRCA_union$cancer<-rep("BRCA",dim(BRCA_union)[1])
  BRCA_union$CNV_group<-ifelse(BRCA_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  BRCA_l[[i]]<-BRCA_union
}
##########


UCS_tpm<-read.table("F:/workplace/UCS/UCS_readcount.genes.tpm.txt")

UCS_tpm<-(log2(UCS_tpm+1))  #表达谱变为log2(TPM+1)
colnames(UCS_tpm)<-str_replace_all(colnames(UCS_tpm),"[.*]",replacement = "-")

UCS_gene<-intersect(rownames(UCS_tpm),path_immu_gene) 
UCS_tpm2<-UCS_tpm[UCS_gene,]  #免疫衰老基因的表达谱

UCS<-fread("f:/workplace/UCS/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(UCS)<-UCS$`Gene Symbol`
UCS<-UCS[,-1]

both_gene<-intersect(path_immu_gene,rownames(UCS))

UCS_CNV<-UCS[both_gene,]
UCS_cnv<-ifelse(UCS_CNV<0,1,ifelse(UCS_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
UCS_cnv<-as.data.frame(UCS_cnv)

UCS_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(UCS_tpm2[path_immu_gene[i],])
  UCS_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  UCS_exp_gene$sample<-as.character(UCS_exp_gene$sample)
  
  x2<-t(UCS_cnv[path_immu_gene[i],])
  UCS_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  UCS_cnv_gene$sample<-as.character(UCS_cnv_gene$sample)
  
  UCS_union<-inner_join(UCS_exp_gene,UCS_cnv_gene)
  UCS_union$cancer<-rep("UCS",dim(UCS_union)[1])
  UCS_union$CNV_group<-ifelse(UCS_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  UCS_l[[i]]<-UCS_union
}
#########


UCEC_tpm<-read.table("F:/workplace/UCEC/UCEC_readcount.genes.tpm.txt")

UCEC_tpm<-(log2(UCEC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(UCEC_tpm)<-str_replace_all(colnames(UCEC_tpm),"[.*]",replacement = "-")

UCEC_gene<-intersect(rownames(UCEC_tpm),path_immu_gene) 
UCEC_tpm2<-UCEC_tpm[UCEC_gene,]  #免疫衰老基因的表达谱

UCEC<-fread("f:/workplace/UCEC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(UCEC)<-UCEC$`Gene Symbol`
UCEC<-UCEC[,-1]

both_gene<-intersect(path_immu_gene,rownames(UCEC))

UCEC_CNV<-UCEC[both_gene,]
UCEC_cnv<-ifelse(UCEC_CNV<0,1,ifelse(UCEC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
UCEC_cnv<-as.data.frame(UCEC_cnv)

UCEC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(UCEC_tpm2[path_immu_gene[i],])
  UCEC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  UCEC_exp_gene$sample<-as.character(UCEC_exp_gene$sample)
  
  x2<-t(UCEC_cnv[path_immu_gene[i],])
  UCEC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  UCEC_cnv_gene$sample<-as.character(UCEC_cnv_gene$sample)
  
  UCEC_union<-inner_join(UCEC_exp_gene,UCEC_cnv_gene)
  UCEC_union$cancer<-rep("UCEC",dim(UCEC_union)[1])
  UCEC_union$CNV_group<-ifelse(UCEC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  UCEC_l[[i]]<-UCEC_union
}
##########


OV_tpm<-read.table("F:/workplace/OV/OV_readcount.genes.tpm.txt")

OV_tpm<-(log2(OV_tpm+1))  #表达谱变为log2(TPM+1)
colnames(OV_tpm)<-str_replace_all(colnames(OV_tpm),"[.*]",replacement = "-")

OV_gene<-intersect(rownames(OV_tpm),path_immu_gene) 
OV_tpm2<-OV_tpm[OV_gene,]  #免疫衰老基因的表达谱

OV<-fread("f:/workplace/OV/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(OV)<-OV$`Gene Symbol`
OV<-OV[,-1]

both_gene<-intersect(path_immu_gene,rownames(OV))

OV_CNV<-OV[both_gene,]
OV_cnv<-ifelse(OV_CNV<0,1,ifelse(OV_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
OV_cnv<-as.data.frame(OV_cnv)

OV_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(OV_tpm2[path_immu_gene[i],])
  OV_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  OV_exp_gene$sample<-as.character(OV_exp_gene$sample)
  
  x2<-t(OV_cnv[path_immu_gene[i],])
  OV_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  OV_cnv_gene$sample<-as.character(OV_cnv_gene$sample)
  
  OV_union<-inner_join(OV_exp_gene,OV_cnv_gene)
  OV_union$cancer<-rep("OV",dim(OV_union)[1])
  OV_union$CNV_group<-ifelse(OV_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  OV_l[[i]]<-OV_union
}
###########


PRAD_tpm<-read.table("F:/workplace/PRAD/PRAD_readcount.genes.tpm.txt")

PRAD_tpm<-(log2(PRAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(PRAD_tpm)<-str_replace_all(colnames(PRAD_tpm),"[.*]",replacement = "-")

PRAD_gene<-intersect(rownames(PRAD_tpm),path_immu_gene) 
PRAD_tpm2<-PRAD_tpm[PRAD_gene,]  #免疫衰老基因的表达谱

PRAD<-fread("f:/workplace/PRAD/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(PRAD)<-PRAD$`Gene Symbol`
PRAD<-PRAD[,-1]

both_gene<-intersect(path_immu_gene,rownames(PRAD))

PRAD_CNV<-PRAD[both_gene,]
PRAD_cnv<-ifelse(PRAD_CNV<0,1,ifelse(PRAD_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
PRAD_cnv<-as.data.frame(PRAD_cnv)

PRAD_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(PRAD_tpm2[path_immu_gene[i],])
  PRAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  PRAD_exp_gene$sample<-as.character(PRAD_exp_gene$sample)
  
  x2<-t(PRAD_cnv[path_immu_gene[i],])
  PRAD_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  PRAD_cnv_gene$sample<-as.character(PRAD_cnv_gene$sample)
  
  PRAD_union<-inner_join(PRAD_exp_gene,PRAD_cnv_gene)
  PRAD_union$cancer<-rep("PRAD",dim(PRAD_union)[1])
  PRAD_union$CNV_group<-ifelse(PRAD_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  PRAD_l[[i]]<-PRAD_union
}
##########


STAD_tpm<-read.table("F:/workplace/STAD/STAD_readcount.genes.tpm.txt")

STAD_tpm<-(log2(STAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(STAD_tpm)<-str_replace_all(colnames(STAD_tpm),"[.*]",replacement = "-")

STAD_gene<-intersect(rownames(STAD_tpm),path_immu_gene) 
STAD_tpm2<-STAD_tpm[STAD_gene,]  #免疫衰老基因的表达谱

STAD<-fread("f:/workplace/STAD/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(STAD)<-STAD$`Gene Symbol`
STAD<-STAD[,-1]

both_gene<-intersect(path_immu_gene,rownames(STAD))

STAD_CNV<-STAD[both_gene,]
STAD_cnv<-ifelse(STAD_CNV<0,1,ifelse(STAD_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
STAD_cnv<-as.data.frame(STAD_cnv)

STAD_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(STAD_tpm2[path_immu_gene[i],])
  STAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  STAD_exp_gene$sample<-as.character(STAD_exp_gene$sample)
  
  x2<-t(STAD_cnv[path_immu_gene[i],])
  STAD_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  STAD_cnv_gene$sample<-as.character(STAD_cnv_gene$sample)
  
  STAD_union<-inner_join(STAD_exp_gene,STAD_cnv_gene)
  STAD_union$cancer<-rep("STAD",dim(STAD_union)[1])
  STAD_union$CNV_group<-ifelse(STAD_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  STAD_l[[i]]<-STAD_union
}
###########


SKCM_tpm<-read.table("F:/workplace/SKCM/SKCM_readcount.genes.tpm.txt")

SKCM_tpm<-(log2(SKCM_tpm+1))  #表达谱变为log2(TPM+1)
colnames(SKCM_tpm)<-str_replace_all(colnames(SKCM_tpm),"[.*]",replacement = "-")

SKCM_gene<-intersect(rownames(SKCM_tpm),path_immu_gene) 
SKCM_tpm2<-SKCM_tpm[SKCM_gene,]  #免疫衰老基因的表达谱

SKCM<-fread("f:/workplace/SKCM/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(SKCM)<-SKCM$`Gene Symbol`
SKCM<-SKCM[,-1]

both_gene<-intersect(path_immu_gene,rownames(SKCM))

SKCM_CNV<-SKCM[both_gene,]
SKCM_cnv<-ifelse(SKCM_CNV<0,1,ifelse(SKCM_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
SKCM_cnv<-as.data.frame(SKCM_cnv)

SKCM_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(SKCM_tpm2[path_immu_gene[i],])
  SKCM_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  SKCM_exp_gene$sample<-as.character(SKCM_exp_gene$sample)
  
  x2<-t(SKCM_cnv[path_immu_gene[i],])
  SKCM_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  SKCM_cnv_gene$sample<-as.character(SKCM_cnv_gene$sample)
  
  SKCM_union<-inner_join(SKCM_exp_gene,SKCM_cnv_gene)
  SKCM_union$cancer<-rep("SKCM",dim(SKCM_union)[1])
  SKCM_union$CNV_group<-ifelse(SKCM_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  SKCM_l[[i]]<-SKCM_union
}
#########


CESC_tpm<-read.table("F:/workplace/CESC/CESC_readcount.genes.tpm.txt")

CESC_tpm<-(log2(CESC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(CESC_tpm)<-str_replace_all(colnames(CESC_tpm),"[.*]",replacement = "-")

CESC_gene<-intersect(rownames(CESC_tpm),path_immu_gene) 
CESC_tpm2<-CESC_tpm[CESC_gene,]  #免疫衰老基因的表达谱

CESC<-fread("f:/workplace/CESC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(CESC)<-CESC$`Gene Symbol`
CESC<-CESC[,-1]

both_gene<-intersect(path_immu_gene,rownames(CESC))

CESC_CNV<-CESC[both_gene,]
CESC_cnv<-ifelse(CESC_CNV<0,1,ifelse(CESC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
CESC_cnv<-as.data.frame(CESC_cnv)

CESC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(CESC_tpm2[path_immu_gene[i],])
  CESC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  CESC_exp_gene$sample<-as.character(CESC_exp_gene$sample)
  
  x2<-t(CESC_cnv[path_immu_gene[i],])
  CESC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  CESC_cnv_gene$sample<-as.character(CESC_cnv_gene$sample)
  
  CESC_union<-inner_join(CESC_exp_gene,CESC_cnv_gene)
  CESC_union$cancer<-rep("CESC",dim(CESC_union)[1])
  CESC_union$CNV_group<-ifelse(CESC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  CESC_l[[i]]<-CESC_union
}
########


ACC_tpm<-read.table("F:/workplace/ACC/ACC_readcount.genes.tpm.txt")

ACC_tpm<-(log2(ACC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(ACC_tpm)<-str_replace_all(colnames(ACC_tpm),"[.*]",replacement = "-")

ACC_gene<-intersect(rownames(ACC_tpm),path_immu_gene) 
ACC_tpm2<-ACC_tpm[ACC_gene,]  #免疫衰老基因的表达谱

ACC<-fread("f:/workplace/ACC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(ACC)<-ACC$`Gene Symbol`
ACC<-ACC[,-1]

both_gene<-intersect(path_immu_gene,rownames(ACC))

ACC_CNV<-ACC[both_gene,]
ACC_cnv<-ifelse(ACC_CNV<0,1,ifelse(ACC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
ACC_cnv<-as.data.frame(ACC_cnv)

ACC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(ACC_tpm2[path_immu_gene[i],])
  ACC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  ACC_exp_gene$sample<-as.character(ACC_exp_gene$sample)
  
  x2<-t(ACC_cnv[path_immu_gene[i],])
  ACC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  ACC_cnv_gene$sample<-as.character(ACC_cnv_gene$sample)
  
  ACC_union<-inner_join(ACC_exp_gene,ACC_cnv_gene)
  ACC_union$cancer<-rep("ACC",dim(ACC_union)[1])
  ACC_union$CNV_group<-ifelse(ACC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  ACC_l[[i]]<-ACC_union
}
#######


PCPG_tpm<-read.table("F:/workplace/PCPG/PCPG_readcount.genes.tpm.txt")

PCPG_tpm<-(log2(PCPG_tpm+1))  #表达谱变为log2(TPM+1)
colnames(PCPG_tpm)<-str_replace_all(colnames(PCPG_tpm),"[.*]",replacement = "-")

PCPG_gene<-intersect(rownames(PCPG_tpm),path_immu_gene) 
PCPG_tpm2<-PCPG_tpm[PCPG_gene,]  #免疫衰老基因的表达谱

PCPG<-fread("f:/workplace/PCPG/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(PCPG)<-PCPG$`Gene Symbol`
PCPG<-PCPG[,-1]

both_gene<-intersect(path_immu_gene,rownames(PCPG))

PCPG_CNV<-PCPG[both_gene,]
PCPG_cnv<-ifelse(PCPG_CNV<0,1,ifelse(PCPG_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
PCPG_cnv<-as.data.frame(PCPG_cnv)

PCPG_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(PCPG_tpm2[path_immu_gene[i],])
  PCPG_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  PCPG_exp_gene$sample<-as.character(PCPG_exp_gene$sample)
  
  x2<-t(PCPG_cnv[path_immu_gene[i],])
  PCPG_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  PCPG_cnv_gene$sample<-as.character(PCPG_cnv_gene$sample)
  
  PCPG_union<-inner_join(PCPG_exp_gene,PCPG_cnv_gene)
  PCPG_union$cancer<-rep("PCPG",dim(PCPG_union)[1])
  PCPG_union$CNV_group<-ifelse(PCPG_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  PCPG_l[[i]]<-PCPG_union
}
##########


SARC_tpm<-read.table("F:/workplace/SARC/SARC_readcount.genes.tpm.txt")

SARC_tpm<-(log2(SARC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(SARC_tpm)<-str_replace_all(colnames(SARC_tpm),"[.*]",replacement = "-")

SARC_gene<-intersect(rownames(SARC_tpm),path_immu_gene) 
SARC_tpm2<-SARC_tpm[SARC_gene,]  #免疫衰老基因的表达谱

SARC<-fread("f:/workplace/SARC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(SARC)<-SARC$`Gene Symbol`
SARC<-SARC[,-1]

both_gene<-intersect(path_immu_gene,rownames(SARC))

SARC_CNV<-SARC[both_gene,]
SARC_cnv<-ifelse(SARC_CNV<0,1,ifelse(SARC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
SARC_cnv<-as.data.frame(SARC_cnv)

SARC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(SARC_tpm2[path_immu_gene[i],])
  SARC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  SARC_exp_gene$sample<-as.character(SARC_exp_gene$sample)
  
  x2<-t(SARC_cnv[path_immu_gene[i],])
  SARC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  SARC_cnv_gene$sample<-as.character(SARC_cnv_gene$sample)
  
  SARC_union<-inner_join(SARC_exp_gene,SARC_cnv_gene)
  SARC_union$cancer<-rep("SARC",dim(SARC_union)[1])
  SARC_union$CNV_group<-ifelse(SARC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  SARC_l[[i]]<-SARC_union
}
###########3


LAML_tpm<-read.table("F:/workplace/LAML/LAML_readcount.genes.tpm.txt")

LAML_tpm<-(log2(LAML_tpm+1))  #表达谱变为log2(TPM+1)
colnames(LAML_tpm)<-str_replace_all(colnames(LAML_tpm),"[.*]",replacement = "-")

LAML_gene<-intersect(rownames(LAML_tpm),path_immu_gene) 
LAML_tpm2<-LAML_tpm[LAML_gene,]  #免疫衰老基因的表达谱

LAML<-fread("f:/workplace/LAML/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(LAML)<-LAML$`Gene Symbol`
LAML<-LAML[,-1]

both_gene<-intersect(path_immu_gene,rownames(LAML))

LAML_CNV<-LAML[both_gene,]
LAML_cnv<-ifelse(LAML_CNV<0,1,ifelse(LAML_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
LAML_cnv<-as.data.frame(LAML_cnv)

LAML_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(LAML_tpm2[path_immu_gene[i],])
  LAML_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  LAML_exp_gene$sample<-as.character(LAML_exp_gene$sample)
  
  x2<-t(LAML_cnv[path_immu_gene[i],])
  LAML_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  LAML_cnv_gene$sample<-as.character(LAML_cnv_gene$sample)
  
  LAML_union<-inner_join(LAML_exp_gene,LAML_cnv_gene)
  LAML_union$cancer<-rep("LAML",dim(LAML_union)[1])
  LAML_union$CNV_group<-ifelse(LAML_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  LAML_l[[i]]<-LAML_union
}
########


PAAD_tpm<-read.table("F:/workplace/PAAD/PAAD_readcount.genes.tpm.txt")

PAAD_tpm<-(log2(PAAD_tpm+1))  #表达谱变为log2(TPM+1)
colnames(PAAD_tpm)<-str_replace_all(colnames(PAAD_tpm),"[.*]",replacement = "-")

PAAD_gene<-intersect(rownames(PAAD_tpm),path_immu_gene) 
PAAD_tpm2<-PAAD_tpm[PAAD_gene,]  #免疫衰老基因的表达谱

PAAD<-fread("f:/workplace/PAAD/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(PAAD)<-PAAD$`Gene Symbol`
PAAD<-PAAD[,-1]

both_gene<-intersect(path_immu_gene,rownames(PAAD))

PAAD_CNV<-PAAD[both_gene,]
PAAD_cnv<-ifelse(PAAD_CNV<0,1,ifelse(PAAD_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
PAAD_cnv<-as.data.frame(PAAD_cnv)

PAAD_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(PAAD_tpm2[path_immu_gene[i],])
  PAAD_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  PAAD_exp_gene$sample<-as.character(PAAD_exp_gene$sample)
  
  x2<-t(PAAD_cnv[path_immu_gene[i],])
  PAAD_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  PAAD_cnv_gene$sample<-as.character(PAAD_cnv_gene$sample)
  
  PAAD_union<-inner_join(PAAD_exp_gene,PAAD_cnv_gene)
  PAAD_union$cancer<-rep("PAAD",dim(PAAD_union)[1])
  PAAD_union$CNV_group<-ifelse(PAAD_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  PAAD_l[[i]]<-PAAD_union
}

####3
ESCA_tpm<-read.table("F:/workplace/ESCA/ESCA_readcount.genes.tpm.txt")

ESCA_tpm<-(log2(ESCA_tpm+1))  #表达谱变为log2(TPM+1)
colnames(ESCA_tpm)<-str_replace_all(colnames(ESCA_tpm),"[.*]",replacement = "-")

ESCA_gene<-intersect(rownames(ESCA_tpm),path_immu_gene) 
ESCA_tpm2<-ESCA_tpm[ESCA_gene,]  #免疫衰老基因的表达谱

ESCA<-fread("f:/workplace/ESCA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(ESCA)<-ESCA$`Gene Symbol`
ESCA<-ESCA[,-1]

both_gene<-intersect(path_immu_gene,rownames(ESCA))

ESCA_CNV<-ESCA[both_gene,]
ESCA_cnv<-ifelse(ESCA_CNV<0,1,ifelse(ESCA_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
ESCA_cnv<-as.data.frame(ESCA_cnv)

ESCA_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(ESCA_tpm2[path_immu_gene[i],])
  ESCA_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  ESCA_exp_gene$sample<-as.character(ESCA_exp_gene$sample)
  
  x2<-t(ESCA_cnv[path_immu_gene[i],])
  ESCA_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  ESCA_cnv_gene$sample<-as.character(ESCA_cnv_gene$sample)
  
  ESCA_union<-inner_join(ESCA_exp_gene,ESCA_cnv_gene)
  ESCA_union$cancer<-rep("ESCA",dim(ESCA_union)[1])
  ESCA_union$CNV_group<-ifelse(ESCA_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  ESCA_l[[i]]<-ESCA_union
}
##########


TGCT_tpm<-read.table("F:/workplace/TGCT/TGCT_readcount.genes.tpm.txt")

TGCT_tpm<-(log2(TGCT_tpm+1))  #表达谱变为log2(TPM+1)
colnames(TGCT_tpm)<-str_replace_all(colnames(TGCT_tpm),"[.*]",replacement = "-")

TGCT_gene<-intersect(rownames(TGCT_tpm),path_immu_gene) 
TGCT_tpm2<-TGCT_tpm[TGCT_gene,]  #免疫衰老基因的表达谱

TGCT<-fread("f:/workplace/TGCT/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(TGCT)<-TGCT$`Gene Symbol`
TGCT<-TGCT[,-1]

both_gene<-intersect(path_immu_gene,rownames(TGCT))

TGCT_CNV<-TGCT[both_gene,]
TGCT_cnv<-ifelse(TGCT_CNV<0,1,ifelse(TGCT_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
TGCT_cnv<-as.data.frame(TGCT_cnv)

TGCT_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(TGCT_tpm2[path_immu_gene[i],])
  TGCT_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  TGCT_exp_gene$sample<-as.character(TGCT_exp_gene$sample)
  
  x2<-t(TGCT_cnv[path_immu_gene[i],])
  TGCT_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  TGCT_cnv_gene$sample<-as.character(TGCT_cnv_gene$sample)
  
  TGCT_union<-inner_join(TGCT_exp_gene,TGCT_cnv_gene)
  TGCT_union$cancer<-rep("TGCT",dim(TGCT_union)[1])
  TGCT_union$CNV_group<-ifelse(TGCT_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  TGCT_l[[i]]<-TGCT_union
}
########


THYM_tpm<-read.table("F:/workplace/THYM/THYM_readcount.genes.tpm.txt")

THYM_tpm<-(log2(THYM_tpm+1))  #表达谱变为log2(TPM+1)
colnames(THYM_tpm)<-str_replace_all(colnames(THYM_tpm),"[.*]",replacement = "-")

THYM_gene<-intersect(rownames(THYM_tpm),path_immu_gene) 
THYM_tpm2<-THYM_tpm[THYM_gene,]  #免疫衰老基因的表达谱

THYM<-fread("f:/workplace/THYM/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(THYM)<-THYM$`Gene Symbol`
THYM<-THYM[,-1]

both_gene<-intersect(path_immu_gene,rownames(THYM))

THYM_CNV<-THYM[both_gene,]
THYM_cnv<-ifelse(THYM_CNV<0,1,ifelse(THYM_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
THYM_cnv<-as.data.frame(THYM_cnv)

THYM_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(THYM_tpm2[path_immu_gene[i],])
  THYM_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  THYM_exp_gene$sample<-as.character(THYM_exp_gene$sample)
  
  x2<-t(THYM_cnv[path_immu_gene[i],])
  THYM_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  THYM_cnv_gene$sample<-as.character(THYM_cnv_gene$sample)
  
  THYM_union<-inner_join(THYM_exp_gene,THYM_cnv_gene)
  THYM_union$cancer<-rep("THYM",dim(THYM_union)[1])
  THYM_union$CNV_group<-ifelse(THYM_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  THYM_l[[i]]<-THYM_union
}
########


MESO_tpm<-read.table("F:/workplace/MESO/MESO_readcount.genes.tpm.txt")

MESO_tpm<-(log2(MESO_tpm+1))  #表达谱变为log2(TPM+1)
colnames(MESO_tpm)<-str_replace_all(colnames(MESO_tpm),"[.*]",replacement = "-")

MESO_gene<-intersect(rownames(MESO_tpm),path_immu_gene) 
MESO_tpm2<-MESO_tpm[MESO_gene,]  #免疫衰老基因的表达谱

MESO<-fread("f:/workplace/MESO/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(MESO)<-MESO$`Gene Symbol`
MESO<-MESO[,-1]

both_gene<-intersect(path_immu_gene,rownames(MESO))

MESO_CNV<-MESO[both_gene,]
MESO_cnv<-ifelse(MESO_CNV<0,1,ifelse(MESO_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
MESO_cnv<-as.data.frame(MESO_cnv)

MESO_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(MESO_tpm2[path_immu_gene[i],])
  MESO_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  MESO_exp_gene$sample<-as.character(MESO_exp_gene$sample)
  
  x2<-t(MESO_cnv[path_immu_gene[i],])
  MESO_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  MESO_cnv_gene$sample<-as.character(MESO_cnv_gene$sample)
  
  MESO_union<-inner_join(MESO_exp_gene,MESO_cnv_gene)
  MESO_union$cancer<-rep("MESO",dim(MESO_union)[1])
  MESO_union$CNV_group<-ifelse(MESO_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  MESO_l[[i]]<-MESO_union
}
########


UVM_tpm<-read.table("F:/workplace/UVM/UVM_readcount.genes.tpm.txt")

UVM_tpm<-(log2(UVM_tpm+1))  #表达谱变为log2(TPM+1)
colnames(UVM_tpm)<-str_replace_all(colnames(UVM_tpm),"[.*]",replacement = "-")

UVM_gene<-intersect(rownames(UVM_tpm),path_immu_gene) 
UVM_tpm2<-UVM_tpm[UVM_gene,]  #免疫衰老基因的表达谱

UVM<-fread("f:/workplace/UVM/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(UVM)<-UVM$`Gene Symbol`
UVM<-UVM[,-1]

both_gene<-intersect(path_immu_gene,rownames(UVM))

UVM_CNV<-UVM[both_gene,]
UVM_cnv<-ifelse(UVM_CNV<0,1,ifelse(UVM_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
UVM_cnv<-as.data.frame(UVM_cnv)

UVM_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(UVM_tpm2[path_immu_gene[i],])
  UVM_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  UVM_exp_gene$sample<-as.character(UVM_exp_gene$sample)
  
  x2<-t(UVM_cnv[path_immu_gene[i],])
  UVM_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  UVM_cnv_gene$sample<-as.character(UVM_cnv_gene$sample)
  
  UVM_union<-inner_join(UVM_exp_gene,UVM_cnv_gene)
  UVM_union$cancer<-rep("UVM",dim(UVM_union)[1])
  UVM_union$CNV_group<-ifelse(UVM_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  UVM_l[[i]]<-UVM_union
}
#########


DLBC_tpm<-read.table("F:/workplace/DLBC/DLBC_readcount.genes.tpm.txt")

DLBC_tpm<-(log2(DLBC_tpm+1))  #表达谱变为log2(TPM+1)
colnames(DLBC_tpm)<-str_replace_all(colnames(DLBC_tpm),"[.*]",replacement = "-")

DLBC_gene<-intersect(rownames(DLBC_tpm),path_immu_gene) 
DLBC_tpm2<-DLBC_tpm[DLBC_gene,]  #免疫衰老基因的表达谱

DLBC<-fread("f:/workplace/DLBC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(DLBC)<-DLBC$`Gene Symbol`
DLBC<-DLBC[,-1]

both_gene<-intersect(path_immu_gene,rownames(DLBC))

DLBC_CNV<-DLBC[both_gene,]
DLBC_cnv<-ifelse(DLBC_CNV<0,1,ifelse(DLBC_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
DLBC_cnv<-as.data.frame(DLBC_cnv)

DLBC_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(DLBC_tpm2[path_immu_gene[i],])
  DLBC_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  DLBC_exp_gene$sample<-as.character(DLBC_exp_gene$sample)
  
  x2<-t(DLBC_cnv[path_immu_gene[i],])
  DLBC_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  DLBC_cnv_gene$sample<-as.character(DLBC_cnv_gene$sample)
  
  DLBC_union<-inner_join(DLBC_exp_gene,DLBC_cnv_gene)
  DLBC_union$cancer<-rep("DLBC",dim(DLBC_union)[1])
  DLBC_union$CNV_group<-ifelse(DLBC_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  DLBC_l[[i]]<-DLBC_union
}
#######


CHOL_tpm<-read.table("F:/workplace/CHOL/CHOL_readcount.genes.tpm.txt")

CHOL_tpm<-(log2(CHOL_tpm+1))  #表达谱变为log2(TPM+1)
colnames(CHOL_tpm)<-str_replace_all(colnames(CHOL_tpm),"[.*]",replacement = "-")

CHOL_gene<-intersect(rownames(CHOL_tpm),path_immu_gene) 
CHOL_tpm2<-CHOL_tpm[CHOL_gene,]  #免疫衰老基因的表达谱

CHOL<-fread("f:/workplace/CHOL/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(CHOL)<-CHOL$`Gene Symbol`
CHOL<-CHOL[,-1]

both_gene<-intersect(path_immu_gene,rownames(CHOL))

CHOL_CNV<-CHOL[both_gene,]
CHOL_cnv<-ifelse(CHOL_CNV<0,1,ifelse(CHOL_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
CHOL_cnv<-as.data.frame(CHOL_cnv)

CHOL_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(CHOL_tpm2[path_immu_gene[i],])
  CHOL_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  CHOL_exp_gene$sample<-as.character(CHOL_exp_gene$sample)
  
  x2<-t(CHOL_cnv[path_immu_gene[i],])
  CHOL_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  CHOL_cnv_gene$sample<-as.character(CHOL_cnv_gene$sample)
  
  CHOL_union<-inner_join(CHOL_exp_gene,CHOL_cnv_gene)
  CHOL_union$cancer<-rep("CHOL",dim(CHOL_union)[1])
  CHOL_union$CNV_group<-ifelse(CHOL_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  CHOL_l[[i]]<-CHOL_union
}
#########


READ_tpm<-read.table("F:/workplace/READ/READ_readcount.genes.tpm.txt")

READ_tpm<-(log2(READ_tpm+1))  #表达谱变为log2(TPM+1)
colnames(READ_tpm)<-str_replace_all(colnames(READ_tpm),"[.*]",replacement = "-")

READ_gene<-intersect(rownames(READ_tpm),path_immu_gene) 
READ_tpm2<-READ_tpm[READ_gene,]  #免疫衰老基因的表达谱

READ<-fread("f:/workplace/READ/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(READ)<-READ$`Gene Symbol`
READ<-READ[,-1]

both_gene<-intersect(path_immu_gene,rownames(READ))

READ_CNV<-READ[both_gene,]
READ_cnv<-ifelse(READ_CNV<0,1,ifelse(READ_CNV==2,1,0))#拷贝数缺失和高水平拷贝数扩增的样本
READ_cnv<-as.data.frame(READ_cnv)

READ_l<-list()

for(i in 1:length(path_immu_gene)){
  
  x<-t(READ_tpm2[path_immu_gene[i],])
  READ_exp_gene<-data.frame(sample=rownames(x),mRNA=x)
  READ_exp_gene$sample<-as.character(READ_exp_gene$sample)
  
  x2<-t(READ_cnv[path_immu_gene[i],])
  READ_cnv_gene<-data.frame(sample=rownames(x2),cnv=x2[,1])
  READ_cnv_gene$sample<-as.character(READ_cnv_gene$sample)
  
  READ_union<-inner_join(READ_exp_gene,READ_cnv_gene)
  READ_union$cancer<-rep("READ",dim(READ_union)[1])
  READ_union$CNV_group<-ifelse(READ_union$cnv==1,"CNV","Non-CNV")%>%as.factor()
  READ_l[[i]]<-READ_union
}
#####

# 47个基因boxplot的输入文件 -------------------------------------------------------

box_l<-list()  #47个基因在泛癌中的mRNA表达
for(i in 1:length(path_immu_gene)){
  box_l[[i]]<-rbind(LUAD_l[[i]],BLCA_l[[i]],HNSC_l[[i]],KIRC_l[[i]],KIRP_l[[i]],LIHC_l[[i]],LUSC_l[[i]],THCA_l[[i]],COAD_l[[i]],GBM_l[[i]],KICH_l[[i]],LGG_l[[i]],BRCA_l[[i]],READ_l[[i]],UCS_l[[i]],UCEC_l[[i]],OV_l[[i]],PRAD_l[[i]],STAD_l[[i]],SKCM_l[[i]],CESC_l[[i]],ACC_l[[i]],PCPG_l[[i]],SARC_l[[i]],LAML_l[[i]],PAAD_l[[i]],ESCA_l[[i]],TGCT_l[[i]],THYM_l[[i]],MESO_l[[i]],UVM_l[[i]],DLBC_l[[i]],CHOL_l[[i]])
  box_l[[i]]$gene=rep(colnames(box_l[[i]])[2],dim(box_l[[i]])[1])
  colnames(box_l[[i]])[2]<-"mRNA_expression"
}

save(box_l,file="box_input_CNVgroup.RData")
load("box_input_CNVgroup.RData")

mRNA<-c()
for(i in 1:length(path_immu_gene)){
  if(!box_l[[i]]$cnv%in%NA){
  mRNA_diff<-compare_means(mRNA_expression ~ CNV_group, data = box_l[[i]],method = "t.test")
  mRNA[i]<-mRNA_diff$p.adj
  }
}
names(mRNA)<-path_immu_gene
which(mRNA<0.05)%>%length
#有36个基因的表达值差异显著，样本按是否CNV分组

cnv_effect_mRNA<-data.frame()
for(i in 1:length(path_immu_gene)){
 cnv_compare<-box_l[[i]] %>%
  group_by(CNV_group) %>%
  dplyr::summarise(mean_mRNA=mean(mRNA_expression))

 effect<-ifelse(cnv_compare$mean_mRNA[1]>cnv_compare$mean_mRNA[2],"activate","inactivate")
 df<-data.frame(gene=box_l[[i]]$gene[1],cnv_effect=effect)
 cnv_effect_mRNA<-rbind(cnv_effect_mRNA,df)
}

mRNA_result<-data.frame(gene=names(mRNA),diff_pvalue=mRNA)
mRNA_cnv_result<-inner_join(mRNA_result,cnv_effect_mRNA)
#45个基因按CNV，non-CNV分组后，mRNA的差异，以及CNV对MRNA的激活或失活统计

write.csv(mRNA_cnv_result,"f:/workplace/结果2/mRNA-boxplot/CNV分组的38个显著基因/45基因激活-失活统计.csv",row.names=F)
# 差异显著的38个基因-mRNA_boxplot ----------------------------------------------------

df38<-data.frame()
for(i in which(mRNA<0.05)%>%as.numeric()){
  
  df38<-rbind(df38,box_l[[i]])
}

gene38<-df38$gene%>%unique()

setwd("f:/workplace/结果2/mRNA-boxplot/CNV分组的38个显著基因/")
for(i in 1:length(gene38)){
  
  df<-subset(df38,df38$gene==gene38[i])
  pic<-ggboxplot(df, x = "CNV_group", y = "mRNA_expression",add.params = list(color = "cancer",size=0.1),add = "jitter",
                 palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                    "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                    "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                    "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(method = "t.test")+ylab(gene38[i])+xlab("")+theme(legend.position = "none")
  ggsave(pic, file=paste(gene38[i],".pdf",sep=""), width=8, height=8)
}


# 激活基因的box图 ---------------------------------------------------------------
setwd("f:/workplace/结果2/mRNA-boxplot/激活6基因box/")
df6<-data.frame()
gene6<-subset(mRNA_cnv_result,mRNA_cnv_result$cnv_effect=="activate")

for(i in rownames(gene6)%>%as.numeric()){
  
  df6<-rbind(df6,box_l[[i]])

}
gene6_name<-df6$gene%>%unique()


plot_list<-list()
for(i in 1:length(gene6_name) ){
  df<-subset(df6,df6$gene==gene6_name[i])
  plot_list[[i]]<-ggboxplot(df, x = "CNV_group", y = "mRNA_expression",outlier.size = 0,add.params = list(color = "cancer",alpha=0.2,size=0.01),add = "jitter",
                palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                   "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                   "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                   "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(label = "p.signif",method = "t.test")+ylab(paste(gene6_name[i]," expression"))+xlab("")+theme(legend.position = "none")
  
}
plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]])

# 36gene ------------------------------------------------------------------

plot_list<-list()
for(i in 1:length(gene38) ){
  df<-subset(df38,df38$gene==gene38[i])
  plot_list[[i]]<-ggboxplot(df, x = "CNV_group", y = "mRNA_expression", fill="CNV_group", palette = c("#00AFBB", "#E7B800"))+stat_compare_means(label = "p.signif",label.x.npc="center",label.y.npc="center",method = "t.test")+ylab(paste(gene38[i]," expression"))+xlab("")+theme(axis.title.y = element_text(size = 6.5),axis.text.x = element_text(size = 6,color="black"),axis.text.y = element_text(size = 6,color="black") ,legend.position = "none")
  
  
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
          plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]],
          plot_list[[19]],plot_list[[20]],plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],plot_list[[26]],
          plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],plot_list[[31]],plot_list[[32]],
          plot_list[[33]],plot_list[[34]],plot_list[[35]],plot_list[[36]])
# 失活的基因box图 ---------------------------------------------------------------

df_inac<-data.frame()
gene_inac<-subset(mRNA_cnv_result,mRNA_cnv_result$cnv_effect=="inactivate")

for(i in rownames(gene_inac)%>%as.numeric()){
  
  df_inac<-rbind(df_inac,box_l[[i]])
  
}
gene_name<-df_inac$gene%>%unique()
gene_name<-gene_name[-32]
gene_name<-gene_name[-32]
symnum <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_list<-list()
for(i in 1:length(gene_name) ){
  df<-subset(df_inac,df_inac$gene==gene_name[i])
  plot_list[[i]]<-ggboxplot(df, x = "CNV_group", y = "mRNA_expression",outlier.size = 0,size=0.01,add.params = list(color = "cancer",alpha=0.2,size=0.001),add = "jitter",
                            palette = mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
                                               "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
                                               "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
                                               "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3"))+stat_compare_means(label = "p.signif", method = "t.test",size=3.8)+ylab(paste(gene_name[i]," expression"))+xlab("")+theme(axis.title.y = element_text(size = 6.5),axis.text.x = element_text(size = 6,color="black"),axis.text.y = element_text(size = 6,color="black") ,legend.position = "none")
  
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],plot_list[[6]],
          plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
          plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],plot_list[[17]],plot_list[[18]],
          plot_list[[19]],plot_list[[20]],plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],plot_list[[26]],
          plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],plot_list[[31]],plot_list[[32]],
          plot_list[[33]],plot_list[[34]],plot_list[[35]],plot_list[[36]])

plot_grid(plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],plot_list[[26]],
          plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],plot_list[[31]],plot_list[[32]],
          plot_list[[33]],plot_list[[34]],plot_list[[35]],plot_list[[36]])
