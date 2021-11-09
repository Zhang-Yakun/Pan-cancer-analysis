


# Activated gene and Repressed genes --------------------------------------------------------


activate_both_gene<-c("FAS","TP53","ATM","BCL2","IGF1","TNF","FASLG","NGF","BCL2L11","IFNG","JAK2","MTOR","TLR2")


supp_both_gene<-c("FASLG","IGF1R","TP53","NGF","IGF1","JAK2","JAK3","IL6","IL7","IL7R","IL2","IL4","BCL2","MTOR","PRKAA2","STAT3")


pathway_both_gene<-union(activate_both_gene,supp_both_gene)



# Multi-omics analysis of immunosenescence genes ----------------------------------------------------------

mut_freq<-read.table("f:/workplace/结果4/mut_freq.txt")
mut_freq2<-t(mut_freq)%>%as.data.frame()
mut_freq2$V1<-str_remove(mut_freq2$V1,"%")%>%as.numeric()
min(mut_freq2$V1)
max(mut_freq2$V1)
quantile(mut_freq2$V1)

mut_top10<-arrange(mut_freq2,desc(V1))[1:10,]


cnv_freq<-read.table("f:/workplace/结果4/cnv_freq.txt")
cnv_freq2<-t(cnv_freq)%>%as.data.frame()
cnv_freq2$V1<-str_remove(cnv_freq2$V1,"%")%>%as.numeric()

cnv_top10<-arrange(cnv_freq2,desc(V1))[1:10,]


mRNA_mutgroup<-read.csv("F:/workplace/结果2/mRNA-boxplot/突变分组的显著基因-boxplot/免疫衰老基因激活-失活统计.csv")
mRNA_mutgroup_gene<-subset(mRNA_mutgroup,mRNA_mutgroup$diff_pvalue<0.05)


mRNA_CNVgroup<-read.csv("f:/workplace/结果2/mRNA-boxplot/CNV分组的38个显著基因/45基因激活-失活统计.csv")
mRNA_CNVgroup_gene<-subset(mRNA_CNVgroup,mRNA_CNVgroup$diff_pvalue<0.05)


protein_gene<-c("ATM","BCL2","BCL2L11","MTOR","TP53")


# Sex-biased immunosenescence genes -----------------------------------------------------------



LUAD_exp<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/LUAD/表达/LUAD_biax_immu_gene_exp.csv")
BLCA_exp<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/BLCA/表达/BLCA_biax_immu_gene_exp.csv")
HNSC_exp<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/HNSC/表达/HNSC_biax_immu_gene_exp.csv")
KIRC_exp<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/KIRC/表达/KIRC_biax_immu_gene_exp.csv")
KIRP_exp<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/KIRP/表达/KIRP_biax_immu_gene_exp.csv")
LIHC_exp<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/LIHC/表达/LIHC_biax_immu_gene_exp.csv")
THCA_exp<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/THCA/表达/THCA_biax_immu_gene_exp.csv")
LUSC_exp<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/LUSC/表达/LUSC_biax_immu_gene_exp.csv")

sex_diff_exp<-rbind(LUAD_exp,BLCA_exp,HNSC_exp,KIRC_exp,KIRP_exp,LIHC_exp,THCA_exp,LUSC_exp)

write.csv(sex_diff_exp,"f:/workplace/结果4/sex_diff_exp.csv",row.names = F)



BLCA_methy<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/BLCA/甲基化/BLCA_甲基化差异免疫衰老基因.csv")
HNSC_methy<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/HNSC/甲基化/HNSC_甲基化差异免疫衰老基因.csv")
KIRC_methy<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/KIRC/甲基化/KIRC_甲基化差异免疫衰老基因.csv")
KIRP_methy<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/KIRP/甲基化/KIRP_甲基化差异免疫衰老基因.csv")
LIHC_methy<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/LIHC/甲基化/LIHC_甲基化差异免疫衰老基因.csv")
LUSC_methy<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/LUSC/甲基化/LUSC_甲基化差异免疫衰老基因.csv")
THCA_methy<-read.csv("C:/Users/Jennifer/Desktop/免疫衰老课题/结果3-性别差异图/THCA/甲基化/THCA_甲基化差异免疫衰老基因.csv")

sex_diff_methy<-rbind(BLCA_methy,HNSC_methy,KIRC_methy,KIRP_methy,LIHC_methy,THCA_methy,LUSC_methy)

write.csv(sex_diff_methy,"f:/workplace/结果4/sex_diff_methy.csv",row.names = F)

sex_diff_exp$X%>%unique()%>%length()
#41 gene
sex_exp_gene<-sex_diff_exp$X%>%unique()

sex_diff_methy$gene%>%unique()%>%length()
#28 gene
sex_methy_gene<-sex_diff_methy$gene%>%unique()



# ISC-pathway prognosis analysis -----------------------------------------------------------------



library(survival)
library(survminer)
library(data.table)
library(dplyr)


surdata<-fread("Survival_SupplementalTable_S1_20171025_xena_sp")%>%as.data.frame()
unique(surdata$`cancer type abbreviation`)

library(survival)
library(survminer)
library(msigdbr)
library(tibble)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)

options(stringsAsFactors = FALSE) 


mysubpathway<-list()

mysubpathway[[1]]<-c("ATM","TP53","SESN1","IGF1","BCL2","FAS")
names(mysubpathway)[1]<-"p53_subpathway"

mysubpathway[[2]]<-c("TNF","HLA-A","HLA-B","KIR3DL1")
names(mysubpathway)[2]<-"antigen_subpathway"

mysubpathway[[3]]<-c("ATM","TP53","FASLG","FAS","TNF","BCL2L11","BCL2")
names(mysubpathway)[3]<-"apoptosis_subpathway"

mysubpathway[[4]]<-c("CX3CL1","CX3CR1","CXCL9","CCR5","CCR7","JAK2","JAK3","STAT3")
names(mysubpathway)[4]<-"chemokine_subpathway"

mysubpathway[[5]]<-c("IL2","IL4","IL6","IL7","JAK2","JAK3","MTOR","STAT3","STAT5A","BCL2")
names(mysubpathway)[5]<-"JAK_STAT_subpathway"

mysubpathway[[6]]<-c("TP53","PRKAA2","KL","MTOR","IGF1R","FOXO1")
names(mysubpathway)[6]<-"longevity_subpathway"

mysubpathway[[7]]<-c("FASLG","FAS","IGF1","IGF1R","MAPK11","TP53","TNF","TAB1")
names(mysubpathway)[7]<-"MAPK_subpathway"

mysubpathway[[8]]<-c("TLR2","TLR9","JAK2","STAT3","CD28","PDCD1","CD247","MTOR","MAPK11")
names(mysubpathway)[8]<-"PDL1_subpathway"

mysubpathway[[9]]<-c("TLR2","IGF1R","IL2","IL4","IL7","TP53","BCL2","BCL2L11","FASLG")
names(mysubpathway)[9]<-"PI3K_subpathway"


# 读入基因表达数据 ----------------------------------------------------------------
BLCA_tpm<-fread("f:/workplace/BLCA/BLCA_readcount.genes.tpm.txt")%>%as.data.frame()
rownames(BLCA_tpm)<-BLCA_tpm$V1
BLCA_tpm<-BLCA_tpm[,-1]

expMatrix_tpm<-(log2(BLCA_tpm+1))  

# GSVA分析 ------------------------------------------------------------------
gsva_BLCA <- gsva(as.matrix(expMatrix_tpm), mysubpathway)
head(gsva_BLCA)  #GSVA result

# survival analysis --------------------------------------------------------------------

BLCA_sample <- intersect(colnames(gsva_BLCA),surdata$sample) 

gsva_BLCA2 <- gsva_BLCA[,BLCA_sample]
survival_BLCA <- subset(surdata,surdata$sample%in%BLCA_sample)

setwd("f:/workplace/结果4/显著的生存图/")
library(cowplot)

BLCA_p<-c()
for(i in 1:dim(gsva_BLCA2)[1]){
  
  Group<-ifelse(gsva_BLCA2[i,]>median(gsva_BLCA2[i,]),"high","low")
  g<-data.frame(sample=names(Group),group=Group)
  sur<-inner_join(g,survival_BLCA)
  
  fit <- survfit(Surv(OS.time, OS) ~ group, data = sur)
  diff <- survdiff(Surv(OS.time, OS) ~ group, data = sur)
  
  BLCA_p[i]<-pchisq(diff$chisq, length(diff$n) - 1, lower.tail = FALSE)
  names(BLCA_p)[i]<-(rownames(gsva_BLCA2)[i])
  
  if(BLCA_p[i]<0.05){
    subpathway <- rownames(gsva_BLCA2)[i]
    gg<-ggsurvplot(fit, pval = TRUE,linetype = c("solid", "dashed"), 
                   palette = c("blue","red"),
                   title = "BLCA",
                   legend.title=subpathway,legend=c(0.9,0.9),
                   conf.int = F)
    
    gg2 <- plot_grid(gg$plot + theme(legend.text = element_text(size = 12)), gg$table, ncol = 1, rel_widths = 1, rel_heights = c(1, .4), align = "hv")
    pdf(paste(paste("BLCA_sur",i,sep=""),".pdf",sep = ""))
    print(gg2)
    dev.off()
  }
}



gg <- ggsurvplot(sfit,
                 risk.table = TRUE,
                 conf.int = F, # 置信区间
                 palette = c("#E63946", "#3A86FF"),
                 risk.table.y.text = F,
                 pval.method = F,
                 pval = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", pval %>% round(digits = 3))),
                 ncensor.plot = F,
                 surv.median.line = "hv",
                 linetype = "strata", xlab = "Time(Months)",
                 legend = c(0.8, 0.95),
                 legend.title = "",
                 legend.labs = c("high-Expr", "low-Expr")
)
library(cowplot)
gg2 <- plot_grid(gg$plot + theme(legend.text = element_text(size = 12)), gg$table, ncol = 1, rel_widths = 1, rel_heights = c(1, .4), align = "hv")
pdf("survplot_TBX5.pdf", width = 4.5, height = 5.5)
print(gg2)
dev.off()
