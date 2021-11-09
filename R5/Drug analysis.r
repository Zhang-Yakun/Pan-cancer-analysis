
# Drug-target analysis --------------------------------------------------------------------

setwd("f:/workplace/结果5/")
sex_df_cancer <-fread("f:/workplace/结果4/sex_gene_cancer.csv")%>%as.data.frame()
sex_HNSC<-subset(sex_df_cancer,sex_df_cancer$cancer=="HNSC")

HNSC_sex_gene<-sex_HNSC$gene%>%unique()


library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
bitr(geneID=HNSC_sex_gene, fromType="SYMBOL", toType="GENENAME", OrgDb="org.Hs.eg.db", drop = TRUE)


# gene expression---------------------------------------------------------------

HNSC_DEG<-fread("f:/workplace/HNSC/HNSC_limma_exp.csv")%>%as.data.frame()

HNSC_diff<-subset(HNSC_DEG,HNSC_DEG$adj.P.Val<0.05)


gene90<-fread("f:/workplace/GeneSymbol_90.txt",header=F)
gene90<-gene90$V1%>%unique()

HNSC_immu_diff <-subset(HNSC_diff,HNSC_diff$V1%in%gene90)


# CMAP analysis, up_list file and down_list file ---------------------------------------------

HNSC_up<-subset(HNSC_immu_diff,HNSC_immu_diff$logFC>0)
HNSC_up_gene<-HNSC_up$V1

HNSC_up_probe<-bitr(geneID=HNSC_up_gene, fromType="SYMBOL", toType="UNIGENE", OrgDb="org.Hs.eg.db", drop = TRUE)
HNSC_up_probe2<-HNSC_up_probe$UNIGENE

HNSC_down<-subset(HNSC_immu_diff,HNSC_immu_diff$logFC<0)
HNSC_down_gene<-HNSC_down$V1

HNSC_down_probe<-bitr(geneID=HNSC_down_gene, fromType="SYMBOL", toType="UNIGENE", OrgDb="org.Hs.eg.db", drop = TRUE)
HNSC_down_probe2<-HNSC_down_probe$UNIGENE

write.table(HNSC_up_probe2,"f:/workplace/结果5/HNSC_up_probe.txt",quote=F,row.names=F,col.names=F)

write.table(HNSC_down_probe2,"f:/workplace/结果5/HNSC_down_probe.txt",quote=F,row.names=F,col.names=F)



# detailed_result ----------------------------------------------

HNSC_cmap_result<-fread("f:/workplace/结果5/HNSC/detailedResults-HNSC.csv")%>%as.data.frame()

HNSC_sensitive_drug<-subset(HNSC_cmap_result,HNSC_cmap_result$score > 0)

HNSC_sen_drug<-HNSC_sensitive_drug$`cmap name`


# Drugbank data------------------------------------------------------

drugbank<-fread("f:/workplace/结果5/HNSC/uniprot_Drugbank.csv")%>%as.data.frame()

library(stringr)
HNSC_sen_drug<-stringr::str_replace(HNSC_sen_drug,"\\b[a-z]",toupper)


drug_target<-subset(drugbank,drugbank$Name%in%HNSC_sen_drug)


target_ID<-bitr(geneID=drug_target$`UniProt ID`, fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)


target_sex_ID<-subset(target_ID,target_ID$SYMBOL%in%HNSC_sex_gene)

colnames(target_sex_ID)[1]<-colnames(drug_target)[4]

drug_sex_target<-subset(drug_target,drug_target$`UniProt ID`%in%target_sex_ID$UNIPROT)
drug_sex_HNSC<-inner_join(target_sex_ID,drug_sex_target)


sex_HNSC2<-subset(sex_HNSC,sex_HNSC$gene%in%drug_sex_HNSC$SYMBOL)
colnames(sex_HNSC2)[1]<-colnames(drug_sex_HNSC)[2]
HNSC_sex_drug<-inner_join(drug_sex_HNSC,sex_HNSC2)


write.csv(HNSC_sex_drug,"f:/workplace/结果5/HNSC/HNSC性别差异药物结果.csv",row.names=F)


# CCLE IC50 value data ------------------------------------------------

ccle<-data.table::fread("f:/workplace/结果5/CCLE.csv")%>%as.data.frame()


ccle_sex<-fread("f:/workplace/结果5/ccle_clinical_patient.txt")%>%as.data.frame()


gender<-unique(ccle_sex$SEX)[1:2]
ccle_sex2<-subset(ccle_sex,ccle_sex$SEX%in%gender)


ccle_sex3<-ccle_sex2[,c("PATIENT_ID","HISTOLOGY","SEX")]
colnames(ccle_sex3)[1]<-colnames(ccle)[1]

cc<-inner_join(ccle_sex3,ccle)



# HNSC cell lines ----------------------------------------------------

cc2<-subset(cc,cc$HISTOLOGY=="Carcinoma") 

cell_line<-cc2$`CCLE Cell Line Name`%>%str_split("_",n=2)


HNSC_cell<-c()

for(i in 1:length(cell_line)){
  
  if(cell_line[[i]][2]=="UPPER_AERODIGESTIVE_TRACT"){
    HNSC_cell<-c(HNSC_cell,i)
    
  }
}

HNSC_cell_line<-cc2[HNSC_cell,]



# sex difference analysis of IC50  -------------------------------------------------------

candicate_drug<-HNSC_cell_line$Compound%>%unique()

df=data.frame()

for(i in 1:length(candicate_drug)){
  drug1<-subset(HNSC_cell_line,HNSC_cell_line$Compound==candicate_drug[i])
  
  drug1_male_IC50<-subset(drug1,drug1$SEX=="Male")[,"IC50 (uM)"]
  drug1_female_IC50<-subset(drug1,drug1$SEX=="Female")[,"IC50 (uM)"]
  
  if( length(drug1_male_IC50)!=0 && length(drug1_female_IC50) !=0){
    p<-wilcox.test(drug1_male_IC50,drug1_female_IC50)

    m<-ifelse(median(drug1_male_IC50)<median(drug1_female_IC50),"male_sensity","female_sensity")
    
    output<-data.frame(drug=candicate_drug[i],wilcox.p=p$p.value,sex_bias=m)
    df<-rbind(df,output)
  }   
}


# GDSC IC50 analysis ------------------------------------------------------



cell_line_sex<-read.csv("model_list_20210324.csv",stringsAsFactors =F)
cell_sex<-cell_line_sex[,c("model_name","COSMIC_ID","gender")]
cell_sex<-subset(cell_sex,cell_sex$gender != "Unknown")
cell_sex$gender<-cell_sex$gender%>%as.character()

#GDSC-dataset1
cell_line_IC50<-fread("GDSC1_fitted_dose_response_25Feb20.csv")%>%as.data.frame()
cell_IC50<-cell_line_IC50[,c("DATASET","CELL_LINE_NAME","COSMIC_ID","TCGA_DESC","DRUG_NAME","PUTATIVE_TARGET","PATHWAY_NAME","LN_IC50")]

#GDSC-dataset2
cell_line_IC50_2<-fread("GDSC2_fitted_dose_response_25Feb20.csv")%>%as.data.frame()
cell_IC50_2<-cell_line_IC50_2[,c("DATASET","CELL_LINE_NAME","COSMIC_ID","TCGA_DESC","DRUG_NAME","PUTATIVE_TARGET","PATHWAY_NAME","LN_IC50")]

#merge dataset1 and dataset2
cell_IC50_all<-rbind(cell_IC50,cell_IC50_2)


cell_IC50_HNSC<-subset(cell_IC50_all,cell_IC50_all$TCGA_DESC=="HNSC")


cell_HNSC_sex<-intersect(cell_IC50_HNSC$CELL_LINE_NAME,cell_sex$model_name)

cell_IC50_HNSC_sex<-subset(cell_IC50_HNSC,cell_IC50_HNSC$CELL_LINE_NAME%in%cell_HNSC_sex)
cell_IC50_HNSC_sex$COSMIC_ID<-cell_IC50_HNSC_sex$COSMIC_ID%>%as.character()
cell_IC50_HNSC_sex<-inner_join(cell_IC50_HNSC_sex,cell_sex)


cell_IC50_HNSC_sex<-cell_IC50_HNSC_sex[,-1]%>%distinct()



candicate_drug<-cell_IC50_HNSC_sex$DRUG_NAME%>%unique()

HNSC_df=data.frame()

for(i in 1:length(candicate_drug)){
  drug1<-subset(cell_IC50_HNSC_sex,cell_IC50_HNSC_sex$DRUG_NAME==candicate_drug[i])
  
  drug1_male_IC50<-subset(drug1,drug1$gender=="Male")[,"LN_IC50"]
  drug1_female_IC50<-subset(drug1,drug1$gender=="Female")[,"LN_IC50"]
  
  if( length(drug1_male_IC50)!=0 && length(drug1_female_IC50) !=0){
    p<-wilcox.test(drug1_male_IC50,drug1_female_IC50)

    m<-ifelse(median(drug1_male_IC50)<median(drug1_female_IC50),"male_sensity","female_sensity")
    
    output<-data.frame(drug=candicate_drug[i],wilcox.p=p$p.value,sex_bias=m)
    HNSC_df<-rbind(HNSC_df,output)
  }   
}

HNSC_drug<-subset(HNSC_df,HNSC_df$wilcox.p<0.05)


write.csv(HNSC_drug,"HNSC_sexdiff_drug.csv",row.names = F)
HNSC_result<-cell_IC50_HNSC_sex%>%subset(cell_IC50_HNSC_sex$DRUG_NAME%in%HNSC_drug$drug)



# GDSC-target drug -----------------------------------------------------


HNSC_sex_gene<-read.table("f:/workplace/结果5/HNSC/HNSC性别差异免疫衰老基因.txt")
HNSC_sex_gene<-HNSC_sex_gene$V1%>%unique()%>%as.character()


library(stringr)
row_ID<-c()
for(i in 1:dim(HNSC_result)[1]){
  
  target<-str_split(HNSC_result$PUTATIVE_TARGET[i],",")%>%unlist()
  both_gene<-intersect(target,HNSC_sex_gene)
  if(length(both_gene)!=0){
    row_ID<-c(row_ID,i)
  }
}

#drug_result<-HNSC_result[row_ID,]
library(clusterProfiler)
library(org.Hs.eg.db)
#bitr(geneID="BCL-XL", fromType="SYMBOL", toType="ALIAS", OrgDb="org.Hs.eg.db", drop = TRUE)

drug_result<-subset(HNSC_result,HNSC_result$PUTATIVE_TARGET=="MTORC1, MTORC2")


# CMAP ------------------------------------------------------


HNSC_cmap_result<-fread("f:/workplace/结果5/cmap/detailedResults_HNSC_sex.csv")%>%as.data.frame()

HNSC_female_drug<-subset(HNSC_cmap_result,HNSC_cmap_result$score > 0)
HNSC_male_drug<-subset(HNSC_cmap_result,HNSC_cmap_result$score < 0)

HNSC_female_drug<-HNSC_female_drug$`cmap name`

HNSC_male_drug<-HNSC_male_drug$`cmap name`


# Drugbank-----------------------------------------------------

drugbank<-fread("f:/workplace/结果5/HNSC/uniprot_Drugbank.csv")%>%as.data.frame()

library(stringr)
HNSC_male_drug<-str_replace(HNSC_male_drug,"\\b[a-z]",toupper)


drug_target<-subset(drugbank,drugbank$Name%in%HNSC_male_drug)


target_ID<-bitr(geneID=drug_target$`UniProt ID`, fromType="UNIPROT", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)


target_male_ID<-subset(target_ID,target_ID$SYMBOL%in%HNSC_sex_gene)

colnames(target_male_ID)[1]<-colnames(drug_target)[4]

drug_male_target<-subset(drug_target,drug_target$`UniProt ID`%in%target_male_ID$`UniProt ID`)
drug_male_HNSC<-inner_join(target_male_ID,drug_male_target)


male_HNSC<-subset(sex_HNSC,sex_HNSC$gene%in%drug_male_HNSC$SYMBOL)
colnames(male_HNSC)[1]<-colnames(drug_male_HNSC)[2]
HNSC_male_drug<-inner_join(drug_male_HNSC,sex_HNSC2)



# CTRP---------------------------------------------------------------



auc<-fread("CTRP_auc.txt")%>%as.data.frame()

auc2<-auc[,c("experiment_id","area_under_curve","master_cpd_id")]


ccle_cellline<-fread("v20.meta.per_cell_line.txt")%>%as.data.frame()

cell<-ccle_cellline[,c("master_ccl_id","ccl_name","ccle_primary_site")]

colnames(cell)[1]<-"experiment_id"

cell_info<-fread("CCLE_clinical_patient.txt")%>%as.data.frame()

a<-cell_info$PATIENT_ID%>%str_split("_",n=2,simplify = TRUE)
cell_info$ccl_name<-a[,1]
cell_sex<-cell_info[,c("ccl_name","SEX")]


cell2<-inner_join(cell,cell_sex)


mydata<-inner_join(cell2,auc2)
HNSC_data<-subset(mydata,mydata$ccle_primary_site=="upper_aerodigestive_tract")


compound<-fread("v20.meta.per_compound.txt")%>%as.data.frame()

com<-compound[,1:2]
HNSC_drug<-inner_join(HNSC_data,com)

drug_name<-HNSC_drug$cpd_name%>%unique()

CTRP_HNSC_df=data.frame()

for(i in 1:length(drug_name)){
  drug1<-subset(HNSC_drug,HNSC_drug$cpd_name==drug_name[i])
  
  drug1_male_AUC<-subset(drug1,drug1$SEX=="Male")[,"area_under_curve"]
  drug1_female_AUC<-subset(drug1,drug1$SEX=="Female")[,"area_under_curve"]
  
  if( length(drug1_male_AUC)!=0 && length(drug1_female_AUC) !=0){
    p<-wilcox.test(drug1_male_AUC,drug1_female_AUC)
   
    m<-ifelse(median(drug1_male_AUC)>median(drug1_female_AUC),"male_sensity","female_sensity")
    
    output<-data.frame(drug=drug_name[i],wilcox.p=p$p.value,sex_bias=m)
    CTRP_HNSC_df<-rbind(CTRP_HNSC_df,output)
  }   
}
CTRP_HNSC_drug<-subset(CTRP_HNSC_df,CTRP_HNSC_df$wilcox.p<0.05)


write.csv(CTRP_HNSC_drug,"CTRP_HNSC_sex_drug.csv",row.names = F)


# GDSC—IC50 boxplot ---------------------------------------------------

OSI_IC50<-subset(cell_IC50_HNSC_sex,cell_IC50_HNSC_sex$DRUG_NAME=="OSI-027")


OSI_IC50$gender<-as.factor(OSI_IC50$gender)

library(ggplot2)
library(ggpubr)
library("RColorBrewer")
display.brewer.all(type = "qual")

ggboxplot(OSI_IC50, x = "gender", y = "LN_IC50",outlier.size = 0,add.params = list(color = "gender"),add = "jitter",
          palette = mycol<-c("#FB8072","#1F78B4"))+
  stat_compare_means(method = "wilcox.test")+
  ylab("IC50 of OSI-027")+xlab("")+
  theme(axis.title.y = element_text(size = 14),axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black") )


# crosstalk-drug-target ------------------------------------------------------
library(devtools)
devtools::install_github("zzwch/crosslinks")
library(crosslinks)
library(ggplot2)

mynode<-fread("f:/workplace/结果5/mynode.csv")%>%as.data.frame()

columns <- list(
  Cancer = mynode$id[mynode$type == "cancer"],
  Drug = mynode$id[mynode$type == "drug"],
  Target = mynode$id[mynode$type == "target"], 
  Pathway = mynode$id[mynode$type == "pathway"]
)

myedge<-fread("f:/workplace/结果5/myedge.csv")%>%as.data.frame()


columnCross2(myedge, mynode, columns,
             height = 1, flank_mult = rep(0.1, length(columns)), segment_shrink = 0.1,
             
             linetype = 1, # 默认值，所有连线都画成直线
             line_color = "color",
             
             pt_shape = 21, # 统一画成实心圆
             pt_size = "size",
             pt_alpha = .8, # 统一设置透明度
             pt_color = "color", pt_fill = "color") # 让空心圆的边粗一些
