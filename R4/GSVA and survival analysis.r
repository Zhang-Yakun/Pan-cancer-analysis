

library(data.table)
library(dplyr)


# pan-cancer survival data ----------------------------------------------------------------

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

# 9 ISC-pathways ---------------------------------------------------------------

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


# gene expression data ----------------------------------------------------------------
BLCA_tpm<-fread("f:/workplace/BLCA/BLCA_readcount.genes.tpm.txt")%>%as.data.frame()
rownames(BLCA_tpm)<-BLCA_tpm$V1
BLCA_tpm<-BLCA_tpm[,-1]

expMatrix_tpm<-(log2(BLCA_tpm+1)) 

# GSVA analysis ------------------------------------------------------------------
gsva_BLCA <- gsva(as.matrix(expMatrix_tpm), mysubpathway)
head(gsva_BLCA)  #GSVA result

# survival analysis --------------------------------------------------------------------


BLCA_sample <- intersect(colnames(gsva_BLCA),surdata$sample) 

gsva_BLCA2 <- gsva_BLCA[,BLCA_sample]
survival_BLCA <- subset(surdata,surdata$sample%in%BLCA_sample)

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
                   palette = c("blue","red"),#线的颜色
                   title = "BLCA",
                   legend.title=subpathway,legend=c(0.9,0.9),
                   conf.int = F)
    
    gg2 <- plot_grid(gg$plot + theme(legend.text = element_text(size = 12)), gg$table, ncol = 1, rel_widths = 1, rel_heights = c(1, .4), align = "hv")
    pdf(paste(paste("BLCA_sur",i,sep=""),".pdf",sep = ""))
    print(gg2)
    dev.off()
  }
}

#data.frame(BLCA_p,HNSC_p,KIRC_p,READ_p,LUAD_p,LUSC_p,KIRP_p,LIHC_p,THCA_p,CHOL_p,COAD_p,GBM_p,KICH_p,LGG_p,BRCA_p,UCS_p,UCEC_p,OV_p,PRAD_p,STAD_p,SKCM_p,CESC_p,ACC_p,PCPG_p,SARC_p,LAML_p,PAAD_p,ESCA_p,TGCT_p,THYM_p,MESO_p,UVM_p,DLBC_p)


# Univariate cox regression analysis ----------------------------------------------------------------

file_names<-list.files(pattern = "*.pdf")

library(stringr)
cancer_type<-str_extract(file_names,"\\b[A-Z]{1,4}")
subpathway_number<-str_extract(file_names,"[1-9]{1}")%>%as.numeric()

library(forestplot)
library(ggplot2)
library(ggstatsplot)
library(survival)
library(stringr)
library(scales)
library(viridis)

marker_gene<-mysubpathway%>%unlist()%>%unique()


tpm<-cbind(BLCA_tpm,HNSC_tpm,KIRC_tpm,READ_tpm,LUAD_tpm,LUSC_tpm,KIRP_tpm,LIHC_tpm,THCA_tpm,CHOL_tpm,COAD_tpm,GBM_tpm,KICH_tpm,LGG_tpm,BRCA_tpm,UCS_tpm,UCEC_tpm,OV_tpm,PRAD_tpm,STAD_tpm,SKCM_tpm,CESC_tpm,ACC_tpm,PCPG_tpm,SARC_tpm,LAML_tpm,PAAD_tpm,ESCA_tpm,TGCT_tpm,THYM_tpm,MESO_tpm,UVM_tpm,DLBC_tpm)
tpm2<-tpm[marker_gene,]


tpm3<-t(tpm2)%>%as.data.frame()
tpm3$sample <- rownames(tpm3)
surdata2<-surdata[,c("sample","OS","OS.time")]

cdata<-inner_join(tpm3,surdata2)
rownames(cdata)<-cdata$sample
cdata<-cdata[,-38] #把sample列变成行名
cdata<-cdata%>%select(OS.time,OS,everything())



# forestplot -------------------------------------------------------------

Coxoutput=data.frame()

for(i in colnames(cdata[,3:ncol(cdata)])){
  cox <- coxph(Surv(OS.time, OS) ~ cdata[,i], data = cdata)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}
Coxoutput <- arrange(Coxoutput,pvalue)   #按照p值排序 [37] pancancer_cox_result
Coxout <- filter(Coxoutput,Coxoutput$pvalue<0.05)  #[28] P_value significant

write.csv(Coxoutput,'f:/workplace/结果4/COX结果图/pancancer_cox_output.csv', row.names = F)

# Univariate cox regression analysis by cancer types -----------------------------------------------------------


cancer33<-unique(surdata$`cancer type abbreviation`)

cox_list<-list()
sig_cox_list<-list()

for(i in 1:length(cancer33)){
  s<-subset(surdata,surdata$`cancer type abbreviation`==cancer33[i])
  sam<-intersect(s$sample,rownames(cdata))
  
  cdat<-cdata[sam,]
  
  Coxoutput=data.frame()
  
  for(j in colnames(cdat[,3:ncol(cdat)])){
    cox <- coxph(Surv(OS.time, OS) ~ cdat[,j], data = cdat)
    coxSummary = summary(cox)
    Coxoutput=rbind(Coxoutput,cbind(gene=j,HR=coxSummary$coefficients[,"exp(coef)"],
                                    z=coxSummary$coefficients[,"z"],
                                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                    lower=coxSummary$conf.int[,3],
                                    upper=coxSummary$conf.int[,4]))
  }
  for(q in c(2:6)){
    Coxoutput[,q] <- as.numeric(as.vector(Coxoutput[,q]))
  }
  Coxoutput <- arrange(Coxoutput,pvalue) 
  write.csv(Coxoutput,file=paste(cancer33[i],"_cox_output.csv",sep=""),row.names = F)
  
  cox_list[[i]]<-Coxoutput
  names(cox_list)[i]=cancer33[i]
  
  
  Coxout <- filter(Coxoutput,Coxoutput$pvalue<0.05)
  write.csv(Coxout,file=paste(cancer33[i],"_cox_significant.csv",sep=""),row.names = F)

  sig_cox_list[[i]]<-Coxout
  names(sig_cox_list)[i]=cancer33[i]
  
}
cox_list%>%head()
sig_cox_list%>%head()

# Pan-cancer forest plot ----------------------------------------------------------


HR_input<-data.frame()
for(i in 1:length(marker_gene)){
  ge=marker_gene[i]
  d=data.frame()#某个基因在33癌症的cox结果
  
  for(j in 1:length(cancer33)){
    cancer1<-subset(cox_list[[j]],cox_list[[j]]$gene==ge)
    cancer1$cancer=cancer33[j]
    d<-rbind(d,cancer1)
  }
  HR_input<-rbind(HR_input,d)
 
}


for(i in 23:length(marker_gene)){
  
  df<-subset(HR_input,HR_input$gene==marker_gene[i])
  #某基因的HR输入数据框
  tabletext<-cbind(c("Cancer Type",NA,df$cancer,NA),
                   c("Pvalue",NA,round(df$pvalue,3),NA),
                   c("Hazard Ratio",NA,format(round(df$HR,2),nsmall=2)))
  pdf(paste(paste(marker_gene[i],"_forest",sep=""),".pdf",sep = ""))#新建PDF文件，准备写入图形
  forestplot(labeltext=tabletext,mean=c(NA,1,round(df$HR,3),NA),title=marker_gene[i],graph.pos=3,#图在表中的列位置
             graphwidth = unit(.3,"npc"),#图在表中的宽度比例
             lower=c(NA,1,round(df$lower,2),NA),upper=c(NA,1,round(df$upper,2),NA),col=fpColors(box="steelblue", lines="black", zero = "black"),#box颜色
             xlab= "Hazard Ratio (95% CI)",zero=1,boxsize=0.6,lwd.ci=1 ,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
             fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
             txt_gp=fpTxtGp(label=gpar(cex=0.8),#各种字体大小设置
                            ticks=gpar(cex=0.8),
                            xlab=gpar(cex = 0.8),
                            title=gpar(cex = 0.8)),
             new_page = F,#是否新页
             hrzl_lines=list("3" = gpar(lwd=1, col="black"),#第三行顶部加黑线，引号内数字标记行位置
                             "36" = gpar(lwd=1, col="black"))#最后一行底部加黑线
  )
  dev.off()
}

#  forest plot ----------------------------------------------------------

pancancer_cox_output<-read.csv("f:/workplace/结果4/COX结果图/pancancer_cox_output.csv")

df<-pancancer_cox_output
tabletext<-cbind(c("Gene Symbol",NA,df$gene,NA),
                 c("Pvalue",NA,round(df$pvalue,3),NA),
                 c("Hazard Ratio",NA,format(round(df$HR,2),nsmall=2)))

pdf("f:/workplace/结果4/HR森林图/pancancer_forest.pdf")#新建PDF文件，准备写入图形
forestplot(labeltext=tabletext,mean=c(NA,1,round(df$HR,3),NA),title="pancancer cox",graph.pos=3,#图在表中的列位置
           graphwidth = unit(.3,"npc"),#图在表中的宽度比例
           lower=c(NA,1,round(df$lower,2),NA),upper=c(NA,1,round(df$upper,2),NA),col=fpColors(box="steelblue", lines="black", zero = "black"),#box颜色
           xlab= "Hazard Ratio (95% CI)",zero=1,boxsize=0.6,lwd.ci=1 ,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
           txt_gp=fpTxtGp(label=gpar(cex=0.8),#各种字体大小设置
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex = 0.8),
                          title=gpar(cex = 0.8)),
           new_page = F,#是否新页
           hrzl_lines=list("3" = gpar(lwd=1, col="black"),#第三行顶部加黑线，引号内数字标记行位置
                           "40" = gpar(lwd=1, col="black"))#最后一行底部加黑线
)
dev.off()


HR_input2=data.frame()
#把sig_cox_list变成一个数据框HR_input2  
for(i in 1:length(cancer33)){
  if(dim(sig_cox_list[[i]])[1]!=0){
    cancer1<-sig_cox_list[[i]]
    cancer1$cancer=names(sig_cox_list)[i]
    HR_input2<-rbind(HR_input2,cancer1) 
  }
}
#32种癌症的森林图

cancer32<-HR_input2$cancer%>%unique()

for(i in 1:length(cancer32)){
  df<-subset(HR_input2 , HR_input2$cancer==cancer32[i])
  #某癌症的HR输入数据框
  tabletext<-cbind(c("Gene Symbol",NA,df$gene,NA),
                   c("Pvalue",NA,round(df$pvalue,3),NA),
                   c("Hazard Ratio",NA,format(round(df$HR,2),nsmall=2)))
  
  ll<-(dim(df)[1]+3)%>%as.character()
  pdf(paste(paste(cancer32[i],"_fore_sig",sep=""),".pdf",sep = ""))#新建PDF文件，准备写入图形
  forestplot(labeltext=tabletext,mean=c(NA,1,round(df$HR,3),NA),title=cancer32[i],graph.pos=3,#图在表中的列位置
             graphwidth = unit(.3,"npc"),#图在表中的宽度比例
             lower=c(NA,1,round(df$lower,2),NA),upper=c(NA,1,round(df$upper,2),NA),col=fpColors(box="steelblue", lines="black", zero = "black"),#box颜色
             xlab= "Hazard Ratio (95% CI)",zero=1,boxsize=0.4,lwd.ci=1 ,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
             fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
             txt_gp=fpTxtGp(label=gpar(cex=0.8),#各种字体大小设置
                            ticks=gpar(cex=0.8),
                            xlab=gpar(cex = 0.8),
                            title=gpar(cex = 0.8)),
             new_page = F,#是否新页
             hrzl_lines=list("3" = gpar(lwd=1, col="black"))#第三行顶部加黑线，引号内数字标记行位置
             
  )
  dev.off()#写入图型并关闭PDF文件
}
write.csv(HR_input,"f:/workplace/结果4/COX结果图/COX分析结果_分癌型.csv",row.names=F)


subpathway_out<-c()
for(i in 1:9){
  subpathway_out[[i]]<-paste(mysubpathway[[i]],collapse = " ")
}

subpathway_output<-data.frame(subpathway=names(mysubpathway),genelist=subpathway_out)
write.csv(subpathway_output,"f:/workplace/结果4/subpathway9.csv",row.names = F)

gsva_output<-cbind(gsva_BLCA,gsva_HNSC,gsva_KIRC,gsva_READ,gsva_LUAD,gsva_LUSC,gsva_KIRP,gsva_LIHC,gsva_THCA,gsva_CHOL,gsva_COAD,gsva_GBM,gsva_KICH,gsva_LGG,gsva_BRCA,gsva_UCS,gsva_UCEC,gsva_OV,gsva_PRAD,gsva_STAD,gsva_SKCM,gsva_CESC,gsva_ACC,gsva_PCPG,gsva_SARC,gsva_LAML,gsva_PAAD,gsva_ESCA,gsva_TGCT,gsva_THYM,gsva_MESO,gsva_UVM,gsva_DLBC)
write.csv(gsva_output,"f:/workplace/结果4/gsva_output.csv")


# survival analysis group by sex  ---------------------------------------------------------------


sex_diff_exp<-fread("sex_diff_exp.csv")%>%as.data.frame()  
colnames(sex_diff_exp)[1]<-"gene"

sex_diff_methy<-fread("sex_diff_methy.csv")%>%as.data.frame()  

sex_diff_exp2<-sex_diff_exp%>%dplyr::arrange(gene)
sex_diff_exp2$data<-rep("expression",dim(sex_diff_exp2)[1])

sex_df1<-sex_diff_exp2[,c("gene","bias","cancer","data")]

sex_diff_methy2<-sex_diff_methy%>%dplyr::arrange(gene)
sex_diff_methy2$data<-rep("methylation",dim(sex_diff_methy2)[1])

sex_df2<-sex_diff_methy2
colnames(sex_df2)[2]<-"bias"

sex_df<-rbind(sex_df1,sex_df2)
sex_df<-rbind(sex_df,c("TP53","Female","LIHC","mutation"))


sex_df_cancer<-arrange(sex_df,data,cancer,gene)


write.csv(sex_df_cancer ,"f:/workplace/结果4/sex_gene_cancer.csv",row.names=F)



sex_df_cancer_exp<-subset(sex_df_cancer,sex_df_cancer$data=="expression")

exp_both_gene<-intersect(sex_df_cancer_exp$gene ,subpathway_genelist)


cancer8<-sex_df_cancer_exp$cancer%>%unique()


surdata8<-subset(surdata,surdata$`cancer type abbreviation`%in%cancer8)

survival_cancer8<-surdata8[,c("sample","cancer type abbreviation","gender","OS","OS.time")]


setwd("f:/workplace/结果4/性别差异生存图/")
for(i in 1:length(cancer8)){
  sur1<-subset(survival_cancer8,survival_cancer8$`cancer type abbreviation`==cancer8[i])
  
  fit1 <- survfit(Surv(OS.time, OS) ~ gender, data = sur1)
  
  gg<-ggsurvplot(fit1, pval = TRUE,linetype = c("solid", "dashed"), 
                 palette = c("orange", "purple"),
                 legend.title=cancer8[i],
                 legend=c(0.2,0.1),
                 conf.int = F)
  
  gg2 <- plot_grid(gg$plot + theme(legend.text = element_text(size = 12)), gg$table, ncol = 1, rel_widths = 1, rel_heights = c(1, .4), align = "hv")
  pdf(paste(paste(cancer8[i],"_sex_survival",sep=""),".pdf",sep = ""))
  print(gg2)
  dev.off()
}



# Forest plot of ISC-pathways  --------------------------------------------------------------



library(data.table)
library(dplyr)


surdata<-fread("Survival_SupplementalTable_S1_20171025_xena_sp")%>%as.data.frame()
unique(surdata$`cancer type abbreviation`)
surdata<-subset(surdata,surdata$`cancer type abbreviation`=="HNSC")

library(survival)
library(survminer)
library(msigdbr)
library(tibble)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)


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


# gene expression ----------------------------------------------------------------
HNSC_tpm<-fread("f:/workplace/HNSC/HNSC_readcount.genes.tpm.txt")%>%as.data.frame()
rownames(HNSC_tpm)<-HNSC_tpm$V1
HNSC_tpm<-HNSC_tpm[,-1]

expMatrix_tpm<-(log2(HNSC_tpm+1))  #所有基因在HNSC的表达log2+1(TPM值)

# GSVA  ------------------------------------------------------------------
gsva_HNSC <- gsva(as.matrix(expMatrix_tpm), mysubpathway)
head(gsva_HNSC)  #GSVA result

gsva_HNSC2<-t(gsva_HNSC)%>%as.data.frame()
gsva_HNSC2$sample <- rownames(gsva_HNSC2)
surdata2<-surdata[,c("sample","OS","OS.time")]

cdata<-inner_join(gsva_HNSC2,surdata2)
rownames(cdata)<-cdata$sample
cdata<-cdata[,-10] #把sample列变成行名
cdata<-cdata%>%select(OS.time,OS,everything())

Coxoutput=data.frame()

for(i in colnames(cdata[,3:ncol(cdata)])){
  cox <- coxph(Surv(OS.time, OS) ~ cdata[,i], data = cdata)
  coxSummary = summary(cox)
  Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                  z=coxSummary$coefficients[,"z"],
                                  pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                  lower=coxSummary$conf.int[,3],
                                  upper=coxSummary$conf.int[,4]))
}
for(i in c(2:6)){
  Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
}
Coxoutput <- arrange(Coxoutput,pvalue)   #按照p值排序 [37] pancancer_cox_result
Coxout <- filter(Coxoutput,Coxoutput$pvalue<0.05)  #[28] P_value significant

write.csv(Coxoutput,'f:/workplace/结果4/HNSC_subpathway_coxoutput.csv', row.names = F)

