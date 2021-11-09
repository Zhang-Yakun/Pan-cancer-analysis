library(stringr)
library(ggplot2)
library("RColorBrewer")
library(data.table)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)


# Pick the immunosenescence genes in the pathway ------------------------------------------------------------

hsa16<-read.csv("f:/workplace/enrichKEGG16.csv",header=T) 
pathway<-hsa16$ID%>%as.character()  #16条通路的ID向量

hsa_gene<-read.csv("f:/workplace/HSA_KEGG.csv",header=T)  

pathway16<-dplyr::filter(hsa_gene,hsa_gene$KEGGID%in%pathway)  #pathway genes
pathway_gene<-bitr(geneID=pathway16$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(pathway_gene)<-c("ENTREZID","gene") %>%unique() #gene ENTRZID

pathway_gene_result<-pathway_gene$gene%>%unique()

gene90<-read.table("f:/workplace/GeneSymbol_90.txt")%>%unique()
gene90<-as.character(gene90$V1)  #90 immunosenescence genes 

path_immu_gene<-intersect(pathway_gene_result,gene90)

write.table(path_immu_gene,"pathway_gene_47.txt",row.names = F,col.names = F,quote = F)

# The expression profiles of 47 genes in generalized carcinoma were screened ----------------------------------------------------------


LUAD_tpm<-read.table("F:/workplace/LUAD/LUAD_readcount.genes.tpm.txt")

LUAD_tpm<-(log2(LUAD_tpm+1))  #log2(TPM+1)
colnames(LUAD_tpm)<-str_replace_all(colnames(LUAD_tpm),"[.*]",replacement = "-")

LUAD_gene<-intersect(rownames(LUAD_tpm),gene90) 
LUAD_tpm2<-LUAD_tpm[LUAD_gene,]  

sample_group<-ifelse(as.numeric(str_sub(colnames(LUAD_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(LUAD_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor")

LUAD_tpm_tumor<-LUAD_tpm2[,as.character(sam_df2$samID)] 
#Immunosenescence gene expression profile after removing normal samples

identical(intersect(path_immu_gene,rownames(LUAD_tpm_tumor)),path_immu_gene)

LUAD_517<- LUAD_tpm_tumor[path_immu_gene,]

LUAD_immune_47<-apply(LUAD_517,2,sum)

LUAD_box_47<-data.frame(expSum=LUAD_immune_47,cancerType=rep("LUAD",length(LUAD_immune_47)))

# 2 BLCA -----------------------------------------------------------------------


BLCA_tpm<-read.table("F:/workplace/BLCA/BLCA_readcount.genes.tpm.txt")

BLCA_tpm<-(log2(BLCA_tpm+1)) 
colnames(BLCA_tpm)<-str_replace_all(colnames(BLCA_tpm),"[.*]",replacement = "-")

BLCA_gene<-intersect(rownames(BLCA_tpm),gene90) 
BLCA_tpm2<-BLCA_tpm[BLCA_gene,]  

sample_group<-ifelse(as.numeric(str_sub(colnames(BLCA_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(BLCA_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

BLCA_tpm_tumor<-BLCA_tpm2[,as.character(sam_df2$samID)]

identical(intersect(path_immu_gene,rownames(BLCA_tpm_tumor)),path_immu_gene)

BLCA_407 <- BLCA_tpm_tumor[path_immu_gene,]


BLCA_immune_47<-apply(BLCA_407,2,sum)

BLCA_box_47<-data.frame(expSum=BLCA_immune_47,cancerType=rep("BLCA",length(BLCA_immune_47)))

# 3HNSC -----------------------------------------------------------------------

HNSC_tpm<-read.table("F:/workplace/HNSC/HNSC_readcount.genes.tpm.txt")

HNSC_tpm<-(log2(HNSC_tpm+1))  
colnames(HNSC_tpm)<-str_replace_all(colnames(HNSC_tpm),"[.*]",replacement = "-")

HNSC_gene<-intersect(rownames(HNSC_tpm),gene90) 
HNSC_tpm2<-HNSC_tpm[HNSC_gene,] 

sample_group<-ifelse(as.numeric(str_sub(colnames(HNSC_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(HNSC_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

HNSC_tpm_tumor<-HNSC_tpm2[,as.character(sam_df2$samID)]

identical(intersect(path_immu_gene,rownames(HNSC_tpm_tumor)),path_immu_gene)

HNSC_522 <- HNSC_tpm_tumor[path_immu_gene,]


HNSC_immune_47<-apply(HNSC_522,2,sum)

HNSC_box_47<-data.frame(expSum=HNSC_immune_47,cancerType=rep("HNSC",length(HNSC_immune_47)))

# 4KIRC -----------------------------------------------------------------------

KIRC_tpm<-read.table("F:/workplace/KIRC/KIRC_readcount.genes.tpm.txt")

KIRC_tpm<-(log2(KIRC_tpm+1))  
colnames(KIRC_tpm)<-str_replace_all(colnames(KIRC_tpm),"[.*]",replacement = "-")

KIRC_gene<-intersect(rownames(KIRC_tpm),gene90) 
KIRC_tpm2<-KIRC_tpm[KIRC_gene,] 

sample_group<-ifelse(as.numeric(str_sub(colnames(KIRC_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(KIRC_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

KIRC_tpm_tumor<-KIRC_tpm2[,as.character(sam_df2$samID)]

identical(intersect(path_immu_gene,rownames(KIRC_tpm_tumor)),path_immu_gene)

KIRC_534 <- KIRC_tpm_tumor[path_immu_gene,]


KIRC_immune_47<-apply(KIRC_534,2,sum)

KIRC_box_47<-data.frame(expSum=KIRC_immune_47,cancerType=rep("KIRC",length(KIRC_immune_47)))

# 5 KIRP --------------------------------------------------------------------

KIRP_tpm<-read.table("F:/workplace/KIRP/KIRP_readcount.genes.tpm.txt")

KIRP_tpm<-(log2(KIRP_tpm+1)) 
colnames(KIRP_tpm)<-str_replace_all(colnames(KIRP_tpm),"[.*]",replacement = "-")

KIRP_gene<-intersect(rownames(KIRP_tpm),gene90) 
KIRP_tpm2<-KIRP_tpm[KIRP_gene,]  

sample_group<-ifelse(as.numeric(str_sub(colnames(KIRP_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(KIRP_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

KIRP_tpm_tumor<-KIRP_tpm2[,as.character(sam_df2$samID)] 

identical(intersect(path_immu_gene,rownames(KIRP_tpm_tumor)),path_immu_gene)

dim(KIRP_tpm_tumor)

KIRP_291 <- KIRP_tpm_tumor[path_immu_gene,]


KIRP_immune_47<-apply(KIRP_291,2,sum)

KIRP_box_47<-data.frame(expSum=KIRP_immune_47,cancerType=rep("KIRP",length(KIRP_immune_47)))

# 6 LIHC ------------------------------------------------------------------

LIHC_tpm<-read.table("F:/workplace/LIHC/LIHC_readcount.genes.tpm.txt")

LIHC_tpm<-(log2(LIHC_tpm+1)) 
colnames(LIHC_tpm)<-str_replace_all(colnames(LIHC_tpm),"[.*]",replacement = "-")

LIHC_gene<-intersect(rownames(LIHC_tpm),gene90) 
LIHC_tpm2<-LIHC_tpm[LIHC_gene,]  

sample_group<-ifelse(as.numeric(str_sub(colnames(LIHC_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(LIHC_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

LIHC_tpm_tumor<-LIHC_tpm2[,as.character(sam_df2$samID)] 


identical(intersect(path_immu_gene,rownames(LIHC_tpm_tumor)),path_immu_gene)

dim(LIHC_tpm_tumor)

LIHC_373 <- LIHC_tpm_tumor[path_immu_gene,]


LIHC_immune_47<-apply(LIHC_373,2,sum)

LIHC_box_47<-data.frame(expSum=LIHC_immune_47,cancerType=rep("LIHC",length(LIHC_immune_47)))

# 7 LUSC ------------------------------------------------------------------

LUSC_tpm<-read.table("F:/workplace/LUSC/LUSC_readcount.genes.tpm.txt")

LUSC_tpm<-(log2(LUSC_tpm+1))  
colnames(LUSC_tpm)<-str_replace_all(colnames(LUSC_tpm),"[.*]",replacement = "-")

LUSC_gene<-intersect(rownames(LUSC_tpm),gene90) 
LUSC_tpm2<-LUSC_tpm[LUSC_gene,]  

sample_group<-ifelse(as.numeric(str_sub(colnames(LUSC_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(LUSC_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

LUSC_tpm_tumor<-LUSC_tpm2[,as.character(sam_df2$samID)] 


identical(intersect(path_immu_gene,rownames(LUSC_tpm_tumor)),path_immu_gene)

dim(LUSC_tpm_tumor)

LUSC_502 <- LUSC_tpm_tumor[path_immu_gene,]

LUSC_immune_47<-apply(LUSC_502,2,sum)

LUSC_box_47<-data.frame(expSum=LUSC_immune_47,cancerType=rep("LUSC",length(LUSC_immune_47)))

# 8 THCA ------------------------------------------------------------------

THCA_tpm<-read.table("F:/workplace/THCA/THCA_readcount.genes.tpm.txt")

THCA_tpm<-(log2(THCA_tpm+1))  
colnames(THCA_tpm)<-str_replace_all(colnames(THCA_tpm),"[.*]",replacement = "-")

THCA_gene<-intersect(rownames(THCA_tpm),gene90) 
THCA_tpm2<-THCA_tpm[THCA_gene,]  

sample_group<-ifelse(as.numeric(str_sub(colnames(THCA_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(THCA_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

THCA_tpm_tumor<-THCA_tpm2[,as.character(sam_df2$samID)] 


identical(intersect(path_immu_gene,rownames(THCA_tpm_tumor)),path_immu_gene)

dim(THCA_tpm_tumor)

THCA_513 <- THCA_tpm_tumor[path_immu_gene,]


THCA_immune_47<-apply(THCA_513,2,sum)

THCA_box_47<-data.frame(expSum=THCA_immune_47,cancerType=rep("THCA",length(THCA_immune_47)))

# 9 COAD ------------------------------------------------------------------

COAD_tpm<-read.table("F:/workplace/COAD/COAD_readcount.genes.tpm.txt")

COAD_tpm<-(log2(COAD_tpm+1))  
colnames(COAD_tpm)<-str_replace_all(colnames(COAD_tpm),"[.*]",replacement = "-")

COAD_gene<-intersect(rownames(COAD_tpm),gene90) 
COAD_tpm2<-COAD_tpm[COAD_gene,] 

sample_group<-ifelse(as.numeric(str_sub(colnames(COAD_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(COAD_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

COAD_tpm_tumor<-COAD_tpm2[,as.character(sam_df2$samID)] 


identical(intersect(path_immu_gene,rownames(COAD_tpm_tumor)),path_immu_gene)

dim(COAD_tpm_tumor)

COAD_288 <- COAD_tpm_tumor[path_immu_gene,]


COAD_immune_47<-apply(COAD_288,2,sum)

COAD_box_47<-data.frame(expSum=COAD_immune_47,cancerType=rep("COAD",length(COAD_immune_47)))

# 10 GBM-----------------------------------------------------------------

GBM_tpm<-read.table("F:/workplace/GBM/GBM_readcount.genes.tpm.txt")

GBM_tpm<-(log2(GBM_tpm+1)) 
colnames(GBM_tpm)<-str_replace_all(colnames(GBM_tpm),"[.*]",replacement = "-")

GBM_gene<-intersect(rownames(GBM_tpm),gene90) 
GBM_tpm2<-GBM_tpm[GBM_gene,] 

sample_group<-ifelse(as.numeric(str_sub(colnames(GBM_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(GBM_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

GBM_tpm_tumor<-GBM_tpm2[,as.character(sam_df2$samID)] 


identical(intersect(path_immu_gene,rownames(GBM_tpm_tumor)),path_immu_gene)

dim(GBM_tpm_tumor)

GBM_167 <- GBM_tpm_tumor[path_immu_gene,]


GBM_immune_47<-apply(GBM_167,2,sum)

GBM_box_47<-data.frame(expSum=GBM_immune_47,cancerType=rep("GBM",length(GBM_immune_47)))



# Immunosenescence gene expression in pan-cancer heatmap -------------------------------------------------------

pancancer_heat<-cbind(LUAD_517,BLCA_407,HNSC_522,KIRC_534,KIRP_291,LIHC_373,LUSC_502,THCA_513,COAD_288,GBM_167,
                      KICH_66,LGG_530,BRCA_1104,READ_95,UCS_57,UCEC_177,OV_308,PRAD_498,STAD_415,SKCM_473,CESC_305,ACC_79,
                      PCPG_184,SARC_263,LAML_173,PAAD_179,ESCA_185,TGCT_156,THYM_120,MESO_87,UVM_80,DLBC_48,CHOL_36)


cormat<-round(cor(pancancer_heat , method = "pearson"),2)

heat_sample<-colnames(pancancer_heat)
heat_group<-rep(c("LUAD-517","BLCA-407","HNSC-522","KIRC-534","KIRP-291","LIHC-373","LUSC-502","THCA-513","COAD-288","GBM-167","KICH-66","LGG-530","BRCA-1104","READ-95","UCS-57","UCEC-177","OV-308","PRAD-498","STAD-415","SKCM-473","CESC-305","ACC-79","PCPG-184","SARC-263","LAML-173","PAAD-179","ESCA-185","TGCT-156","THYM-120","MESO-87","UVM-80","DLBC-48","CHOL-36"),
                c(517,407,522,534,291,373,502,513,288,167,66,530,1104,95,57,177,308,498,415,473,305,79,184,263,173,179,185,156,120,87,80,48,36))
heat_cov<-data.frame(sample=heat_sample,CancerType=heat_group)

annotation_col<-data.frame(CancerType = as.factor(heat_cov$CancerType))
rownames(annotation_col) = colnames(pancancer_heat)                          

mycol<-c(brewer.pal(12,"Set3"),brewer.pal(9,"Set1"),brewer.pal(12,"Paired"))

ann_colors = list(
  CancerType = c(`LUAD-517`="#FFFF99",`BLCA-407`="#386CB0",`HNSC-522`="#F0027F",`KIRC-534`="#BF5B17",`KIRP-291`="#A6CEE3",`LIHC-373`="#1F78B4",`LUSC-502`="#B2DF8A",`THCA-513`="#33A02C",`COAD-288`="#FB9A99",`GBM-167`="#E31A1C",`KICH-66`="#FDBF6F",`LGG-530`="#FF7F00",
                 `BRCA-1104`="#CAB2D6",`READ-95`="#6A3D9A",`UCS-57`="#FFFF99",`UCEC-177`="#B15928",`OV-308`="#1B9E77",`PRAD-498`="#D95F02",`STAD-415`="#7570B3",`SKCM-473`="#E7298A",`CESC-305`="#66A61E",`ACC-79`="#E6AB02",`PCPG-184`="#A6761D",`SARC-263`="#666666",`LAML-173`="#E41A1C",`PAAD-179`="#377EB8",`ESCA-185`="#4DAF4A",`TGCT-156`="#984EA3",`THYM-120`="#FF7F00",`MESO-87`="#FFFF33",`UVM-80`="#A65628",`DLBC-48`="#F781BF",`CHOL-36`="#999999")
)

n <- t(scale(t(pancancer_heat)))
n[n > 2] <- 2
n[n < -2] <- -2
df <- n
rownames(df)<-rownames(pancancer_heat)


pheatmap(df,cellwidth =0.05, cellheight = 6, fontsize = 8,fontsize_row=6,
         method="pearson", #"pearson" (default), "kendall", or "spearman"
         scale="row", #scale
         cluster_rows=T,
         cluster_cols=F,
         color = colorRampPalette(c("green", "black", "red"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         treeheight_col = "0",
         treeheight_row = "0",
         #legend = F,
         border_color = "NA")



# Immunosenescence score heatmap----------------------------------------------------------------

immune_all_score<-apply(pancancer_heat,2,sum)%>%as.data.frame()%>%t()
identical(colnames(immune_all_score),colnames(pancancer_heat))
rownames(immune_all_score)<-"immunosenescence"

n <- t(scale(t(immune_all_score)))
n[n > 2] <- 2
n[n < -2] <- -2
df <- n
rownames(df)<-rownames(immune_all_score)

pheatmap(df,cellwidth =0.05, cellheight = 8, fontsize = 8,fontsize_row=6,
         method="pearson", #"pearson" (default), "kendall", or "spearman"
         scale="row", #scale
         cluster_rows=F,
         cluster_cols=F,
         color = colorRampPalette(c("blue", "black", "yellow"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         treeheight_col = "0",
         treeheight_row = "0",
         #legend = F,
         border_color = "NA")


# Immunosenescence score boxplot ---------------------------------------------------------

pancancer_box_47<-rbind(LUAD_box_47,BLCA_box_47,HNSC_box_47,KIRC_box_47,KIRP_box_47,LIHC_box_47,LUSC_box_47,THCA_box_47,COAD_box_47,GBM_box_47,KICH_box_47,LGG_box_47,BRCA_box_47,READ_box_47,UCS_box_47,UCEC_box_47,OV_box_47,PRAD_box_47,STAD_box_47,SKCM_box_47,CESC_box_47,ACC_box_47,PCPG_box_47,SARC_box_47,LAML_box_47,PAAD_box_47,ESCA_box_47,TGCT_box_47,THYM_box_47,MESO_box_47,UVM_box_47,DLBC_box_47,CHOL_box_47)


sortp<-pancancer_box_47 %>%
  group_by(cancerType) %>%
  dplyr::summarise(expSum=median(expSum))

cancersort<-sortp$cancerType[order(sortp$expSum)]%>%as.character() #按照中值将癌型排序

mycol<-c(brewer.pal(12,"Set3"),brewer.pal(9,"Set1"),brewer.pal(12,"Paired"))
#33 color

ggplot(pancancer_box_47,aes(cancerType,expSum,fill=cancerType))+geom_boxplot()+theme_minimal()+# 不要背景
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 14),
        axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
  #axis.text.y = element_text(size=8))+
  ylab("Immunosenescence score")+
  scale_fill_manual(values = mycol) +
  scale_x_discrete(limits = cancersort)

