
library(stringr)
library(ggplot2)
library("RColorBrewer")
display.brewer.all()
brewer.pal.info

LUAD_tpm<-read.table("F:/workplace/LUAD/LUAD_readcount.genes.tpm.txt")

LUAD_tpm<-(log2(LUAD_tpm+1))  #log2(TPM+1)
colnames(LUAD_tpm)<-str_replace_all(colnames(LUAD_tpm),"[.*]",replacement = "-")

gene90<-read.table("f:/workplace/GeneSymbol_90.txt")%>%unique()
gene90<-as.character(gene90$V1)  

LUAD_gene<-intersect(rownames(LUAD_tpm),gene90) 
LUAD_tpm2<-LUAD_tpm[LUAD_gene,] 

sample_group<-ifelse(as.numeric(str_sub(colnames(LUAD_tpm),14,15))<10,"Tumor","Normal")
sam_df<-data.frame(samID=colnames(LUAD_tpm),samGP=sample_group)
sam_df2<-subset(sam_df ,sam_df$samGP=="Tumor") 

LUAD_tpm_tumor<-LUAD_tpm2[,as.character(sam_df2$samID)] 
LUAD_immune_all<-apply(LUAD_tpm_tumor,2,sum)

LUAD_box<-data.frame(expSum=LUAD_immune_all,cancerType=rep("LUAD",length(LUAD_immune_all)))
#Sum of gene expressions



pancancer_box<-rbind(LUAD_box,BLCA_box,HNSC_box,KIRC_box,KIRP_box,LIHC_box,LUSC_box,THCA_box,COAD_box,GBM_box,KICH_box,LGG_box,BRCA_box,READ_box,UCS_box,UCEC_box,OV_box,PRAD_box,STAD_box,SKCM_box,CESC_box,ACC_box,PCPG_box,SARC_box,LAML_box,PAAD_box,ESCA_box,TGCT_box,THYM_box,MESO_box,UVM_box,DLBC_box,CHOL_box)


sortp<-pancancer_box %>%
  group_by(cancerType) %>%
  dplyr::summarise(expSum=median(expSum))

cancersort<-sortp$cancerType[order(sortp$expSum)]%>%as.character() #按照中值将癌型排序

mycol<-c(brewer.pal(12,"Set3"),brewer.pal(9,"Set1"),brewer.pal(12,"Paired"))
#33 color

ggplot(pancancer_box,aes(cancerType,expSum,fill=cancerType))+geom_boxplot()+theme_minimal()+# 不要背景
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12),
        axis.text.x = element_text(size=10,angle = 90, hjust = 1),
        axis.text.y = element_text(size=8))+
  ylab("Expression of immunosenescence gene")+
  scale_fill_manual(values = mycol) +
  scale_x_discrete(limits = cancersort)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

