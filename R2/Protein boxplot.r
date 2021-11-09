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



# Pathway genes -------------------------------------------------------------

gene47<-read.table("f:/workplace/结果2/pathway_gene_47.txt")
path_immu_gene<-gene47$V1%>%as.character()


# 有RPPA数据的免疫衰老基因 ----------------------------------------------------------
BLCA_RPPA<-fread("f:/workplace/BLCA/BLCA_RPPA_RBN")%>%as.data.frame()%>%distinct()
RPPA_symbol<-bitr(geneID=BLCA_RPPA$Sample_description, fromType="ALIAS", toType=c("GENENAME","SYMBOL"), OrgDb="org.Hs.eg.db", drop = TRUE)
immu_protein<-intersect(RPPA_symbol$SYMBOL,path_immu_gene)
##"ATM","BCL2", "BCL2L11","MTOR", "TP53"
subset(RPPA_symbol,RPPA_symbol$SYMBOL%in%immu_protein)
####"ATM"  "BCL2" "BIM"  "MTOR" "P53"

# Protein data of 11 cancer types -------------------------------------------------------------

BLCA_RPPA<-fread("f:/workplace/BLCA/BLCA_RPPA_RBN")%>%as.data.frame()%>%distinct()
# 131*128
rownames(BLCA_RPPA)<-BLCA_RPPA$Sample_description
BLCA_RPPA<-BLCA_RPPA[,-1]
BLCA_protein<-t(BLCA_RPPA)

BRCA_RPPA<-fread("f:/workplace/BRCA/BRCA_RPPA_RBN")%>%as.data.frame()%>%distinct()
# 131*748
rownames(BRCA_RPPA)<-BRCA_RPPA$Sample_description
BRCA_RPPA<-BRCA_RPPA[,-1]
BRCA_protein<-t(BRCA_RPPA)


COAD_RPPA<-fread("f:/workplace/COAD/COAD_RPPA_RBN")%>%as.data.frame()%>%distinct()
# 131*335
rownames(COAD_RPPA)<-COAD_RPPA$Sample_description
COAD_RPPA<-COAD_RPPA[,-1]
COAD_protein<-t(COAD_RPPA)


GBM_RPPA<-fread("f:/workplace/GBM/GBM_RPPA_RBN")%>%as.data.frame()%>%distinct()
# 131*216
rownames(GBM_RPPA)<-GBM_RPPA$Sample_description
GBM_RPPA<-GBM_RPPA[,-1]
GBM_protein<-t(GBM_RPPA)


HNSC_RPPA<-fread("f:/workplace/HNSC/HNSC_RPPA_RBN")%>%as.data.frame()%>%distinct()
#131*213
rownames(HNSC_RPPA)<-HNSC_RPPA$Sample_description
HNSC_RPPA<-HNSC_RPPA[,-1]
HNSC_protein<-t(HNSC_RPPA)


KIRC_RPPA<-fread("f:/workplace/KIRC/KIRC_RPPA_RBN")%>%as.data.frame()%>%distinct()
#131*455
rownames(KIRC_RPPA)<-KIRC_RPPA$Sample_description
KIRC_RPPA<-KIRC_RPPA[,-1]
KIRC_protein<-t(KIRC_RPPA)


LUAD_RPPA<-fread("f:/workplace/LUAD/LUAD_RPPA_RBN")%>%as.data.frame()%>%distinct()
#131*238
rownames(LUAD_RPPA)<-LUAD_RPPA$Sample_description
LUAD_RPPA<-LUAD_RPPA[,-1]
LUAD_protein<-t(LUAD_RPPA)


LUSC_RPPA<-fread("f:/workplace/LUSC/LUSC_RPPA_RBN")%>%as.data.frame()%>%distinct()
#131*196
rownames(LUSC_RPPA)<-LUSC_RPPA$Sample_description
LUSC_RPPA<-LUSC_RPPA[,-1]
LUSC_protein<-t(LUSC_RPPA)


OV_RPPA<-fread("f:/workplace/OV/OV_RPPA_RBN")%>%as.data.frame()%>%distinct()
#131*413
rownames(OV_RPPA)<-OV_RPPA$Sample_description
OV_RPPA<-OV_RPPA[,-1]
OV_protein<-t(OV_RPPA)


READ_RPPA<-fread("f:/workplace/READ/READ_RPPA_RBN")%>%as.data.frame()%>%distinct()
#131*131
rownames(READ_RPPA)<-READ_RPPA$Sample_description
READ_RPPA<-READ_RPPA[,-1]
READ_protein<-t(READ_RPPA)


UCEC_RPPA<-fread("f:/workplace/UCEC/UCEC_RPPA_RBN")%>%as.data.frame()%>%distinct()
#131*405
rownames(UCEC_RPPA)<-UCEC_RPPA$Sample_description
UCEC_RPPA<-UCEC_RPPA[,-1]
UCEC_protein<-t(UCEC_RPPA)


# CNA group of 11 cancer types -------------------------------------------------------------

BLCA<-fread("f:/workplace/BLCA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(BLCA)<-BLCA$`Gene Symbol`
BLCA<-BLCA[,-1]

both_gene<-intersect(path_immu_gene,rownames(BLCA))

BLCA_CNV<-BLCA[both_gene,]
BLCA_cnv<-ifelse(BLCA_CNV<0,1,ifelse(BLCA_CNV==2,1,0))

########
BRCA<-fread("f:/workplace/BRCA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(BRCA)<-BRCA$`Gene Symbol`
BRCA<-BRCA[,-1]

both_gene<-intersect(path_immu_gene,rownames(BRCA))

BRCA_CNV<-BRCA[both_gene,]

BRCA_cnv<-ifelse(BRCA_CNV<0,1,ifelse(BRCA_CNV==2,1,0))

########
COAD<-fread("f:/workplace/COAD/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(COAD)<-COAD$`Gene Symbol`
COAD<-COAD[,-1]

both_gene<-intersect(path_immu_gene,rownames(COAD))

COAD_CNV<-COAD[both_gene,]

COAD_cnv<-ifelse(COAD_CNV<0,1,ifelse(COAD_CNV==2,1,0))

###########
GBM<-fread("f:/workplace/GBM/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(GBM)<-GBM$`Gene Symbol`
GBM<-GBM[,-1]

both_gene<-intersect(path_immu_gene,rownames(GBM))

GBM_CNV<-GBM[both_gene,]

GBM_cnv<-ifelse(GBM_CNV<0,1,ifelse(GBM_CNV==2,1,0))

#########
HNSC<-fread("f:/workplace/HNSC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(HNSC)<-HNSC$`Gene Symbol`
HNSC<-HNSC[,-1]

both_gene<-intersect(path_immu_gene,rownames(HNSC))

HNSC_CNV<-HNSC[both_gene,]

HNSC_cnv<-ifelse(HNSC_CNV<0,1,ifelse(HNSC_CNV==2,1,0))

######
KIRC<-fread("f:/workplace/KIRC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(KIRC)<-KIRC$`Gene Symbol`
KIRC<-KIRC[,-1]

both_gene<-intersect(path_immu_gene,rownames(KIRC))

KIRC_CNV<-KIRC[both_gene,]

KIRC_cnv<-ifelse(KIRC_CNV<0,1,ifelse(KIRC_CNV==2,1,0))

######
LUAD<-fread("f:/workplace/LUAD/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(LUAD)<-LUAD$`Gene Symbol`
LUAD<-LUAD[,-1]

both_gene<-intersect(path_immu_gene,rownames(LUAD))

LUAD_CNV<-LUAD[both_gene,]

LUAD_cnv<-ifelse(LUAD_CNV<0,1,ifelse(LUAD_CNV==2,1,0))

#######
LUSC<-fread("f:/workplace/LUSC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(LUSC)<-LUSC$`Gene Symbol`
LUSC<-LUSC[,-1]

both_gene<-intersect(path_immu_gene,rownames(LUSC))

LUSC_CNV<-LUSC[both_gene,]

LUSC_cnv<-ifelse(LUSC_CNV<0,1,ifelse(LUSC_CNV==2,1,0))

#####
OV<-fread("f:/workplace/OV/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(OV)<-OV$`Gene Symbol`
OV<-OV[,-1]

both_gene<-intersect(path_immu_gene,rownames(OV))

OV_CNV<-OV[both_gene,]

OV_cnv<-ifelse(OV_CNV<0,1,ifelse(OV_CNV==2,1,0))
#####
READ<-fread("f:/workplace/READ/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(READ)<-READ$`Gene Symbol`
READ<-READ[,-1]

both_gene<-intersect(path_immu_gene,rownames(READ))

READ_CNV<-READ[both_gene,]

READ_cnv<-ifelse(READ_CNV<0,1,ifelse(READ_CNV==2,1,0))

######
UCEC<-fread("f:/workplace/UCEC/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(UCEC)<-UCEC$`Gene Symbol`
UCEC<-UCEC[,-1]

both_gene<-intersect(path_immu_gene,rownames(UCEC))

UCEC_CNV<-UCEC[both_gene,]

UCEC_cnv<-ifelse(UCEC_CNV<0,1,ifelse(UCEC_CNV==2,1,0))


# ATM protein boxplot -----------------------------------------------------------

##### ATM protein expression
BLCA_ATM<-data.frame(sample=rownames(BLCA_protein),protein=BLCA_protein[,"ATM"])
BRCA_ATM<-data.frame(sample=rownames(BRCA_protein),protein=BRCA_protein[,"ATM"])
COAD_ATM<-data.frame(sample=rownames(COAD_protein),protein=COAD_protein[,"ATM"])
GBM_ATM<-data.frame(sample=rownames(GBM_protein),protein=GBM_protein[,"ATM"])
HNSC_ATM<-data.frame(sample=rownames(HNSC_protein),protein=HNSC_protein[,"ATM"])
KIRC_ATM<-data.frame(sample=rownames(KIRC_protein),protein=KIRC_protein[,"ATM"])
LUAD_ATM<-data.frame(sample=rownames(LUAD_protein),protein=LUAD_protein[,"ATM"])
LUSC_ATM<-data.frame(sample=rownames(LUSC_protein),protein=LUSC_protein[,"ATM"])
OV_ATM<-data.frame(sample=rownames(OV_protein),protein=OV_protein[,"ATM"])
READ_ATM<-data.frame(sample=rownames(READ_protein),protein=READ_protein[,"ATM"])
UCEC_ATM<-data.frame(sample=rownames(UCEC_protein),protein=UCEC_protein[,"ATM"])

##### ATM CNA group
BLCA_ATM_CNA<-data.frame(sample=colnames(BLCA_cnv),CNA_group=BLCA_cnv["ATM",])
BRCA_ATM_CNA<-data.frame(sample=colnames(BRCA_cnv),CNA_group=BRCA_cnv["ATM",])
COAD_ATM_CNA<-data.frame(sample=colnames(COAD_cnv),CNA_group=COAD_cnv["ATM",])
GBM_ATM_CNA<-data.frame(sample=colnames(GBM_cnv),CNA_group=GBM_cnv["ATM",])
HNSC_ATM_CNA<-data.frame(sample=colnames(HNSC_cnv),CNA_group=HNSC_cnv["ATM",])
KIRC_ATM_CNA<-data.frame(sample=colnames(KIRC_cnv),CNA_group=KIRC_cnv["ATM",])
LUAD_ATM_CNA<-data.frame(sample=colnames(LUAD_cnv),CNA_group=LUAD_cnv["ATM",])
LUSC_ATM_CNA<-data.frame(sample=colnames(LUSC_cnv),CNA_group=LUSC_cnv["ATM",])
OV_ATM_CNA<-data.frame(sample=colnames(OV_cnv),CNA_group=OV_cnv["ATM",])
READ_ATM_CNA<-data.frame(sample=colnames(READ_cnv),CNA_group=READ_cnv["ATM",])
UCEC_ATM_CNA<-data.frame(sample=colnames(UCEC_cnv),CNA_group=UCEC_cnv["ATM",])


###BLCA
BLCA_ATM_union<-inner_join(BLCA_ATM,BLCA_ATM_CNA)
BLCA_ATM_union$cancer<-rep("BLCA",dim(BLCA_ATM_union)[1])
BLCA_ATM_union$CNA_group<-ifelse(BLCA_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###BRCA
BRCA_ATM_union<-inner_join(BRCA_ATM,BRCA_ATM_CNA)
BRCA_ATM_union$cancer<-rep("BRCA",dim(BRCA_ATM_union)[1])
BRCA_ATM_union$CNA_group<-ifelse(BRCA_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###COAD
COAD_ATM_union<-inner_join(COAD_ATM,COAD_ATM_CNA)
COAD_ATM_union$cancer<-rep("COAD",dim(COAD_ATM_union)[1])
COAD_ATM_union$CNA_group<-ifelse(COAD_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###GBM
GBM_ATM_union<-inner_join(GBM_ATM,GBM_ATM_CNA)
GBM_ATM_union$cancer<-rep("GBM",dim(GBM_ATM_union)[1])
GBM_ATM_union$CNA_group<-ifelse(GBM_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###HNSC
HNSC_ATM_union<-inner_join(HNSC_ATM,HNSC_ATM_CNA)
HNSC_ATM_union$cancer<-rep("HNSC",dim(HNSC_ATM_union)[1])
HNSC_ATM_union$CNA_group<-ifelse(HNSC_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

####KIRC
KIRC_ATM_union<-inner_join(KIRC_ATM,KIRC_ATM_CNA)
KIRC_ATM_union$cancer<-rep("KIRC",dim(KIRC_ATM_union)[1])
KIRC_ATM_union$CNA_group<-ifelse(KIRC_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###LUAD
LUAD_ATM_union<-inner_join(LUAD_ATM,LUAD_ATM_CNA)
LUAD_ATM_union$cancer<-rep("LUAD",dim(LUAD_ATM_union)[1])
LUAD_ATM_union$CNA_group<-ifelse(LUAD_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###LUSC
LUSC_ATM_union<-inner_join(LUSC_ATM,LUSC_ATM_CNA)
LUSC_ATM_union$cancer<-rep("LUSC",dim(LUSC_ATM_union)[1])
LUSC_ATM_union$CNA_group<-ifelse(LUSC_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###OV
OV_ATM_union<-inner_join(OV_ATM,OV_ATM_CNA)
OV_ATM_union$cancer<-rep("OV",dim(OV_ATM_union)[1])
OV_ATM_union$CNA_group<-ifelse(OV_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###READ
READ_ATM_union<-inner_join(READ_ATM,READ_ATM_CNA)
READ_ATM_union$cancer<-rep("READ",dim(READ_ATM_union)[1])
READ_ATM_union$CNA_group<-ifelse(READ_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

###UCEC
UCEC_ATM_union<-inner_join(UCEC_ATM,UCEC_ATM_CNA)
UCEC_ATM_union$cancer<-rep("UCEC",dim(UCEC_ATM_union)[1])
UCEC_ATM_union$CNA_group<-ifelse(UCEC_ATM_union$CNA_group==1,"Amplified","Non-amplified")%>%as.factor()

ATM_box_input<-rbind(BLCA_ATM_union,BRCA_ATM_union,COAD_ATM_union,GBM_ATM_union,HNSC_ATM_union,KIRC_ATM_union,LUAD_ATM_union,LUSC_ATM_union,OV_ATM_union,READ_ATM_union,UCEC_ATM_union)
###3387sample

compare_means(protein ~ CNA_group, data = ATM_box_input,method = "t.test")

p <- ggboxplot(ATM_box_input, x = "CNA_group", y = "protein",
               add.params = list(color = "cancer",size=0.5), palette =  c(brewer.pal(9,"Set1"),"darkred","blue"),
               add = "jitter")
p + stat_compare_means(method = "t.test")+xlab("protein VS CNA (ATM)")

##### ATM mutation group


# 1 BLCA ------------------------------------------------------------------

BLCA<-fread("f:/workplace/BLCA/BLCA_mc3.txt")%>%as.data.frame()%>%distinct()
BLCA_mut<-subset(BLCA,effect!="Silent") #Remove silent mutation

BLCA_mut$sample%>%unique()%>%length()


BLCA_mut_47<-subset(BLCA_mut,BLCA_mut$gene %in% path_immu_gene)

BLCA_m<-table(BLCA_mut_47$sample,BLCA_mut_47$gene)

BLCA_ATM_mut<-data.frame(sample=rownames(BLCA_m),mut_group=BLCA_m[,"ATM"])
BLCA_ATM_mut$mut_group<-ifelse(BLCA_ATM_mut$mut_group>1,"Mutation","Non-mutation")


# 2 BRCA ------------------------------------------------------------------

BRCA<-fread("f:/workplace/BRCA/BRCA_mc3.txt")%>%as.data.frame()%>%distinct()
BRCA_mut<-subset(BRCA,effect!="Silent") 

BRCA_mut$sample%>%unique()%>%length()


BRCA_mut_47<-subset(BRCA_mut,BRCA_mut$gene %in% path_immu_gene)

BRCA_m<-table(BRCA_mut_47$sample,BRCA_mut_47$gene)

BRCA_ATM_mut<-data.frame(sample=rownames(BRCA_m),mut_group=BRCA_m[,"ATM"])
BRCA_ATM_mut$mut_group<-ifelse(BRCA_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 3 COAD ------------------------------------------------------------------

COAD<-fread("f:/workplace/COAD/COAD_mc3.txt")%>%as.data.frame()%>%distinct()
COAD_mut<-subset(COAD,effect!="Silent")

COAD_mut$sample%>%unique()%>%length()


COAD_mut_47<-subset(COAD_mut,COAD_mut$gene %in% path_immu_gene)

COAD_m<-table(COAD_mut_47$sample,COAD_mut_47$gene)

COAD_ATM_mut<-data.frame(sample=rownames(COAD_m),mut_group=COAD_m[,"ATM"])
COAD_ATM_mut$mut_group<-ifelse(COAD_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 4 GBM ------------------------------------------------------------------

GBM<-fread("f:/workplace/GBM/GBM_mc3.txt")%>%as.data.frame()%>%distinct()
GBM_mut<-subset(GBM,effect!="Silent") 

GBM_mut$sample%>%unique()%>%length()


GBM_mut_47<-subset(GBM_mut,GBM_mut$gene %in% path_immu_gene)

GBM_m<-table(GBM_mut_47$sample,GBM_mut_47$gene)

GBM_ATM_mut<-data.frame(sample=rownames(GBM_m),mut_group=GBM_m[,"ATM"])
GBM_ATM_mut$mut_group<-ifelse(GBM_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 5 HNSC ------------------------------------------------------------------

HNSC<-fread("f:/workplace/HNSC/HNSC_mc3.txt")%>%as.data.frame()%>%distinct()
HNSC_mut<-subset(HNSC,effect!="Silent") 

HNSC_mut$sample%>%unique()%>%length()


HNSC_mut_47<-subset(HNSC_mut,HNSC_mut$gene %in% path_immu_gene)

HNSC_m<-table(HNSC_mut_47$sample,HNSC_mut_47$gene)

HNSC_ATM_mut<-data.frame(sample=rownames(HNSC_m),mut_group=HNSC_m[,"ATM"])
HNSC_ATM_mut$mut_group<-ifelse(HNSC_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 6 KIRC ------------------------------------------------------------------

KIRC<-fread("f:/workplace/KIRC/KIRC_mc3.txt")%>%as.data.frame()%>%distinct()
KIRC_mut<-subset(KIRC,effect!="Silent") 

KIRC_mut$sample%>%unique()%>%length()


KIRC_mut_47<-subset(KIRC_mut,KIRC_mut$gene %in% path_immu_gene)

KIRC_m<-table(KIRC_mut_47$sample,KIRC_mut_47$gene)

KIRC_ATM_mut<-data.frame(sample=rownames(KIRC_m),mut_group=KIRC_m[,"ATM"])
KIRC_ATM_mut$mut_group<-ifelse(KIRC_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 7 LUAD ------------------------------------------------------------------

LUAD<-fread("f:/workplace/LUAD/LUAD_mc3.txt")%>%as.data.frame()%>%distinct()
LUAD_mut<-subset(LUAD,effect!="Silent") 

LUAD_mut$sample%>%unique()%>%length()


LUAD_mut_47<-subset(LUAD_mut,LUAD_mut$gene %in% path_immu_gene)

LUAD_m<-table(LUAD_mut_47$sample,LUAD_mut_47$gene)

LUAD_ATM_mut<-data.frame(sample=rownames(LUAD_m),mut_group=LUAD_m[,"ATM"])
LUAD_ATM_mut$mut_group<-ifelse(LUAD_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 8 LUSC ------------------------------------------------------------------

LUSC<-fread("f:/workplace/LUSC/LUSC_mc3.txt")%>%as.data.frame()%>%distinct()
LUSC_mut<-subset(LUSC,effect!="Silent") 

LUSC_mut$sample%>%unique()%>%length()


LUSC_mut_47<-subset(LUSC_mut,LUSC_mut$gene %in% path_immu_gene)

LUSC_m<-table(LUSC_mut_47$sample,LUSC_mut_47$gene)

LUSC_ATM_mut<-data.frame(sample=rownames(LUSC_m),mut_group=LUSC_m[,"ATM"])
LUSC_ATM_mut$mut_group<-ifelse(LUSC_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 9 OV ------------------------------------------------------------------

OV<-fread("f:/workplace/OV/OV_mc3.txt")%>%as.data.frame()%>%distinct()
OV_mut<-subset(OV,effect!="Silent") 

OV_mut$sample%>%unique()%>%length()


OV_mut_47<-subset(OV_mut,OV_mut$gene %in% path_immu_gene)

OV_m<-table(OV_mut_47$sample,OV_mut_47$gene)

OV_ATM_mut<-data.frame(sample=rownames(OV_m),mut_group=OV_m[,"ATM"])
OV_ATM_mut$mut_group<-ifelse(OV_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 10 READ ------------------------------------------------------------------

READ<-fread("f:/workplace/READ/READ_mc3.txt")%>%as.data.frame()%>%distinct()
READ_mut<-subset(READ,effect!="Silent") 

READ_mut$sample%>%unique()%>%length()


READ_mut_47<-subset(READ_mut,READ_mut$gene %in% path_immu_gene)

READ_m<-table(READ_mut_47$sample,READ_mut_47$gene)

READ_ATM_mut<-data.frame(sample=rownames(READ_m),mut_group=READ_m[,"ATM"])
READ_ATM_mut$mut_group<-ifelse(READ_ATM_mut$mut_group>1,"Mutation","Non-mutation")

# 11 UCEC ------------------------------------------------------------------

UCEC<-fread("f:/workplace/UCEC/UCEC_mc3.txt")%>%as.data.frame()%>%distinct()
UCEC_mut<-subset(UCEC,effect!="Silent")

UCEC_mut$sample%>%unique()%>%length()


UCEC_mut_47<-subset(UCEC_mut,UCEC_mut$gene %in% path_immu_gene)

UCEC_m<-table(UCEC_mut_47$sample,UCEC_mut_47$gene)

UCEC_ATM_mut<-data.frame(sample=rownames(UCEC_m),mut_group=UCEC_m[,"ATM"])
UCEC_ATM_mut$mut_group<-ifelse(UCEC_ATM_mut$mut_group>1,"Mutation","Non-mutation")

#######Protein data & mutation data

###BLCA
BLCA_ATM_union<-inner_join(BLCA_ATM,BLCA_ATM_mut)
BLCA_ATM_union$cancer<-rep("BLCA",dim(BLCA_ATM_union)[1])

###BRCA
BRCA_ATM_union<-inner_join(BRCA_ATM,BRCA_ATM_mut)
BRCA_ATM_union$cancer<-rep("BRCA",dim(BRCA_ATM_union)[1])

###COAD
COAD_ATM_union<-inner_join(COAD_ATM,COAD_ATM_mut)
COAD_ATM_union$cancer<-rep("COAD",dim(COAD_ATM_union)[1])

###GBM
GBM_ATM_union<-inner_join(GBM_ATM,GBM_ATM_mut)
GBM_ATM_union$cancer<-rep("GBM",dim(GBM_ATM_union)[1])


###HNSC
HNSC_ATM_union<-inner_join(HNSC_ATM,HNSC_ATM_mut)
HNSC_ATM_union$cancer<-rep("HNSC",dim(HNSC_ATM_union)[1])

####KIRC
KIRC_ATM_union<-inner_join(KIRC_ATM,KIRC_ATM_mut)
KIRC_ATM_union$cancer<-rep("KIRC",dim(KIRC_ATM_union)[1])

###LUAD
LUAD_ATM_union<-inner_join(LUAD_ATM,LUAD_ATM_mut)
LUAD_ATM_union$cancer<-rep("LUAD",dim(LUAD_ATM_union)[1])

###LUSC
LUSC_ATM_union<-inner_join(LUSC_ATM,LUSC_ATM_mut)
LUSC_ATM_union$cancer<-rep("LUSC",dim(LUSC_ATM_union)[1])

###OV
OV_ATM_union<-inner_join(OV_ATM,OV_ATM_mut)
OV_ATM_union$cancer<-rep("OV",dim(OV_ATM_union)[1])

###READ
READ_ATM_union<-inner_join(READ_ATM,READ_ATM_mut)
READ_ATM_union$cancer<-rep("READ",dim(READ_ATM_union)[1])

###UCEC
UCEC_ATM_union<-inner_join(UCEC_ATM,UCEC_ATM_mut)
UCEC_ATM_union$cancer<-rep("UCEC",dim(UCEC_ATM_union)[1])


ATM_box_input<-rbind(BLCA_ATM_union,BRCA_ATM_union,COAD_ATM_union,GBM_ATM_union,HNSC_ATM_union,KIRC_ATM_union,LUAD_ATM_union,LUSC_ATM_union,OV_ATM_union,READ_ATM_union,UCEC_ATM_union)
###3387sample

compare_means(protein ~ mut_group, data = ATM_box_input,method = "t.test")

p <- ggboxplot(ATM_box_input, x = "mut_group", y = "protein",
               add.params = list(color = "cancer",size=0.5), palette =  c(brewer.pal(9,"Set1"),"darkred","blue"),
               add = "jitter")
p + stat_compare_means(method = "t.test")+xlab("(ATM)")
