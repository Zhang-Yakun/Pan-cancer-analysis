
library(stringr)
library(ggplot2)
library("RColorBrewer")
library(data.table)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# 47 pathway immunosenescence genes -------------------------------------------------------------

gene47<-read.table("f:/workplace/结果2/pathway_gene_47.txt")
path_immu_gene<-gene47$V1%>%as.character()

# # Copy number variation data of 33 cancer types---------------------------------------------------------

# 1 LUAD ------------------------------------------------------------------
LUAD<-fread("f:/workplace/LUAD/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(LUAD)<-LUAD$`Gene Symbol`
LUAD<-LUAD[,-1]

both_gene<-intersect(path_immu_gene,rownames(LUAD))

LUAD_CNV<-LUAD[both_gene,]

LUAD_amp<-ifelse(LUAD_CNV>0,1,0) #Copy number amplification
LUAD_loss<-ifelse(LUAD_CNV<0,1,0) #Missing copy number

LUAD_high_amp<-ifelse(LUAD_CNV>1,1,0)
LUAD_high_amp2<-data.frame(gene=rownames(LUAD_high_amp),LUAD_high_CNV=apply(LUAD_high_amp,1,sum))

LUAD_amp2<-data.frame(gene=rownames(LUAD_amp),LUAD_CNV=apply(LUAD_amp,1,sum))

LUAD_loss2<-data.frame(gene=rownames(LUAD_loss),LUAD_CNV=apply(LUAD_loss,1, sum))


# 2 BLCA -----------------------------------------------------------------------

BLCA<-fread("f:/workplace/BLCA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()

rownames(BLCA)<-BLCA$`Gene Symbol`
BLCA<-BLCA[,-1]

both_gene<-intersect(path_immu_gene,rownames(BLCA))

BLCA_CNV<-BLCA[both_gene,]

BLCA_amp<-ifelse(BLCA_CNV>0,1,0)
BLCA_loss<-ifelse(BLCA_CNV<0,1,0)

BLCA_high_amp<-ifelse(BLCA_CNV>1,1,0)
BLCA_high_amp2<-data.frame(gene=rownames(BLCA_high_amp),BLCA_high_CNV=apply(BLCA_high_amp,1,sum))

BLCA_amp2<-data.frame(gene=rownames(BLCA_amp),BLCA_CNV=apply(BLCA_amp,1,sum))

BLCA_loss2<-data.frame(gene=rownames(BLCA_loss),BLCA_CNV=apply(BLCA_loss,1, sum))

l_amp<-list(LUAD_amp2,BLCA_amp2,HNSC_amp2,KIRC_amp2,KIRP_amp2,LIHC_amp2,LUSC_amp2,THCA_amp2,COAD_amp2,GBM_amp2,KICH_amp2,LGG_amp2,BRCA_amp2,READ_amp2,UCS_amp2,UCEC_amp2,OV_amp2,PRAD_amp2,STAD_amp2,SKCM_amp2,CESC_amp2,ACC_amp2,PCPG_amp2,SARC_amp2,LAML_amp2,PAAD_amp2,ESCA_amp2,TGCT_amp2,THYM_amp2,MESO_amp2,UVM_amp2,DLBC_amp2,CHOL_amp2)
l_loss<-list(LUAD_loss2,BLCA_loss2,HNSC_loss2,KIRC_loss2,KIRP_loss2,LIHC_loss2,LUSC_loss2,THCA_loss2,COAD_loss2,GBM_loss2,KICH_loss2,LGG_loss2,BRCA_loss2,READ_loss2,UCS_loss2,UCEC_loss2,OV_loss2,PRAD_loss2,STAD_loss2,SKCM_loss2,CESC_loss2,ACC_loss2,PCPG_loss2,SARC_loss2,LAML_loss2,PAAD_loss2,ESCA_loss2,TGCT_loss2,THYM_loss2,MESO_loss2,UVM_loss2,DLBC_loss2,CHOL_loss2)
l_high_amp<-list(LUAD_high_amp2,BLCA_high_amp2,HNSC_high_amp2,KIRC_high_amp2,KIRP_high_amp2,LIHC_high_amp2,LUSC_high_amp2,THCA_high_amp2,COAD_high_amp2,GBM_high_amp2,KICH_high_amp2,LGG_high_amp2,BRCA_high_amp2,READ_high_amp2,UCS_high_amp2,UCEC_high_amp2,OV_high_amp2,PRAD_high_amp2,STAD_high_amp2,SKCM_high_amp2,CESC_high_amp2,ACC_high_amp2,PCPG_high_amp2,SARC_high_amp2,LAML_high_amp2,PAAD_high_amp2,ESCA_high_amp2,TGCT_high_amp2,THYM_high_amp2,MESO_high_amp2,UVM_high_amp2,DLBC_high_amp2,CHOL_high_amp2)

a<-data.frame(gene=both_gene)
for(i in 1:length(l_amp)){
  a<-inner_join(a,l_amp[[i]])
}
rownames(a)<-a$gene
a2<-a[,-1]
cnv_amp<-t(a2)%>%as.data.frame()
amp_sum<-apply(cnv_amp,2,sum)
cnv_amp[34,]=amp_sum
rownames(cnv_amp)[34]<-"Allcancer"


s<-data.frame(gene=both_gene)
for(i in 1:length(l_loss)){
  s<-inner_join(s,l_loss[[i]])
}
rownames(s)<-s$gene
s2<-s[,-1]
cnv_loss<-t(s2)%>%as.data.frame()
loss_sum<-apply(cnv_loss,2,sum)
cnv_loss[34,]=loss_sum
rownames(cnv_loss)[34]<-"Allcancer"


g<-data.frame(gene=both_gene)
for(i in 1:length(l_high_amp)){
  g<-inner_join(g,l_high_amp[[i]])
}
rownames(g)<-g$gene
g2<-g[,-1]
cnv_high_amp<-t(g2)%>%as.data.frame()
high_amp_sum<-apply(cnv_high_amp,2,sum)
cnv_high_amp[34,]=high_amp_sum
rownames(cnv_high_amp)[34]<-"Allcancer"


cnv_sum<-fread("f:/workplace/结果2/结果2-CNV频数统计.csv")%>%as.data.frame()
cnv_sum[34,]<-c("Allcancer",sum(cnv_sum$cases))


rownames(cnv_amp)<-paste(cnv_sum$cancer,cnv_sum$cases,sep = "-")
rownames(cnv_loss)<-paste(cnv_sum$cancer,cnv_sum$cases,sep = "-")
rownames(cnv_high_amp)<-paste(cnv_sum$cancer,cnv_sum$cases,sep = "-")


# CNV geom_bar --------------------------------------------------------

mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
         "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
         "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
         "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3")

library(reshape2)

bar1<-melt(a,id.vars="gene",value.name = "value")

ggplot(bar1,aes(gene,weight=value,fill=variable))+
  geom_bar(position="stack")+
  theme(legend.position = "none",
        axis.text.x = element_text(size=10,angle = 90, hjust = 1),
        axis.title.x=element_blank())+
  scale_fill_manual(values = mycol)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

bar2<-melt(s,id.vars="gene",value.name = "value")

ggplot(bar2,aes(gene,weight=value,fill=variable))+
  geom_bar(position="stack")+
  theme(legend.position = "none",
        axis.text.x = element_text(size=10,angle = 90, hjust = 1),
        axis.title.x=element_blank())+
  scale_fill_manual(values = mycol)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))


# The frequency of CNV amplification and deletion for each gene in each cancer -------------------------------------------------

sum_num<-cnv_sum$cases%>%as.numeric()
amp_freq_sum<-c()
for(i in 1:dim(cnv_amp)[2]){
  
  amp_freq_sum[i]<-round((cnv_amp[34,i]/sum_num[34])*100,2)
}

cnv_amp[35,]<-paste(amp_freq_sum,"%",sep = "")
rownames(cnv_amp)[35]<-"Frequency"

write.csv(cnv_amp,"f:/workplace/结果2/泛癌-CNV扩增频率表.csv")

loss_freq_sum<-c()
for(i in 1:dim(cnv_loss)[2]){
  
  loss_freq_sum[i]<-round((cnv_loss[34,i]/sum_num[34])*100,2)
}

cnv_loss[35,]<-paste(loss_freq_sum,"%",sep = "")
rownames(cnv_loss)[35]<-"Frequency"

write.csv(cnv_loss,"f:/workplace/结果2/泛癌-CNV缺失频率表.csv")

# CNV frequency-heatmap---------------------------------------------------------

heat_amp<-cnv_amp[1:34,]

amp_freq=matrix(0,nrow = 34,ncol = 45)

for(i in 1:dim(heat_amp)[1]){
  for( j in 1:dim(heat_amp)[2]){
    amp_freq[i,j]<-round((as.numeric(as.character(heat_amp[i,j]))/sum_num[i])*100,2)
  }
}

rownames(amp_freq)<-rownames(heat_amp)
colnames(amp_freq)<-colnames(heat_amp)

amp_freq2<-ifelse(amp_freq<1,0,ifelse(amp_freq<2,1.5,ifelse(amp_freq<5,3,ifelse(amp_freq<10,6,ifelse(amp_freq<20,12,ifelse(amp_freq<30,22,35))))))

bk <- c(0,1,2,5,10,20,30,40)

pheatmap(amp_freq2,cellwidth =8, cellheight = 8, fontsize = 8,fontsize_row=8,
         method="pearson",
         scale="none", 
         cluster_rows=F,
         cluster_cols=T,
         color = colorRampPalette(colors = c("white","orange","black"))(7),
         show_colnames=T,show_rownames =T,
         #annotation_col = annotation_col,
         #annotation_colors = ann_colors,
         treeheight_col = "0",
         treeheight_row = "0",
         legend_breaks = c(0,1,2,5,10,20,30,40),
         breaks=bk,
         border_color = "NA")

amp_freq3<-matrix(0,nrow = 34,ncol = 45)

for(i in 1:dim(amp_freq)[1]){
  amp_freq3[i,]<-paste(amp_freq[i,],"%",sep = "")
}

rownames(amp_freq3)<-rownames(amp_freq)
colnames(amp_freq3)<-colnames(amp_freq)

write.csv(amp_freq3,"f:/workplace/结果2/热图CNV_扩增频率表.csv")


heat_loss<-cnv_loss[1:34,]

loss_freq=matrix(0,nrow = 34,ncol = 45)

for(i in 1:dim(heat_loss)[1]){
  for( j in 1:dim(heat_loss)[2]){
    loss_freq[i,j]<-round((as.numeric(as.character(heat_loss[i,j]))/sum_num[i])*100,2)
  }
}

rownames(loss_freq)<-rownames(heat_loss)
colnames(loss_freq)<-colnames(heat_loss)

loss_freq2<-ifelse(loss_freq<1,0,ifelse(loss_freq<2,1.5,ifelse(loss_freq<5,3,ifelse(loss_freq<10,6,ifelse(loss_freq<20,12,ifelse(loss_freq<30,22,35))))))

bk <- c(0,1,2,5,10,20,30,40)

pheatmap(loss_freq2,cellwidth =8, cellheight = 8, fontsize = 8,fontsize_row=8,
         method="pearson", 
         scale="none", 
         cluster_rows=F,
         cluster_cols=T,
         color = colorRampPalette(colors = c("white","orange","black"))(7),
         show_colnames=T,show_rownames =T,
         #annotation_col = annotation_col,
         #annotation_colors = ann_colors,
         treeheight_col = "0",
         treeheight_row = "0",
         legend_breaks = c(0,1,2,5,10,20,30,40),
         breaks=bk,
         border_color = "NA")

loss_freq3<-matrix(0,nrow = 34,ncol = 45)

for(i in 1:dim(loss_freq)[1]){
  loss_freq3[i,]<-paste(loss_freq[i,],"%",sep = "")
}

rownames(loss_freq3)<-rownames(loss_freq)
colnames(loss_freq3)<-colnames(loss_freq)

write.csv(loss_freq3,"f:/workplace/结果2/热图CNV_缺失频率表.csv")
