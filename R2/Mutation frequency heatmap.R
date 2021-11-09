
library(stringr)
library(ggplot2)
library("RColorBrewer")
library(data.table)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)


gene47<-read.table("pathway_gene_47.txt")
path_immu_gene<-gene47$V1%>%as.character()


# Mutation data preprocessing for 33 cancers -----------------------------------------------------------

# 1 LUAD ------------------------------------------------------------------

LUAD<-fread("f:/workplace/LUAD/LUAD_mc3.txt")%>%as.data.frame()%>%distinct()
LUAD_mut<-subset(LUAD,effect!="Silent")

LUAD_mut$sample%>%unique()%>%length()

LUAD_mut_47<-subset(LUAD_mut,LUAD_mut$gene %in% path_immu_gene)

t1<-table(LUAD_mut_47$gene)

LUAD_muts<-data.frame(gene=names(t1),LUAD_mut_number=as.numeric(t1))

# 2 BLCA -----------------------------------------------------------------------

BLCA<-fread("f:/workplace/BLCA/BLCA_mc3.txt")%>%as.data.frame()%>%distinct()
BLCA_mut<-subset(BLCA,effect!="Silent") 

BLCA_mut$sample%>%unique()%>%length()


BLCA_mut_47<-subset(BLCA_mut,BLCA_mut$gene %in% path_immu_gene)


t2<-table(BLCA_mut_47$gene)

BLCA_muts<-data.frame(gene=names(t2),BLCA_mut_number=as.numeric(t2))


cancer33<-c("LUAD_muts","BLCA_muts","HNSC_muts","KIRC_muts","KIRP_muts","LIHC_muts","LUSC_muts","THCA_muts","COAD_muts","GBM_muts","KICH_muts","LGG_muts","BRCA_muts","READ_muts","UCS_muts","UCEC_muts","OV_muts","PRAD_muts","STAD_muts","SKCM_muts","CESC_muts","ACC_muts","PCPG_muts","SARC_muts","LAML_muts","PAAD_muts","ESCA_muts","TGCT_muts","THYM_muts","MESO_muts","UVM_muts","DLBC_muts","CHOL_muts")

l<-list(LUAD_muts,BLCA_muts,HNSC_muts,KIRC_muts,KIRP_muts,LIHC_muts,LUSC_muts,THCA_muts,COAD_muts,GBM_muts,KICH_muts,LGG_muts,BRCA_muts,READ_muts,UCS_muts,UCEC_muts,OV_muts,PRAD_muts,STAD_muts,SKCM_muts,CESC_muts,ACC_muts,PCPG_muts,SARC_muts,LAML_muts,PAAD_muts,ESCA_muts,TGCT_muts,THYM_muts,MESO_muts,UVM_muts,DLBC_muts,CHOL_muts)

d<-data.frame(gene=path_immu_gene)
for(i in 1:length(l)){
  
  if(dim(l[[i]])[1] < 47){
    r1<-setdiff(path_immu_gene,as.character(l[[i]]$gene))
    r2<-rep(0,length(r1))
    dif_g<-data.frame(r1,r2)
    colnames(dif_g)<-colnames(l[[i]])
    l[[i]]<-rbind(l[[i]],dif_g)
  }
  d<-inner_join(d,l[[i]])
}

rownames(d)<-d$gene

# gene mutation frequency geom_bar------------------------------------------------------------

mycol<-c(brewer.pal(12,"Paired"),brewer.pal(8,"Accent"),brewer.pal(9,"Set1"),"black",brewer.pal(3,"Dark2"))

mycol<-c("#A6CEE3","#1F78B4","#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00","#CAB2D6",
         "#6A3D9A","#FFFF99","#B15928","#7FC97F","#BEAED4","#FDC086","#F0027F","#386CB0","#FFFF99",
         "#BF5B17","#666666","#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628",
         "#F781BF","black","red", "#1B9E77","#D95F02", "#7570B3")

library(reshape2)
d2<-d
d2[9,]<-d2[8,]
d2$gene[9]<-"TP53"
df2<-melt(d,id.vars="gene",value.name = "value")

ggplot(df2,aes(reorder(gene,value),weight=value,fill=variable))+
  geom_bar(position="stack")+
  theme(legend.position = "none",
        axis.text.x = element_text(size=10,angle = 90, hjust = 1),
        axis.title.x=element_blank())+
  scale_fill_manual(values = mycol)+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))


# Statistical analysis of mutation frequency -----------------------------------------------------------------
d<-d[,-1]
gene47_mut_num<-t(d)%>%as.data.frame()
rownames(gene47_mut_num) <- str_split(rownames(gene47_mut_num),"_",simplify=TRUE)[,1]

write.csv(gene47_mut_num,"f:/workplace/结果2/泛癌-47基因突变频数.csv")

pan_mut_num<-read.csv("f:/workplace/结果2/泛癌-47基因突变频数.csv")
pan_mut_num<-pan_mut_num[1:33,]
rownames(pan_mut_num)<-pan_mut_num$cancer
pan_mut_num<-pan_mut_num[,-1]
a<-apply(pan_mut_num,2,sum)
pan_mut_num2<-rbind(pan_mut_num,a)
rownames(pan_mut_num2)[34]<-"Sum"

for(i in 1:dim(pan_mut_num2)[2]){
  
  mut_sum[i]<-round((pan_mut_num2[34,i]/pan_mut_num2[34,1])*100,2)
}

pan_mut_num2[35,]<-paste(mut_sum,"%",sep = "")
rownames(pan_mut_num2)[35]<-"Frequency"

write.csv(pan_mut_num2,"f:/workplace/结果2/泛癌-47基因突变频率.csv")

# Heatmap of gene mutation frequency----------------------------------------------------------------

pan_mut_num2<-read.csv("f:/workplace/结果2/泛癌-47基因突变频率.csv")
rownames(pan_mut_num2)<-pan_mut_num2$X
pan_mut_num2<-pan_mut_num2[,-1]

heat_mut<-pan_mut_num2[1:34,2:48]
rownames(heat_mut)<-paste(rownames(pan_mut_num2),pan_mut_num2[,1],sep = "-")[-35]
cancer_sum<-as.character(pan_mut_num2$cases[-35])%>%as.numeric()

h=matrix(0,nrow = 34,ncol = 47)

for(i in 1:dim(heat_mut)[1]){
  for( j in 1:dim(heat_mut)[2]){
    h[i,j]<-round((as.numeric(as.character(heat_mut[i,j]))/cancer_sum[i])*100,2)
  }
}

rownames(h)<-rownames(heat_mut)
colnames(h)<-colnames(heat_mut)

h2<-ifelse(h<1,0,ifelse(h<2,1.5,ifelse(h<5,3,ifelse(h<10,6,ifelse(h<20,12,ifelse(h<30,22,35))))))

bk <- c(0,1,2,5,10,20,30,40)

pheatmap(h2,cellwidth =8, cellheight = 8, fontsize = 8,fontsize_row=8,
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


h3<-matrix(0,nrow = 34,ncol = 47)

for(i in 1:dim(h)[1]){
  h3[i,]<-paste(h[i,],"%",sep = "")
}

rownames(h3)<-rownames(h)
colnames(h3)<-colnames(h)

write.csv(h3,"f:/workplace/结果2/热图34，47-突变频率表.csv")
