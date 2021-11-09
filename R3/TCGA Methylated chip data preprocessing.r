
# TCGA Methylated chip data loading and collation --------------------------------------------------


library(data.table)
library(impute)
library(ChAMP)
library(stringr)
library(tibble)
library(dplyr)
options(stringsAsFactors = F)
#if(!dir.exists("raw_data"))dir.create("raw_data")
#if(!dir.exists("Rdata"))dir.create("Rdata")
#if(!dir.exists("figure"))dir.create("figure")


methy<-fread("HumanMethylation450")%>%as.data.frame() #492样本，485578行探针


methyl <- column_to_rownames(methy,"sample")


sig<-substr(colnames(methyl),14,14)%>%as.numeric() 
methyl<-methyl[,order(sig)] 

group_list<-ifelse(as.numeric(str_sub(colnames(methyl),14,15))<10,"Tumor","Normal") 

table(group_list)
pd <- data.frame( sampleID = colnames(methyl), group_list = group_list)


beta=as.matrix(methyl)

beta=impute.knn(beta)  

sum(is.na(beta))

beta=beta$data
beta=beta+0.00001   

myLoad=champ.filter(beta = beta ,pd = pd) 
dim(myLoad$beta)  
save(myLoad,file = './Rdata/step1_myLoad.Rdata')


# Differential analysis-DMP -----------------------------------------------------------------


group_list <- pd$group_list
myDMP <- champ.DMP(beta = beta,pheno=group_list) 

head(myDMP$Tumor_to_Normal)

df_DMP <- myDMP$Tumor_to_Normal  
df_DMP=df_DMP[df_DMP$gene!="",]  

save(df_DMP,file = "./Rdata/step3.df_DMP.Rdata")
write.csv(df_DMP,"LUAD_methy_DMP_limmaresult.csv")

df_DMP<-read.csv("LUAD_methy_DMP_limmaresult.csv")

logFC_t <- 0.2     
P.Value_t <- 0.0000005  

df_DMP$change <- ifelse(df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$logFC) > logFC_t,
                        ifelse(df_DMP$logFC > logFC_t ,'UP','DOWN'),'NOT') 

table(df_DMP$change) 

# Differential DMP volcanic map ----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(tibble)
dat  = rownames_to_column(df_DMP)
for_label <- dat%>% head(3)
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(adj.P.Val))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p

ggsave()


# Pathway score tables -----------------------------------------------------------------

pathway_score_methyl$suppression_score=paste("-",pathway_score_methyl$suppression_score,sep="")%>%as.numeric()

c0=data.frame()

for(i in 1:length(pathway_score_methyl$activate_score)){
  c1=data.frame(X1=rownames(pathway_score_methyl)[i],X2="activate_score",X3=pathway_score_methyl$activate_score[i])
  c2=data.frame(X1=rownames(pathway_score_methyl)[i],X2="suppression_score",X3=pathway_score_methyl$suppression_score[i])
  cc=rbind(c1,c2)
  c0<-rbind(c0,cc)
}


s<-subset(hsa_gene,hsa_gene$KEGGID%in%c0$X1)
s2<-data.frame(s$KEGGID,s$DESCRIPTION)%>%distinct()
colnames(s2)<-c("pathwayID","description")
colnames(c0)<-c("pathwayID","A","S")

ps<-inner_join(s2,c0)


ggplot(data=ps,aes(A,description)) + 
  geom_point(aes(color = S), size=5)+
  scale_color_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        axis.text.y = element_text(size = 8))+
  
  labs(color ="Pathway Score")
ggsave("pathwayscore_methylation.pdf",width=6)
