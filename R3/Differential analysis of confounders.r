setwd("f:/workplace/pathway multi data/")


LUAD<-data.frame(Cancer=rep("LUAD",6),
                 Clinical=c("Age","Histological\nType","Stage","Purity","Smoking","Grade"),
                 Pvalue=c(0.9,0.19,0.04,0.03,0.66,NA))
LUAD$group<-ifelse(LUAD$Pvalue<0.05,"yellow","white")
BLCA<-data.frame(Cancer=rep("BLCA",6),
                 Clinical=c("Age","Histological\nType","Stage","Purity","Smoking","Grade"),
                 Pvalue=c(0.01,0.9,0.8,NA,0.02,NA))
BLCA$group<-ifelse(BLCA$Pvalue<0.05,"yellow","white")
HNSC<-data.frame(Cancer=rep("HNSC",6),
                 Clinical=c("Age","Histological\nType","Stage","Purity","Smoking","Grade"),
                 Pvalue=c(0.0005,0.26,0.84,NA,0.000003,0.01))
HNSC$group<-ifelse(HNSC$Pvalue<0.05,"yellow","white")
KIRC<-data.frame(Cancer=rep("KIRC",6),
                 Clinical=c("Age","Histological\nType","Stage","Purity","Smoking","Grade"),
                 Pvalue=c(0.006,NA,0.13,NA,0.05,0.40))
KIRC$group<-ifelse(KIRC$Pvalue<0.05,"yellow","white")
KIRP<-data.frame(Cancer=rep("KIRP",6),
                 Clinical=c("Age","Histological\nType","Stage","Purity","Smoking","Grade"),
                 Pvalue=c(0.22,NA,0.86,NA,0.37,NA))
KIRP$group<-ifelse(KIRP$Pvalue<0.05,"yellow","white")
LIHC<-data.frame(Cancer=rep("LIHC",6),
                 Clinical=c("Age","Histological\nType","Stage","Purity","Smoking","Grade"),
                 Pvalue=c(0.02,0.57,0.04,NA,NA,0.49))
LIHC$group<-ifelse(LIHC$Pvalue<0.05,"yellow","white")
LUSC<-data.frame(Cancer=rep("LUSC",6),
                 Clinical=c("Age","Histological\nType","Stage","Purity","Smoking","Grade"),
                 Pvalue=c(0.86,0.41,0.05,NA,0.22,NA))
LUSC$group<-ifelse(LUSC$Pvalue<0.05,"yellow","white")
THCA<-data.frame(Cancer=rep("THCA",6),
                 Clinical=c("Age","Histological\nType","Stage","Purity","Smoking","Grade"),
                 Pvalue=c(0.23,0.12,0.03,NA,NA,NA))
THCA$group<-ifelse(THCA$Pvalue<0.05,"yellow","white")

df<-rbind(LUAD,BLCA,HNSC,KIRC,KIRP,LIHC,LUSC,THCA)

ggplot(df, aes(Clinical, Cancer)) + 
  geom_tile(aes(fill = group), colour = "gray90",size=1)+
  scale_fill_manual(values = c("white","yellow"),na.value = "gray90")+
  theme_minimal()+# 不要背景
  geom_text(aes(label=Pvalue),size=5)+
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(size = 12),# 调整x轴文字
        legend.position = "none",
        axis.text.y = element_text(size = 12)) #调整y轴文字

