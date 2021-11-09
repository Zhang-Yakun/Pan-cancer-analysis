library(ggplot2)

cancer <- c("BLCA","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","THCA")

molecule1 <- rep("Mutation",8)
number1 <- c(NA,2,NA,2,1,4,NA,1)
d1 <- data.frame(CancerType=cancer, Molecules=molecule1, Numbers=number1 )

molecule2 <- rep("SCNA",8)
number2 <- c(10,33,15,15,16,27,27,10)
d2 <- data.frame(CancerType=cancer, Molecules=molecule2, Numbers=number2 )

molecule3 <- rep("Methylation",8)
number3 <- c(213,633,232,63,47,2,193,127)
d3 <- data.frame(CancerType=cancer, Molecules=molecule3, Numbers=number3 )

molecule4 <- rep("mRNA",8)
number4 <- c(119,445,165,400,260,307,281,44)
d4 <- data.frame(CancerType=cancer, Molecules=molecule4, Numbers=number4 )

molecule5 <- rep("Protein",8)
number5 <- c(4,NA,1,NA,NA,NA,1,NA)
d5<-data.frame(CancerType=cancer, Molecules=molecule5, Numbers=number5 )

sex_point<-rbind(d1,d2,d3,d4,d5)

gg<-ggplot(data=sex_point,aes(Molecules,CancerType )) + 
  geom_tile(fill="gray90", colour = "white",size=0.5)+
  geom_point(aes(color = Molecules,size=Numbers))+
  scale_size_continuous(range = c(0, 14))+
  theme_minimal()+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size = 14,color=c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3")),# 调整x轴文字
        axis.text.y = element_text(size = 14, color = "#893737"))
gg

yhist<-data.frame(cancer=c("THCA","LUSC","LUAD","LIHC","KIRP","KIRC","HNSC","BLCA"),
                  his<-c(182,502,340,324,480,413,1113,346))

yplot = ggplot(data = yhist, aes(x =cancer,y=his)) + 
  geom_bar(stat="identity",position = "dodge",width = 0.6,fill="#893737")+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()
        )+
  ylab("total sexual-dimorphism features")+
  coord_flip()+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line.x = element_line(colour = "black"))

ggpubr::ggarrange(gg,yplot,widths = c(5,3),heights = c(1,3), align = "hv")


