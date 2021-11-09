
library(stringr)
library(dplyr)
library(data.table)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(ggplot2)
library(limma)

# Expression matrix count2TPM ------------------------------------------------------------------

setwd("f:/workplace/BLCA/")

expression <- fread("BLCA_HiSeqV2")%>%as.data.frame()%>%distinct() #log2(count+1)
rownames(expression)<-expression$sample

expMatrix<-round((2^expression[,-1])-1)

rownames(expMatrix)<-rownames(expression)

# Length of genes ------------------------------------------------------------------


eff_length2 <- read.csv("f:/workplace/eff_length_symbol.csv", row.names = 1, header = T)
eff_length2$gene_id <- rownames(eff_length2) 

# read count2TPM ----------------------------------------------------------


feature_ids <- rownames(expMatrix)


if (!all(feature_ids %in% rownames(eff_length2))) {
  tbl <- table(feature_ids %in% rownames(eff_length2))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]], tbl[[1]])
  warning(msg1)
}

if (!identical(feature_ids, rownames(eff_length2))) {
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}

expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2), ]
mm <- match(rownames(expMatrix), rownames(eff_length2))
eff_length2 <- eff_length2[mm, ]

if (identical(rownames(eff_length2), rownames(expMatrix))) {
  print("GTF and expression matix now have the same gene and gene in same order")
}

x <- expMatrix / eff_length2$eff_length
expMatrix_tpm <- t(t(x) / colSums(x)) * 1e6   #TPM

expMatrix_tpm[1:3, 1:5]
write.table(expMatrix_tpm, "BLCA_readcount.genes.tpm.txt", sep = "\t", quote = F, row.names = T)

#log2(TPM+1)
expMatrix_tpm<-(log2(expMatrix_tpm+1)) 

# Limma 【Tumor VS normal】 ----------------------------------------


limma_input<-as.data.frame(expMatrix_tpm)  

sig<-substr(colnames(limma_input),14,14)%>%as.numeric() 
limma_input2<-limma_input[,order(sig)] 


group_list=c(rep('Tumor',table(sig)[1]%>%as.numeric()),rep('Normal',table(sig)[2]%>%as.numeric()))

group_list <- factor(group_list,levels = c("Tumor","Normal")) 


design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
fit <- lmFit(limma_input2 , design)
contrast.matrix <- makeContrasts(Tumor - Normal,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
DEG <- na.omit(tempOutput)  #limma results

write.csv(DEG,"BLCA_limma_exp.csv",row.names = T)

diff_gene<-subset(DEG,DEG$P.Value<0.05)%>%rownames  

DEG$gene<-rownames(DEG)


keytypes(org.Hs.eg.db)
DEG_Entrz<-bitr(geneID=DEG$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(DEG_Entrz)<-c("gene","EntrzID") %>%unique()
DEG2<-left_join(DEG,DEG_Entrz) #geneID EntrzID

# Hypergeometric test - Differentially expressed upregulation genes and pathway genes in cancer -------------------------------------------------

upper_gene<-dplyr::filter( DEG2 , DEG2$logFC > 1 & DEG2$P.Value < 0.05) 

hsa16<-read.csv("f:/workplace/enrichKEGG16.csv",header=T) 
pathway<-hsa16$ID%>%as.character()  

hsa_gene<-read.csv("f:/workplace/HSA_KEGG.csv",header=T)   

pathway_up_score=c()
p1=c()

for(i in 1:length(pathway)){
  
  x<-hsa_gene[hsa_gene$KEGGID==pathway[i],2]
  n<-rep(pathway[i],length(x))
  d<-data.frame(ID=n,gene=x) 
  
  
  N=30000  #up_gene$EntrzID%>%length()
  q <- intersect(upper_gene$EntrzID,d$gene)%>%length() 
  m <- length(d$gene) 
  n <- N-m   
  k <- length(upper_gene$EntrzID) 
  
  p_value <- 1 - phyper( q, m, n, k ) 
  p1[i] <- p_value          
  pathway_up_score[i] <- (-log10(p_value)) 
  
}

names(pathway_up_score)=pathway

# Hypergeometric test - Differentially expressed downregulation genes and pathway genes in cancer -------------------------------------------------

down_gene<-dplyr::filter( DEG2 , DEG2$logFC < -1 & DEG2$P.Value < 0.05) 


pathway_down_score=c()
p2=c()

for(i in 1:length(pathway)){
  
  x<-hsa_gene[hsa_gene$KEGGID==pathway[i],2]
  s<-rep(pathway[i],length(x))
  d<-data.frame(ID=s,gene=x) 
  
  N=30000  #do_gene$EntrzID%>%length()
  q <- intersect(down_gene$EntrzID,d$gene)%>%length()
  m <- length(d$gene)
  n <- N-m
  k <- length(down_gene$EntrzID) 
  p_value <- 1 - phyper( q, m, n, k ) 
  p2[i] <- p_value     
  pathway_down_score[i] <- (-log10(p_value))
  
}

names(pathway_down_score)=pathway

pathway_score_exp<-data.frame(activate_score = pathway_up_score, 
                              suppression_score = pathway_down_score,
                              activate_pvalue = p1,
                              suppression_pvalue = p2)

rownames(pathway_score_exp)<-pathway   

write.csv(pathway_score_exp,"BLCA_pathwayScore_expression.csv")
pathway_score_exp<-read.csv("BLCA_pathwayScore_expression.csv")
rownames(pathway_score_exp)<-pathway_score_exp$X

# Statistical chart of pathway scores -----------------------------------------------------------------

pathway_score_exp$suppression_score=paste("-",pathway_score_exp$suppression_score,sep="")%>%as.numeric()

c0=data.frame()

for(i in 1:length(pathway_score_exp$activate_score)){
  c1=data.frame(X1=rownames(pathway_score_exp)[i],X2="activate_score",X3=pathway_score_exp$activate_score[i])
  c2=data.frame(X1=rownames(pathway_score_exp)[i],X2="suppression_score",X3=pathway_score_exp$suppression_score[i])
  cc=rbind(c1,c2)
  c0<-rbind(c0,cc)
}
s<-subset(hsa_gene,hsa_gene$KEGGID%in%c0$X1)
s2<-data.frame(s$KEGGID,s$DESCRIPTION)%>%distinct()
colnames(s2)<-c("pathwayID","description")
colnames(c0)<-c("pathwayID","A","S")

ps<-inner_join(s2,c0)
write.csv(ps,"BLCA_ps_expression.csv",row.names = F)

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
ggsave("BLCA_pathwayscore_expression.pdf",width=6)