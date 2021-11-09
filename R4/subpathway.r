setwd("f:/workplace/subpathway/")

library(dplyr)
library(data.table)
library(tibble)
library(clusterProfiler) #转换基因ID的包
library(org.Hs.eg.db)
library(stringr)

# 通路基因 --------------------------------------------------------------------


hsa16<-read.csv("f:/workplace/LUAD/enrichKEGG16.csv",header=T) #16条免疫衰老相关通路
pathway<-hsa16$ID%>%as.character()  #16条通路的ID向量

hsa_gene<-read.csv("f:/workplace/LUAD/HSA_KEGG.csv",header=T)   #KEGG通路和基因EntrzID对应列表

p53_pathway<-subset(hsa_gene,hsa_gene$DESCRIPTION=="p53 signaling pathway")  #p53信号通路的基因
p53_pathway_gene<-bitr(geneID=p53_pathway$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(p53_pathway_gene)<-c("ENTREZID","gene") %>%unique() #p53信号通路的基因和ENTRZID


# 免疫衰老基因 ------------------------------------------------------------------


immune_gene90<-read.table("f:/workplace/LUAD/GeneSymbol_90.txt",header=F) #90个免疫衰老基因列表
immune_gene90<-immune_gene90$V1%>%unique()

#把90个免疫衰老基因的symbol换成EntrzID

keytypes(org.Hs.eg.db)
gene90<-bitr(geneID=immune_gene90, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(gene90)<-c("gene","ENTREZID") %>%unique()


# LUAD表达数据来挖掘p53的子通路 ------------------------------------------------------


LUAD_exp<-fread( "f:/workplace/LUAD/LUAD_readcount.genes.tpm.txt")%>%as.data.frame() #LUAD的TPM基因表达谱
rownames(LUAD_exp)<-LUAD_exp$V1
LUAD_exp<-LUAD_exp[,-1]

#log2+1处理TPM表达谱
LUAD_exp<-log2(LUAD_exp+1) # log2(TPM+1)表达数据【17557*576】

p53_pathway_gene<-bitr(geneID=p53_pathway$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(p53_pathway_gene)<-c("ENTREZID","gene") %>%unique()

p53_pathway_exp<-LUAD_exp[p53_pathway_gene$gene,]


# p53通路中基因在癌症中的表达差异 -----------------------------------------------

# 用Limma包找癌症LUAD中差异表达的基因【Tumor VS normal】 ----------------------------------------


limma_input<-as.data.frame(p53_pathway_exp)  #行是基因，列是样本）的表达谱，之后作为limma输入做差异

sig<-substr(colnames(limma_input),14,14)%>%as.numeric() # "0"是癌症样本[407]，"1"是正常样本[19]

limma_input2<-limma_input[,order(sig)] #把limma的输入表达谱按 癌症，正常 排序【0,1】


group_list=c(rep('Tumor',table(sig)[1]%>%as.numeric()),rep('Normal',table(sig)[2]%>%as.numeric()))
## 强制限定顺序
group_list <- factor(group_list,levels = c("Tumor","Normal")) 

library(limma)
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
fit <- lmFit(limma_input2 , design)
contrast.matrix <- makeContrasts(Tumor - Normal,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
p53_DEG <- na.omit(tempOutput)  #limma结果——p53通路中癌症差异表达基因的P值和logFC值
p53_DEG$gene<-rownames(p53_DEG)

# 通路中免疫衰老基因与非免疫衰老基因的差异-wilcox.test ----------------------------------------

pathway_immune_gene<-p53_DEG[intersect(gene90$gene,p53_pathway_gene$gene),]%>%distinct()%>%na.omit() 
#通路中的免疫衰老基因的logFC值

pathway_other_gene<-p53_DEG[setdiff(rownames(p53_DEG),gene90$gene),]%>%distinct()%>%na.omit()
#通路中非免疫衰老基因的logFC值

subpathway_wt<-wilcox.test(pathway_immune_gene$logFC,pathway_other_gene$logFC)
subpathway_pvalue<-data.frame(cancer="LUAD",wilcox.test.pvalue=subpathway_wt$p.value) 
#p53通路中免疫衰老子通路显著在癌症中差异
write.csv(subpathway_pvalue,"subpathway_p53_pvalue.csv",row.names = F)

pathway_immune_uniprotid<-bitr(geneID=pathway_immune_gene$gene, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db", drop = TRUE)
pathway_immune_uniprotid$UNIPROT

write.csv(pathway_immune_uniprotid$UNIPROT,"subpathway_p53_immunegene.csv",row.names = F)
#p53通路中免疫衰老基因的UNIPROTID输出文件,用来在KEGG中注释通路基因