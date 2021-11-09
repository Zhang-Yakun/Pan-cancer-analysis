

library(dplyr)
library(data.table)
library(tibble)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(stringr)
library(limma)

# immunosenescence pathway genes--------------------------------------------------------------------


hsa16<-read.csv("f:/workplace/LUAD/enrichKEGG16.csv",header=T) 
pathway<-hsa16$ID%>%as.character()  

hsa_gene<-read.csv("f:/workplace/LUAD/HSA_KEGG.csv",header=T) # EntrzID

p53_pathway<-subset(hsa_gene,hsa_gene$DESCRIPTION=="p53 signaling pathway")  
p53_pathway_gene<-bitr(geneID=p53_pathway$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(p53_pathway_gene)<-c("ENTREZID","gene") %>%unique()

immune_gene90<-read.table("f:/workplace/LUAD/GeneSymbol_90.txt",header=F) 
immune_gene90<-immune_gene90$V1%>%unique()

keytypes(org.Hs.eg.db)
gene90<-bitr(geneID=immune_gene90, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(gene90)<-c("gene","ENTREZID") %>%unique()


# Reconstruction subpathway ------------------------------------------------------


LUAD_exp<-fread( "f:/workplace/LUAD/LUAD_readcount.genes.tpm.txt")%>%as.data.frame() #LUAD的TPM基因表达谱
rownames(LUAD_exp)<-LUAD_exp$V1
LUAD_exp<-LUAD_exp[,-1]

LUAD_exp<-log2(LUAD_exp+1) 

p53_pathway_gene<-bitr(geneID=p53_pathway$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(p53_pathway_gene)<-c("ENTREZID","gene") %>%unique()

p53_pathway_exp<-LUAD_exp[p53_pathway_gene$gene,]



limma_input<-as.data.frame(p53_pathway_exp)  

sig<-substr(colnames(limma_input),14,14)%>%as.numeric() 

limma_input2<-limma_input[,order(sig)]


group_list=c(rep('Tumor',table(sig)[1]%>%as.numeric()),rep('Normal',table(sig)[2]%>%as.numeric()))

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
p53_DEG <- na.omit(tempOutput)  
p53_DEG$gene<-rownames(p53_DEG)

# Differences between immunosenescence genes and non-immunosenescence genes in the pathway-wilcox.test ----------------------------------------

pathway_immune_gene<-p53_DEG[intersect(gene90$gene,p53_pathway_gene$gene),]%>%distinct()%>%na.omit() 

pathway_other_gene<-p53_DEG[setdiff(rownames(p53_DEG),gene90$gene),]%>%distinct()%>%na.omit()

subpathway_wt<-wilcox.test(pathway_immune_gene$logFC,pathway_other_gene$logFC)
subpathway_pvalue<-data.frame(cancer="LUAD",wilcox.test.pvalue=subpathway_wt$p.value) 

write.csv(subpathway_pvalue,"subpathway_p53_pvalue.csv",row.names = F)

pathway_immune_uniprotid<-bitr(geneID=pathway_immune_gene$gene, fromType="SYMBOL", toType="UNIPROT", OrgDb="org.Hs.eg.db", drop = TRUE)
pathway_immune_uniprotid$UNIPROT

write.csv(pathway_immune_uniprotid$UNIPROT,"subpathway_p53_immunegene.csv",row.names = F)
