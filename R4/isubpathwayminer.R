library(iSubpathwayMiner)
library(clusterProfiler) 
library(org.Hs.eg.db)
library(data.table)
library(tibble)
library(dplyr)

gene126<-read.table("genclip.txt",header = T)
gene2<-gene126$GeneSymbol
gene90<-read.table("GeneSymbol_90.txt")%>%unique
gene1<-gene90$V1

intersect(gene1,gene2)%>%length() 
gene176<-union(gene1,gene2)


keytypes(org.Hs.eg.db)
all_gene<-bitr(geneID=gene176, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(all_gene)<-c("gene","EntrzID") %>%unique()

moleculeList<-all_gene$EntrzID 

reGM<-SubpathwayGM(moleculeList,n=5,s=5) #subpathways
result<-printGraph(reGM$ann) 
result[1:10,]
#result


result1<-printGraph(reGM$ann,detail=TRUE)
write.csv(result1,file="subpathway_result.csv",row.names=FALSE)


plotAnnGraph("path:00230_1",reGM$subGraphList,reGM$ann,displayInR=TRUE,gotoKEGG=TRUE)
#Visualized subpathways at R and KEGG sites


