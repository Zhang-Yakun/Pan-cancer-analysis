# Mutation data preprocessing--------------------------------------------------------------

setwd("f:/workplace/LUAD/")

library(dplyr)
library(data.table)

mutation<-fread("LUAD_mc3.txt")%>%as.data.frame()%>%distinct()

LUAD_mut<-subset(mutation,effect!="Silent") #After silencing mutations were removed, 16 mutation types remained

# Samples with excessive number of mutations were screened--------------------------------------------------------------

sample_mut_number<-table(LUAD_mut$sample) 
mut_number<-data.frame(sample=names(sample_mut_number),number=as.numeric(sample_mut_number))
mut_num<-subset(mut_number,number<1000) 

LUAD_mut2<-subset(LUAD_mut, sample %in% mut_num$sample ) 


# Highly mutated genes with a mutation frequency of 1% of non-silent mutations are retained ---------------------------------------------------

gene_mut_freq<-table(LUAD_mut2$sample,LUAD_mut2$gene)
gene_mut_freq<-as.matrix(gene_mut_freq)

zero<-function(x){
  ifelse(x>0,1,0)%>%sum()
}

a<-apply(gene_mut_freq, MARGIN=2, FUN = zero) 
n<-a[(a/dim(gene_mut_freq)[1])>=0.01]%>%names() 

LUAD_mut3<-subset(LUAD_mut2,gene%in%n)


mut_488<-gene_mut_freq[unique(LUAD_mut3$sample),n] 


write.csv(mut_488,"LUAD_mutation_preprocess_01.csv")

mutation_preprocess<-subset(mutation,mutation$sample%in%rownames(mut_488))
mutation_preprocess<-subset(mutation_preprocess,mutation_preprocess$gene%in%colnames(mut_488))

write.csv(mutation_preprocess,"LUAD_mutation_preprocess.csv",row.names = F)
