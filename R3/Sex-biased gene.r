library(dplyr)
library(ggplot2)
library(tibble)
library(data.table)
library(stringr)

gene90<-fread("f:/workplace/GeneSymbol_90.txt",header=F)%>%as.data.frame()%>%distinct()
gene90<-gene90$V1

HNSC_sexdiff_methy<-read.csv("F:/workplace/pathway multi data/HNSC_limma_pathway_methylation.csv")


logFC_t <- 0.01    
P.Value_t <- 0.05   

HNSC_sexdiff_methy$change <- ifelse(HNSC_sexdiff_methy$P.Value < P.Value_t & abs(HNSC_sexdiff_methy$logFC) > logFC_t,
                      ifelse(HNSC_sexdiff_methy$logFC > logFC_t ,'High','Low'),'NOT') 

table(HNSC_sexdiff_methy$change) 


# Volcano plot --------------------------------------------------------------

methy_gene<-intersect(HNSC_sexdiff_methy$X,gene90)

for_label <- HNSC_sexdiff_methy %>% dplyr::filter(HNSC_sexdiff_methy$X %in% methy_gene & HNSC_sexdiff_methy$change != "NOT")


p <- ggplot(data = HNSC_sexdiff_methy, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("red","blue","grey"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()



# heatmap ---------------------------------------------------------------


HNSC_methy_matrix<-fread("f:/workplace/HNSC/HNSC_multidata_sex_diff/HNSC_甲基化匹配后的样本.csv")%>%as.data.frame()
rownames(HNSC_methy_matrix)<-HNSC_methy_matrix[,1]
colnames(HNSC_methy_matrix)[1]<-"sample"

both_gene<-intersect(colnames(HNSC_methy_matrix),gene90)


HNSC_methy_matrix2<-HNSC_methy_matrix[,both_gene]


diff_methyl<-fread("f:/workplace/HNSC/HNSC_multidata_sex_diff/HNSC_diff_methyl.csv")%>%as.data.frame()


both_gene2<-intersect(colnames(HNSC_methy_matrix2),diff_methyl[,1])

HNSC_methy_matrix3<-HNSC_methy_matrix2[,both_gene2]

HNSC_heatinput<-t(HNSC_methy_matrix3)



library(pheatmap)

cormat<-round(cor(HNSC_heatinput , method = "pearson"),2)

heat_sample<-data.frame(sample=colnames(HNSC_heatinput))
identical(heat_sample$sample,HNSC_methy_matrix$sample)

heat_cov<-cbind(heat_sample,HNSC_methy_matrix[,c("sex","histology_subtype","age_at_diagnosis","stage","smoking","Grade")])  

heat_cov$age=ifelse(heat_cov$age_at_diagnosis>=50,">=50","<50")

annotation_col<-data.frame(Gender = factor(rep(c("Female","Male"),c(dim(HNSC_heatinput)[2]/2,dim(HNSC_heatinput)[2]/2))),
                           Age = as.factor(heat_cov$age),
                           Histology_subtype = as.factor(heat_cov$histology_subtype),
                           Stage = as.factor(heat_cov$stage),
                           Smoking = as.factor(heat_cov$smoking),
                           Grade = as.factor(heat_cov$Grade)) #按正常疾病分组
rownames(annotation_col) = colnames(HNSC_heatinput)

ann_colors = list(
  gender = c(Female = "#E83140", Male = "#3C7DAF")
)

n <- t(scale(t(HNSC_heatinput)))
n[n > 2] <- 2
n[n < -2] <- -2
df <- n
rownames(df)<-rownames(HNSC_heatinput)


library(RColorBrewer)  

pheatmap(HNSC_heatinput,cellwidth =1, cellheight = 4, fontsize = 6,fontsize_row=4,
         method="pearson", 
         scale="row", 
         cluster_rows=T,
         cluster_cols=F,
         color = colorRampPalette(c("#3C7DAF", "white", "#E83140"))(20),
         show_colnames=F,show_rownames =T,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         treeheight_col = "0",
         border_color = "NA")


# The pie chart----------------------------------------------------------

library("RColorBrewer")
display.brewer.all()

HNSC_me_pie<-c(11,16,20)
label<-c("Male-biased methy","Famale-biased methy","other methy")
piepercent<- round(100*HNSC_me_pie/sum(HNSC_me_pie), 2)
piepercent2<-paste(piepercent,rep("%",3),sep = "")

pie(HNSC_me_pie,labels = piepercent2,main = "number of sex-biased gene",col=c("#92C5DE","#FC8D59","#FFFFBF"))
legend("topright", c("Male-biased methy","Famale-biased methy","other methy"), cex = 0.6,fill = c("#92C5DE","#FC8D59","#FFFFBF"))
dev.off()



# circos plot--------------------------------------------------------

library(dplyr)
library(ggplot2)
library(tibble)
library(data.table)
library(stringr)

HNSC_methy_matrix<-fread("f:/workplace/HNSC/HNSC_multidata_sex_diff/HNSC_甲基化匹配后的样本.csv")%>%as.data.frame()

rownames(HNSC_methy_matrix)<-HNSC_methy_matrix[,1]
colnames(HNSC_methy_matrix)[1]<-"sample"

both_gene<-intersect(colnames(HNSC_methy_matrix),gene90)

HNSC_methy_matrix2<-HNSC_methy_matrix[,both_gene]

diff_methyl<-fread("f:/workplace/HNSC/HNSC_multidata_sex_diff/HNSC_diff_methyl.csv")%>%as.data.frame()

both_gene2<-intersect(colnames(HNSC_methy_matrix2),diff_methyl[,1])
#显著性别差异的通路基因

HNSC_methy_matrix3<-HNSC_methy_matrix2[,both_gene2]


HNSC_heatinput<-t(HNSC_methy_matrix3)


sexdiff_immu_gene<-rownames(HNSC_heatinput)


diff_methyl2<-subset(diff_methyl,diff_methyl$V1 %in% sexdiff_immu_gene)
diff_methyl2$group <- ifelse(diff_methyl2$logFC>0,"Female","Male")

HNSC_circos<-data.frame(gene=diff_methyl2$V1,group=diff_methyl2$group)


write.csv(HNSC_circos,"f:/workplace/HNSC/HNSC_multidata_sex_diff/HNSC_性别差异免疫衰老基因.csv",row.names=F)


genePos<-fread("genePos.txt")%>%as.data.frame()%>%distinct() 

diff_pos<-subset( genePos , genePos$Gene %in% diff_gene ) 
library(stringr)
geneLabel<-data.frame(Chromosome = str_c("chr",diff_pos$Chr,sep=""),
                      chromStart = diff_pos$Start,
                      chromEnd = diff_pos$End,
                      Gene = diff_pos$Gene)   

refer<-fread("refer.txt")%>%as.data.frame()%>%distinct() 

scatter<-fread("SNP6_nocnv_genomicSegment")%>%as.data.frame() 
scatter<-subset(scatter , scatter$sample %in% colnames(limma_input2)) 



  