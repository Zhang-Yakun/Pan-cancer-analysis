library(dplyr)
library(data.table)




tumor<-fread("THCA_clinicalMatrix",header = T)%>%as.data.frame() 
survival<-fread("THCA_survival.txt",header = T)%>%as.data.frame() 

tumor2<-data.frame(sample = tumor$sampleID,
                   sex = tumor$gender,
                   age_at_diagnosis = tumor$age_at_initial_pathologic_diagnosis,
                   histology_subtype = tumor$histological_type,
                   stage = tumor$pathologic_stage,
                   smoking = tumor$tobacco_smoking_history,
                   Grade = tumor$neoplasm_histologic_grade,
                   purity = tumor$ABSOLUTE_Purity
)

survival2<-survival[,c(1,3,4)]
colnames(survival2)<-c("sample","OS","os.time")

ts<-inner_join(tumor2,survival2)%>%distinct() 

write.csv(ts,"THCA_clinical_survival.csv",row.names = F)
