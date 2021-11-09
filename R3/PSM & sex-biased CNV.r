

library(data.table)
library(dplyr)
library(data.table)

cnv<-fread("f:/workplace/THCA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")%>%as.data.frame()%>%distinct()
#cnv[24776*517] -2,-1,0,1,2: 2 copy del, 1 copy del, no change, amplification, high-amplification

rownames(cnv)<-cnv$`Gene Symbol`
cnv<-cnv[,-1]

# copy number data ---------------------------------------------------------------

hsa16<-read.csv("f:/workplace/enrichKEGG16.csv",header=T) 
pathway<-hsa16$ID%>%as.character()  

hsa_gene<-read.csv("f:/workplace/HSA_KEGG.csv",header=T)  

library(clusterProfiler)
library(org.Hs.eg.db)

pathway16<-dplyr::filter(hsa_gene,hsa_gene$KEGGID==pathway)  
pathway_gene<-bitr(geneID=pathway16$ENTREZID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db", drop = TRUE)
colnames(pathway_gene)<-c("ENTREZID","gene") %>%unique() 

both_gene<-intersect(pathway_gene$gene,rownames(cnv)) 

pathway_cnv<-cnv[both_gene,]


# Clinical and survival data -------------------------------------------------------------

ts<-fread("f:/workplace/THCA/THCA_clinical_survival.csv")%>%as.data.frame()


cnv<-t(pathway_cnv)
cnv2<-as.data.frame(cnv)

rownames(cnv2)<-rownames(cnv)
colnames(cnv2)<-colnames(cnv)

cnv3<-cnv2%>%mutate(sample=rownames(cnv2))%>%dplyr::select(sample,everything())

clinical<-subset(ts,ts$sample%in%intersect(ts$sample,cnv3$sample)) 
cnv_clinical<-inner_join(cnv3,clinical)  


# Propensity score matching-PSM ------------------------------------------------------

library(dplyr)
library(ggplot2)
library(MatchIt)

cnv_clinical$histology_subtype<-as.factor(cnv_clinical$histology_subtype)
cnv_clinical$histology_subtype=ifelse(cnv_clinical$histology_subtype== "Thyroid Papillary Carcinoma - Classical/usual",1,cnv_clinical$histology_subtype)


cnv_clinical$stage<-as.factor(cnv_clinical$stage)
cnv_clinical$stage=ifelse(cnv_clinical$stage== "",1,cnv_clinical$stage)


cnv_clinical$Grade<-as.factor(cnv_clinical$Grade)
cnv_clinical$Grade=ifelse(cnv_clinical$Grade== "G1",1,cnv_clinical$Grade)

cnv_clinical$smoking<-as.factor(cnv_clinical$smoking)
cnv_clinical$smoking=ifelse(cnv_clinical$smoking == "[Discrepancy]",1,cnv_clinical$smoking)

cnv_clinical %>% group_by(sex) %>% summarise(sample_number = n(),
                                             age = mean(age_at_diagnosis),
                                             histology_subtype = median(histology_subtype),
                                             survival_time = mean(os.time),
                                             stage=median(stage)
                                             #purity = median(purity),
                                             #grade = median(Grade)
                                             #smoking = median(smoking)
)

cnv_clinical<-cnv_clinical %>% mutate(sex=ifelse(sex=="MALE",1,0)) 

mc_ps <- glm(sex ~ age_at_diagnosis  + histology_subtype +stage ,
             family = binomial(), data = cnv_clinical) 
summary(mc_ps)

pre_mc_ps<- data.frame(pr_score = predict(mc_ps, type = "response"),
                       sex = mc_ps$model$sex) 
head(pre_mc_ps)

library(ggplot2)
pre_mc_ps %>%
  ggplot(aes(x = pr_score)) +
  geom_histogram(color = "white") +
  facet_wrap(~sex) +
  xlab("Probability") +
  theme_bw() 


# The propensity score is matched with algorithm ------------------------------------------------------------

clinical_cov<-c("age_at_diagnosis","histology_subtype","stage") 
gene<-colnames(cnv3)[-1] 
cnv_clinical_nomiss <- cnv_clinical %>%  # MatchIt-NA
  dplyr::select(sample,gene, sex, one_of(clinical_cov)) %>%
  na.omit() 



rownames(cnv_clinical_nomiss)<-cnv_clinical_nomiss$sample
cnv_clinical_nomiss<-cnv_clinical_nomiss[,-1]  

mod_match <- matchit(sex ~ age_at_diagnosis + histology_subtype + stage,
                     method = "nearest", data = cnv_clinical_nomiss) 
dta_m <- match.data(mod_match) 
dim(dta_m) 

# Check the balance of covariables in the matching sample -----------------------------------------------------------

fn_bal <- function(dta, variable) {
  dta$variable <- dta[, variable]
  dta$sex <- as.factor(dta$sex)
  support <- c(min(dta$variable), max(dta$variable))
  ggplot(dta, aes(x = distance, y = variable, color = sex)) +
    geom_point(alpha = 0.2, size = 1.3) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(variable) +
    theme_bw() +
    ylim(support)
} 

library(gridExtra) 

grid.arrange(
  fn_bal(dta_m, "age_at_diagnosis") + theme(legend.position = "none"),
  fn_bal(dta_m, "histology_subtype"),
  #fn_bal(dta_m, "smoking"),
  fn_bal(dta_m, "stage"),
  nrow = 2, widths = c(1, 0.8)
)

dta_m %>%
  group_by(sex) %>%
  dplyr::select(one_of(clinical_cov)) %>%
  summarise_all(funs(mean))  #psm之后两组样本协变量的均值差异

lapply(clinical_cov, function(v) {
  t.test(dta_m[, v] ~ dta_m$sex)
}) 

# After the covariates are balanced, assess the differences in the expressed values between the two groups of samples ----------------------------------------

dta_m2<-dta_m[,1:length(both_gene)]

limma_input<-t(dta_m2)
limma_input<-as.data.frame(limma_input)  

sex<-dta_m$sex
limma_input2<-limma_input[,order(sex)]


group_list=c(rep('Female',(dim(limma_input)[2])/2),rep('Male',(dim(limma_input)[2])/2))

group_list <- factor(group_list,levels = c("Female","Male"),ordered = F) 

library(limma)
design <- model.matrix(~ 0 + group_list)
colnames(design) <- levels(group_list)
fit <- lmFit(limma_input2 , design)
contrast.matrix <- makeContrasts(Female - Male,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG86 <- na.omit(tempOutput)  

write.csv(nrDEG86,"f:/workplace/THCA/THCA_multidata_sex_diff/THCA_limma_cnv.csv",row.names = T)

diff_gene<-subset(nrDEG86,nrDEG86$P.Value<0.05)%>%rownames()

diff<-subset(nrDEG86,nrDEG86$P.Value<0.05)
diff<-data.frame(Gene=rownames(diff),Ratio1=diff$logFC,Ratio2=diff$AveExpr,Stat=diff$t,Pvalue=diff$P.Value,adjPvalue=diff$adj.P.Val)

write.csv(diff,"f:/workplace/THCA/THCA_multidata_sex_diff/THCA_diff_CNV.csv",row.names = F)

