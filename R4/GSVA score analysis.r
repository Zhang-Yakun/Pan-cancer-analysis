

library(data.table)
library(dplyr)
library(stringr)
library(limma)
library(ggplot2)

# GSVA-output
gsva<-fread("gsva_output.csv")%>%as.data.frame()
rownames(gsva)<-gsva$V1
gsva<-gsva[,-1]

#sample group
sample_group<-ifelse(as.numeric(str_sub(colnames(gsva),14,15))<10,"Tumor","Normal")


# Differential analysis of ISC-pathways----------------------------------------------------------------

group_list<-data.frame(sample=colnames(gsva),group=sample_group)


design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva)


contrast.matrix <- makeContrasts(Tumor-Normal, levels = design)

# Tumor vs. Normal
fit <- lmFit(gsva, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)

#limma result
write.csv(x, "gsva_pancancer_limma.csv", quote = F)


df <- data.frame(ID = rownames(x), score = x$t)


# GSVA plot ------------------------------------------------------------------

#Group by GSVA score
cutoff <- 2
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

#Sort by GSVA score
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('palegreen3', 'dodgerblue4'), guide = FALSE) + 
  

  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, 
             size = 0.3) + 
  scale_y_continuous(limits = c(-22,15),breaks = c(-20,-15,-10,-5,0, 5, 10, 15))+
  
  #label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),
            size = 5, 
            hjust = "outward" ) +  
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=ID, color = group),
            size = 5, hjust = "inward") +  
  scale_colour_manual(values = c("black","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score, tumor \n versus non-malignant")+
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 0.6)) + 
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) #去除y轴
