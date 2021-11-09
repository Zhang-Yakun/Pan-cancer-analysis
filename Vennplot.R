library(VennDiagram)
library(colorfulVennPlot)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(purrr)
library(RColorBrewer)
library(grDevices)
library(Cairo)
library(stringr)
library(tibble)
library(tidyr)
library(stringr)
options(stringsAsFactors = FALSE) #禁止chr转成factor


# 被激活的通路交集维恩图 -------------------------------------------------------------

#####输入

kegg<-read.csv("f:/workplace/免疫衰老基因富集KEGG通路结果-180.csv")

#p53信号通路

k<-subset(kegg,kegg$ID=="hsa04115")
p53_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

#PI3K-AKT信号通路
k<-subset(kegg,kegg$ID=="hsa04151")
PI3K_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

# Apoptosis号通路
k<-subset(kegg,kegg$ID=="hsa04210")
apop_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

# antigen-通路
k<-subset(kegg,kegg$ID=="hsa04612")
antigen_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

# PD-L1 pathway
k<-subset(kegg,kegg$ID=="hsa05235")
PD1_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

ls<-list(p53_gene,PI3K_gene,apop_gene,antigen_gene)
names(ls)<-c("P53 pathway","PI3K-AKT pathway","Apoptosis","Antigen pathway")
str(ls)

venn.diagram(x = ls, filename = "f:/workplace/结果1/venn_test.png")

# 子区域数量
number_area <- 2^length(ls) - 1

# 子区域编号
## 自定义一个函数，将整数转换成二进制字符串
intToBin <- function(x){
  if (x == 1)
    1
  else if (x == 0)
    NULL
  else {
    mod <- x %% 2
    c(intToBin((x-mod) %/% 2), mod)
  }
}

x_area <- seq(number_area) %>% map(., ~intToBin(.x)) %>%  # 转换成二进制字符串
  map_chr(., ~paste0(.x, collapse = "")) %>% 
  map_chr(., ~str_pad(.x, width = length(ls), side = "left", pad = "0"))

# 取ls中的向量取并集
G_union <- ls$`P53 pathway` %>% union(ls$`PI3K-AKT pathway`) %>% 
  union(ls$Apoptosis) %>% union(ls$`Antigen pathway`)

# 自定义一个函数，计算子区域中元素的数量
area_calculate <- function(data_ls, character_area){
  character_num <- 1:4 %>% map_chr(., ~substr(character_area, .x, .x)) %>% 
    as.integer() %>% as.logical()
  
  element_alone <- G_union
  for (i in 1:4) {
    element_alone <- 
      if (character_num[i]) {
        intersect(element_alone, data_ls[[i]])
      } else {
        setdiff(element_alone, data_ls[[i]])
      }
  }
  return(element_alone)
}

# 调用自定义的函数，求各个子区域的元素
element_ls <- map(x_area, ~area_calculate(data_ls = ls, character_area = .x))

# 计算各个子区域中元素的数量：
quantity_area <- map_int(element_ls, length)

# 计算各个子区域中元素数量的百分比
percent_area <- (quantity_area / sum(quantity_area)) %>% round(3) # 四舍五入保留3位小数
percent_area <- (percent_area * 100) %>% paste0("%")

# 色板长度
length_pallete <- max(quantity_area) - min(quantity_area) + 1
color_area <- colorRampPalette(brewer.pal(n = 9, name = "Reds"))(length_pallete) 
#其中YlGn可以换成其他颜色组合，例如Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd
color_tb <- tibble(quantity = seq(min(quantity_area), max(quantity_area), by = 1),
                   color = color_area) 

nest1 <- tibble(quantity = quantity_area, percent = percent_area, area = x_area) %>% 
  group_by(quantity) %>% nest() %>% 
  left_join(color_tb, by = "quantity") %>% 
  arrange(quantity) %>% 
  unnest()

# Venn图中显示数量
regions <- nest1$quantity
names(regions) <- nest1$area

setwd("f:/workplace/结果1")
CairoPDF(file = "venn_num2.pdf", width = 10, height = 10)
plot.new()
plotVenn4d(regions, Colors = nest1$color, Title = "", 
           labels = c("P53 pathway","PI3K-AKT pathway","Apoptosis","Antigen pathway")) #从左至右与输入文件顺序对应
dev.off() # 关闭绘图以保存

# 被抑制的通路交集维恩图 -------------------------------------------------------------

#PI3K-AKT信号通路
k<-subset(kegg,kegg$ID=="hsa04151")
PI3K_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

#MAPK信号通路
k<-subset(kegg,kegg$ID=="hsa04010")
MAPK_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

#chemokine信号通路
k<-subset(kegg,kegg$ID=="hsa04062")
chemo_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

#JAK-STAT信号通路
k<-subset(kegg,kegg$ID=="hsa04630")
JAK_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

#longevity信号通路
k<-subset(kegg,kegg$ID=="hsa04211")
longe_gene<-str_split(k$geneID,pattern = "/")%>%unlist()

ls<-list(PI3K_gene,MAPK_gene,chemo_gene,JAK_gene)
names(ls)<-c("PI3K-AKT pathway","MAPK pathway","chemokine pathway","JAK-STAT pathway")
str(ls)

venn.diagram(x = ls, filename = "f:/workplace/结果1/venn_test2.png")


# 子区域数量
number_area <- 2^length(ls) - 1

# 子区域编号
## 自定义一个函数，将整数转换成二进制字符串
intToBin <- function(x){
  if (x == 1)
    1
  else if (x == 0)
    NULL
  else {
    mod <- x %% 2
    c(intToBin((x-mod) %/% 2), mod)
  }
}

x_area <- seq(number_area) %>% map(., ~intToBin(.x)) %>%  # 转换成二进制字符串
  map_chr(., ~paste0(.x, collapse = "")) %>% 
  map_chr(., ~str_pad(.x, width = length(ls), side = "left", pad = "0"))

# 取ls中的向量取并集
G_union <- ls$`PI3K-AKT pathway` %>% union(ls$`MAPK pathway`) %>% 
  union(ls$`chemokine pathway`) %>% union(ls$`JAK-STAT pathway`)

# 自定义一个函数，计算子区域中元素的数量
area_calculate <- function(data_ls, character_area){
  character_num <- 1:4 %>% map_chr(., ~substr(character_area, .x, .x)) %>% 
    as.integer() %>% as.logical()
  
  element_alone <- G_union
  for (i in 1:4) {
    element_alone <- 
      if (character_num[i]) {
        intersect(element_alone, data_ls[[i]])
      } else {
        setdiff(element_alone, data_ls[[i]])
      }
  }
  return(element_alone)
}

# 调用自定义的函数，求各个子区域的元素
element_ls <- map(x_area, ~area_calculate(data_ls = ls, character_area = .x))

# 计算各个子区域中元素的数量：
quantity_area <- map_int(element_ls, length)

# 计算各个子区域中元素数量的百分比
percent_area <- (quantity_area / sum(quantity_area)) %>% round(3) # 四舍五入保留3位小数
percent_area <- (percent_area * 100) %>% paste0("%")

# 色板长度
length_pallete <- max(quantity_area) - min(quantity_area) + 1
color_area <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(length_pallete) 
#其中YlGn可以换成其他颜色组合，例如Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd
color_tb <- tibble(quantity = seq(min(quantity_area), max(quantity_area), by = 1),
                   color = color_area) 

nest1 <- tibble(quantity = quantity_area, percent = percent_area, area = x_area) %>% 
  group_by(quantity) %>% nest() %>% 
  left_join(color_tb, by = "quantity") %>% 
  arrange(quantity) %>% 
  unnest()

# Venn图中显示数量
regions <- nest1$quantity
names(regions) <- nest1$area

setwd("f:/workplace/结果1")
CairoPDF(file = "venn_num3.pdf", width = 10, height = 10)
plot.new()
plotVenn4d(regions, Colors = nest1$color, Title = "", 
           labels = c("PI3K-AKT pathway","MAPK pathway","chemokine pathway","JAK-STAT pathway")) #从左至右与输入文件顺序对应
dev.off() # 关闭绘图以保存
