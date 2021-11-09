library(RCurl)
library(org.Hs.eg.db)
library(tidyverse)
gene90 <- read.csv("gene90.csv",check.names = T,header=T)
ids <- gene90$Gene.Symbol %>% unique()
data <- clusterProfiler::bitr(ids, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ko <- c()
for (g in data$ENTREZID) {
  query <- paste("https://www.kegg.jp/entry/hsa:", g, sep = "")
  result <- getURL(query)
  k <- result %>% str_extract("ko:K\\d+")
  ko <- c(ko, k)
}
data$KO <- ko
head(data)

