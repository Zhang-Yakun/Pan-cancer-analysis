
BiocManager::install("KEGGgraph")
library(KEGGgraph)


toyKGML <- system.file("extdata/hsa00020.xml", package="KEGGgraph")
toyGraph <- parseKGML2Graph(toyKGML, genesOnly=TRUE,expandGenes=TRUE)

mapkNodes<-nodes(toyGraph)
mapkEdges<-edges(toyGraph)


# Visualize the network   -----------------------------------------------------------

mapkNodes<-mapkNodes[sapply(mapkNodes,length)>0]

res<-lapply(1:length(mapkEdges),function(t){
  name<-names(mapkEdges)[t]
  len<-length(mapkEdges[[t]])
  do.call(rbind,lapply(1:len,function(n){
    c(name,mapkEdges[[t]][n])
  }))
  
})

result<-data.frame(do.call(rbind,res))
write.table(result,"KEGG_graph_edges.txt",sep = "\t",row.names = F,col.names = F,quote = F)
write.table(mapkNodes,"KEGG_graph_nodes.txt",sep = "\t",row.names = F,col.names = F,quote = F)


# Node degree -----------------------------------------------------------

mapkGoutdegree <- sapply(edges(toyGraph),length)
mapkGindegree <- sapply(inEdges(toyGraph),length)

degrees <- data.frame(indegrees = mapkGindegree,outdegrees = mapkGoutdegree)
head(degrees)



# Key genes in network ---------------------------------------------------------------

library(RBGL)


toyKGML <- system.file("extdata/hsa00020.xml", package="KEGGgraph")
toyGraph <- parseKGML2Graph(toyKGML, genesOnly=TRUE,expandGenes=TRUE)

bcc <- brandes.betweenness.centrality(toyGraph) #计算网络的中心度
rbccs <- bcc$relative.betweenness.centrality.vertices[1L,]
toprbccs <- sort(rbccs,decreasing = TRUE)[1:4]

toprbccs #hub
