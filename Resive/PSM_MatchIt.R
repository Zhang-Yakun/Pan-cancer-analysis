install.packages("MatchIt")
library(MatchIt)

#倾向性得分的包MatchIt

setwd("f:/Rwork/")

data<-read.csv("data.csv")  #读入文件包含协变量信息
data$group2 <- as.logical(data$group == 'uncensored')  #把分组的组别变为逻辑变量TRUE和FALSE

data_match <- matchit(group2~sex+age+TNM, data = data, method="nearest", ratio=1, caliper = 0.02)
#matchit函数用来倾向性得分匹配

summary(data_match)  #展示匹配前后的各项数据

plot(data_match, type = "jitter") #绘制匹配前后倾向值评分的分布图

plot(data_match, type = "hist") #绘制匹配前后倾向值评分的直方图


dta_m <- match.data(data_match) #match.data() 输出从matchit()得到的匹配数据集，以表格形式展示


###安装软件包
###统计检验：采用'tableone'包
###vars是给定的特征向量，分类变量是factor类型的变量，连续型变量当作numeric变量，
#factorVar是factor类型的变量，strata为分组

install.packages("tableone")
library(tableone)

table1 <- CreateTableOne(vars = c('sex', 'age', 'TNM'),
                         data = dta_m,
                         factorVars = c('sex', 'TNM'),
                         strata = 'group')


###数据导出
write.csv(dta_m, file = "data.match.csv")
#结果文件中，distance是倾向性得分的分值，weights是权重。