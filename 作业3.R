####作业3####
####LUAD ESTIMATE####
#计算患者免疫评分与肿瘤纯度#
setwd("TCGA-LUAD")
setwd("ESTIMATE")  #设置工作目录
#安装包
library(utils) #这个包应该不用下载，自带的
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

#计算免疫评分
filterCommonGenes(input.f = "tpms01A_log2.txt",   #输入文件名
                  output.f = "tpms01A_log2.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("tpms01A_log2.gct",   #刚才的输出文件名
              "tpms01A_log2_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#提取结果
ESTIMATE_result <- read.table("tpms01A_log2_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ESTIMATE_result <- ESTIMATE_result[,-1]   
colnames(ESTIMATE_result) <- ESTIMATE_result[1,]   
ESTIMATE_result <- as.data.frame(t(ESTIMATE_result[-1,]))
rownames(ESTIMATE_result) <- colnames(exp)
write.table(ESTIMATE_result, file = "ESTIMATE_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F) 

#LUAD 计算肿瘤样本cibersort 画彩虹图 根据BTK分组画图
####LUAD cibersort####
setwd("TCGA-LUAD")
setwd("CIBERSORT")   
#install.packages('e1071')
#install.packages('parallel')
#install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.17")
library(e1071)
library(parallel)
library(preprocessCore)
library(tidyverse)
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'tpms01A_log2.txt'   #肿瘤患者表达谱
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)

res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
ciber.res <- as.data.frame(ciber.res)
write.table(ciber.res,"ciber.res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####cibersort彩虹图 无需掌握####
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板

####BTK 分组比较图####
#读取肿瘤患者01A表达谱
exp <- read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
med=median(as.numeric(exp["BTK",]))
exp <- exp %>% t() %>% as.data.frame()
exp <- exp %>% mutate(group=factor(ifelse(exp$BTK>med,"high","low"),levels = c("low","high")))
class(exp$group)
a <- ciber.res
identical(rownames(a),rownames(exp))
a$group <- exp$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=CIBERSORT,value = Fraction,-c(group,sample))
ggboxplot(b, x = "CIBERSORT", y = "Fraction",
          fill = "group", palette = "jco")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

####LUAD cibersort 自己和自己做相关性 绘图####
setwd("TCGA-LUAD")
setwd("COR")
#install.packages("corrplot")
library(corrplot)
library(ggcorrplot)
library(tidyverse)
ciber <- read.table("ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

cor <-sapply(ciber,function(x,y) cor(x,y,method="spearman"),ciber)
rownames(cor)<-colnames(ciber)

corrplot(cor,
         method="pie",
         col=colorRampPalette(c("#01468b","white","#ee0000"))(100),
         type="upper",
         addCoef.col = "black",
         number.cex = 0.35,
         tl.col="black", tl.srt=45,
         tl.cex = 0.7,
         diag=FALSE)

#LUAD BTK 和 cibersort里22个基因做散点图 画22张图
####BTK与cibersort相关性散点图####
setwd("TCGA-LUAD")
setwd("sandian")
exp = read.table("tpms01A_log2.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber = read.table("ciber.res.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
ciber <- ciber %>% t() %>% as.data.frame()
rownames(ciber) <- gsub(" ",".",rownames(ciber))
identical(colnames(ciber),colnames(exp))
exp_ciber <- rbind(exp,ciber)
exp_ciber <- exp_ciber %>% t() %>% as.data.frame()
#install.packages("ggstatsplot")
library(ggstatsplot)
library(tidyverse)
ggscatterstats(data = exp_ciber, #要分析的数据
               y = BTK, #设置Y轴
               x = B.cells.naive,#设置X轴
               type = "nonparametric", 
               margins = "both",#是否显示 边缘，默认为true                                      
               xfill = "#01468b", #x轴边缘图形的颜色
               yfill = "#ee0000", #y轴边缘图形的颜色
               marginal.type = "densigram")#在图片坐标轴边缘添加图形类型


