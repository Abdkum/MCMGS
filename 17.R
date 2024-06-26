
#install.packages("reshape2")
#install.packages("ggpubr")


#引用包
library(reshape2)
library(ggpubr)
riskFile="tcgarisk.txt"           #风险文件
estimateFile="TMEscores.txt"      #肿瘤微环境打分文件
setwd("E:\\SCI\\NKcell\\NKCMGS\\新建文件夹\\tme47.estimateVioplot")      #设置工作目录

#读取风险文件
Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
Risk$risk=factor(Risk$risk, levels=c("low","high"))

#读取肿瘤微环境打分文件
score=read.table(estimateFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]
rownames(score)=gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
score=score[row.names(Risk),,drop=F]

#数据合并
rt=cbind(Risk[,"risk",drop=F], score)

#将合并后的数据转换为ggplot2的输入文件
data=melt(rt, id.vars=c("risk"))
colnames(data)=c("Risk", "scoreType", "Score")

#绘制小提琴图
p=ggviolin(data, x="scoreType", y="Score", fill = "Risk",
	     xlab="",
	     ylab="TME score",
	     legend.title="Risk",
	     add = "boxplot", add.params = list(color="white"),
	     palette = c("blue","red"), width=1)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Risk),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出图形
pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()

