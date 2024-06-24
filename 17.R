
#install.packages("reshape2")
#install.packages("ggpubr")


#���ð�
library(reshape2)
library(ggpubr)
riskFile="tcgarisk.txt"           #�����ļ�
estimateFile="TMEscores.txt"      #����΢��������ļ�
setwd("E:\\SCI\\NKcell\\NKCMGS\\�½��ļ���\\tme47.estimateVioplot")      #���ù���Ŀ¼

#��ȡ�����ļ�
Risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
Risk$risk=factor(Risk$risk, levels=c("low","high"))

#��ȡ����΢��������ļ�
score=read.table(estimateFile, header=T, sep="\t", check.names=F, row.names=1)
score=score[,1:3]
rownames(score)=gsub("(.*?)\\_(.*?)", "\\2", rownames(score))
score=score[row.names(Risk),,drop=F]

#���ݺϲ�
rt=cbind(Risk[,"risk",drop=F], score)

#���ϲ��������ת��Ϊggplot2�������ļ�
data=melt(rt, id.vars=c("risk"))
colnames(data)=c("Risk", "scoreType", "Score")

#����С����ͼ
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

#���ͼ��
pdf(file="vioplot.pdf", width=6, height=5)
print(p1)
dev.off()
