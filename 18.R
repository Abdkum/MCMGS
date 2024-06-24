
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)
riskFile="risk.txt"      #风险输入文件
immFile="infiltration_estimation_for_tcga.csv"     #免疫细胞浸润文件
setwd("D:\\biowolf\\FerrLnc\\27.immHeatmap")     #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)

#读取免疫细胞浸润文件
immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)

#病人风险值和免疫细胞合并
sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, c("risk", "riskScore")]
immune=immune[sameSample,]
data=cbind(risk, immune)

#高低风险组免疫差异分析
outTab=data.frame()
sigCell=c("risk","riskScore")
for(i in colnames(data)[3:ncol(data)]){
	if(sd(data[,i])<0.001){next}
	wilcoxTest=wilcox.test(data[,i] ~ data[,"risk"])
	pvalue=wilcoxTest$p.value
	if(wilcoxTest$p.value<0.05){
		outTab=rbind(outTab,cbind(immune=i, pvalue))
		sigCell=c(sigCell, i)
	}
}
write.table(file="immuneCor.txt", outTab, sep="\t", quote=F, row.names=F)

#热图数据
data=data[,sigCell]
data=data[order(data[,"riskScore"]),]
annCol=data[,1:2]
annCol[,"risk"]=factor(annCol[,"risk"], unique(annCol[,"risk"]))
data=t(data[,(3:ncol(data))])
annRow=sapply(strsplit(rownames(data),"_"), '[', 2)
annRow=as.data.frame(annRow)
row.names(annRow)=row.names(data)
colnames(annRow)=c("Methods")
annRow[,"Methods"]=factor(annRow[,"Methods"], unique(annRow[,"Methods"]))
gapCol=as.vector(cumsum(table(annCol[,"risk"])))
gapRow=as.vector(cumsum(table(annRow[,"Methods"])))

#定义颜色
risk=c("blue", "red")
names(risk)=c("low", "high")
ann_colors=list(risk=risk)

#热图可视化
pdf("immHeatmap.pdf", width=9, height=6)
pheatmap(data,
         annotation=annCol,
         annotation_row=annRow,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =F,
         gaps_row=gapRow,
         gaps_col=gapCol,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=6,
         fontsize_row=5,
         fontsize_col=6)
dev.off()

