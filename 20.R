
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSEABase")

#install.packages("pheatmap")
#install.packages("reshape2")


#引用包
library(limma)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(reshape2)

expFile="symbol.txt"         #表达数据文件
gmtFile="immune.gmt"         #免疫功能的基因集文件
riskFile="tcgarisk.txt"      #风险文件
setwd("E:\\SCI\\NKcell\\NK\\29.immFunction")       #设置工作目录

#读取表达输入文件，并对输入文件处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=mat[rowMeans(mat)>0,]
	
#读取数据集文件
geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssgsea分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
	return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
data=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(data), data)
write.table(ssgseaOut, file="immFunScore.txt", sep="\t", quote=F, col.names=F)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

#读取风险文件,获取高低风险
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
lowSample=row.names(risk[risk$risk=="low",])
highSample=row.names(risk[risk$risk=="high",])
lowData=data[,lowSample]
highData=data[,highSample]
data=cbind(lowData, highData)
conNum=ncol(lowData)        #低风险组样品数目
treatNum=ncol(highData)     #高风险组样品数目
sampleType=c(rep(1,conNum), rep(2,treatNum))

#差异分析
sigVec=c()
for(i in row.names(data)){
	test=wilcox.test(data[i,] ~ sampleType)
	pvalue=test$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	sigVec=c(sigVec, paste0(i, Sig))
}
row.names(data)=sigVec

#热图可视化
Type=c(rep("Low risk",conNum), rep("High risk",treatNum))
Type=factor(Type, levels=c("Low risk", "High risk"))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf("heatmap.pdf", width=8, height=4.6)
pheatmap(data,
         annotation=Type,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=T,
         fontsize=7,
         fontsize_row=7,
         fontsize_col=7)
dev.off()


