
#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      #p值过滤条件
qvalueFilter=0.05      #矫正后的p值过滤条件

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd("E:\\SCI\\NKcell\\NKCMGS\\10.KEGG富集结果的圆圈图")            #设置工作目录
rt=read.table("diffgene.txt", header=F, sep="\t", check.names=F)     #读取输入文件

#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存富集结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=30, color=colorSel)
dev.off()


