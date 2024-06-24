
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("vioplot")


#引用包
library(limma)
library(vioplot)
immuneFile="CIBERSORT-Results.txt"      #免疫细胞浸润文件
cluFile="risk.all.txt"                   #分型结果文件
pFilter=0.05          #免疫细胞浸润结果的过滤条件
setwd("E:\\SCI\\NKcell\\NKCMGS\\12.FigureYa230immunelandscape\\imm\\21.clusterImmVioplot")     #设置工作目录

#读取免疫细胞浸润结果文件，并对数据进行整理
immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
immune=as.matrix(immune[,1:(ncol(immune)-3)])

#删掉正常样品
group=sapply(strsplit(row.names(immune),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
immune=immune[group==0,]
row.names(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", row.names(immune))
immune=avereps(immune)

#读取分型的文件
cluster=read.table(cluFile, header=F, sep="\t", row.names=1, check.names=F)

#按照分型对样品分组
lowName=row.names(cluster)[cluster[,1]==1]
highName=row.names(cluster)[cluster[,1]==2]

#提取不同分型的免疫细胞含量
lowImm=intersect(row.names(immune), lowName)
highImm=intersect(row.names(immune), highName)
rt=rbind(immune[lowImm,], immune[highImm,])
lowNum=length(lowImm)
highNum=length(highImm)

#绘制小提琴图
outTab=data.frame()
pdf("vioplot.pdf", width=13, height=9)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63), ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，低表达用蓝色表示，高表达用红色表示
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
for(i in 1:ncol(rt)){
	  if(sd(rt[1:lowNum,i])==0){
	    rt[1,i]=0.00001
	  }
	  if(sd(rt[(lowNum+1):(lowNum+highNum),i])==0){
	    rt[(lowNum+1),i]=0.00001
	  }
	  lowData=rt[1:lowNum,i]
	  highData=rt[(lowNum+1):(lowNum+highNum),i]
	  vioplot(lowData,at=3*(i-1),lty=1,add = T,col=bioCol[1])
	  vioplot(highData,at=3*(i-1)+1,lty=1,add = T,col=bioCol[2])
	  wilcoxTest=wilcox.test(lowData,highData)
	  p=wilcoxTest$p.value
	  if(p<pFilter){
	      cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
		  outTab=rbind(outTab,cellPvalue)
	  }
	  mx=max(c(lowData,highData))
	  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
	  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Low", "High"),
       lwd=4.5,bty="n",cex=1.5,
       col=bioCol[1:2])
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 0.9,srt = 45,pos=2)
dev.off()

#输出免疫细胞和p值表格文件
write.table(outTab,file="diff.result.txt",sep="\t",row.names=F,quote=F)

