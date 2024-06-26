
#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


#引用包
library(survival)
library(survminer)
library(timeROC)

riskFile="tcgarisk.txt"     #风险文件
cliFile="TCGACLINICAL.txt"      #临床数据文件
setwd("E:\\SCI\\NKcell\\NKCMGS\\7.3ROC")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制1 3 5年的ROC曲线######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
	           marker=risk$riskScore,cause=1,
	           weighting='aalen',
	           times=c(1,3,5),ROC=TRUE)
pdf(file="ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
	   c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	     paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	     paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	   col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()


######绘制临床的ROC曲线######
predictTime=1     #定义预测年限
aucText=c()
pdf(file="cliROC.pdf", width=5, height=5)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#对临床数据进行循环，绘制临床数据的ROC曲线
for(i in 4:ncol(rt)){
	ROC_rt=timeROC(T=rt$futime,
				   delta=rt$fustat,
				   marker=rt[,i], cause=1,
				   weighting='aalen',
				   times=c(predictTime),ROC=TRUE)
	plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
	aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()

