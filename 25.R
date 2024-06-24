
#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#���ð�
library(limma)
library(survival)
library(survminer)
library(timeROC)
riskFile="allRisk.txt"     #�����ļ�
expFile="symbol.txt"       #���������ļ�
cliFile="time.txt"         #�ٴ������ļ�
geneFiles=c("Ran ZhangSig.txt", "Guoda SongSig.txt","Zhuofan MouSig.txt", "Chenghao WenSig.txt","Ding HuSig.txt", "Zhengtong LvSig.txt","NingShaoSig.txt", "WangliMeiSig.txt","DaixingHuSig.txt", "HuanLiuSig.txt")       #ģ�ͻ����ļ�
setwd("E:\\SCI\\NKcell\\NKCMGS\\5.1modelcompare\\33.modelCompare")     #���ù���Ŀ¼

#��ȡ�����ļ������������ļ�����
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
data=data[,group==0]

#��ȡ�����ļ�
riskRT=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
riskRT=riskRT[,c("futime","fustat","riskScore")]
colnames(riskRT)=c("futime","fustat","GILncSig")

for(i in geneFiles){
  #��ȡ�����б�
  header=unlist(strsplit(i, "\\."))
  gene=read.table(i, header=F, sep="\t", check.names=F)
  sameGene=intersect(as.vector(gene[,1]), row.names(data))
  data1=data[sameGene,]
  colnames(data1)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data1))
  data1=t(data1)
  data1=avereps(data1)
  
  #��ȡ��������
  cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)
  cli$futime=cli$futime/1
  
  #���ݺϲ���������
  sameSample=intersect(row.names(data1), row.names(cli))
  data1=data1[sameSample,]
  cli=cli[sameSample,]
  data1=cbind(cli,data1)
  
  #����÷�
  multiCox=coxph(Surv(futime, fustat) ~ ., data = data1)
  riskScore=predict(multiCox,type="risk", newdata=data1)
  data1=cbind(data1, riskScore)
  data1=data1[row.names(riskRT),]
  riskRT=cbind(riskRT, data1[,"riskScore"])
  colnames(riskRT)[ncol(riskRT)]=header[[1]]
}

#����ROC����
predictTime=1
bioCol=rainbow(ncol(riskRT)-2)
aucText=c()
pdf(file="ROC.pdf", width=6, height=6)
i=3
ROC_rt=timeROC(T=riskRT$futime,delta=riskRT$fustat,marker=riskRT[,i],cause=1,weighting='aalen',times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0(colnames(riskRT)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
for(i in 4:ncol(riskRT)){
  ROC_rt=timeROC(T=riskRT$futime,delta=riskRT$fustat,marker=riskRT[,i],cause=1,weighting='aalen',times=c(predictTime),ROC=TRUE)
  plot(ROC_rt,time=predictTime,col=bioCol[i-2],title=FALSE,lwd=2,add=TRUE)
  aucText=c(aucText,paste0(colnames(riskRT)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-2)])
dev.off()