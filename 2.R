### 输入
包含基因表达水平、生存时间、追踪情况等三列的文件，测试用文件为dsg1.txt

```{r,message=FALSE,warning=FALSE}
#install.packages(c("survival","survminer","ggplot2"))
library(survival)
library("survminer")
library(ggplot2)

#读入测试文件dsg1.txt
svdata<-read.table("tcga.txt",header=T,as.is=T)
head(svdata)
#对行按照基因表达水平排序，默认从低到高
sortsv<-svdata[order(svdata$expression),]
```
### 输出
中位值分组的生存曲线、最佳分组生存曲线、遍历所有分组情况下的P值和Hazard Ratio的分布情况

### Median separation
```{r,message=FALSE,warning=FALSE}
#先根据表达水平的中位值分组，画生存曲线，保存
ssdf<-cbind(sortsv,data.frame(gp=ifelse(sortsv$expression>median(sortsv$expression),"high","low")))
fit<-survfit(Surv(futime,fustat)~gp,data=ssdf)
pdf(file="medianSeparation.pdf")
sc_median<-ggsurvplot(fit, linetype = "strata", conf.int = F, pval = TRUE,palette = c("#D95F02","#1B9E77"),legend.title="",legend=c(0.7,0.9),legend.labs=c("High-expression","low-expression"))
sc_median
dev.off()
sc_median
```

遍历所有分组情况，计算P值和Hazard Ratio，p值用于判断分组之间差异是否显著，而Hazard Ratio用于衡量分组之间的差异程度
```{r,message=FALSE,warning=FALSE}
pvals<-c()
hrs<-c()
for(n in 1:(nrow(sortsv)-1)){
ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(n,nrow(sortsv)-n))))
diff<-survdiff(Surv(futime,fustat)~gp,data=ssdf,rho = 0)
pv<-pchisq(diff$chisq,length(diff$n)-1,lower.tail=FALSE)
pvals<-c(pvals,pv)
hr<-diff$obs[1]*diff$exp[2]/(diff$obs[2]*diff$exp[1])
hrs<-c(hrs,hr)
}
```

展示所有分组情况下的P值和Hazard Ratio的分布情况，水平虚线标记位置的P值为0.05，两条竖直虚线标记的HR为0.5和2
```{r,message=FALSE,warning=FALSE}
fd<-data.frame(Tag=1:(nrow(sortsv)-1),HR=hrs,Pvalue=pvals)
head(fd)
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+geom_point(shape=21,aes(fill=Tag))+scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+geom_vline(xintercept=c(-1,1),linetype="dashed")+annotate("text",y=-log10(pvals[which.min(pvals)]),x=log2(hrs[which.min(pvals)]),label="min-Pvalue")
ggsave(file="Pvalue_Hazard-Ratio.pdf")
```

### Best separation
虚线左上角区域的点p值最小，是最佳的分组方式，分组情况如下
```{r,message=FALSE,warning=FALSE}
ggplot(fd,aes(x=log2(HR),y=-log10(Pvalue)))+geom_point(shape=21,aes(fill=Tag))+scale_fill_gradientn(colours=c("#4477AA","#77AADD","#FFFFFF","#EE9988","#BB4444"),guide="legend")+geom_hline(yintercept=(-log10(0.05)),linetype="dashed")+geom_vline(xintercept=c(-1,1),linetype="dashed")+geom_text(aes(label=paste(Tag,nrow(sortsv)-Tag,sep=":")),vjust=-1)
```

画出最佳分组的生存曲线
```{r,message=FALSE,warning=FALSE}
ssdf<-cbind(sortsv,data.frame(gp=rep(c("low","high"),c(which.min(pvals),nrow(sortsv)-which.min(pvals)))))
fit<-survfit(Surv(futime,fustat)~gp,data=ssdf)
pdf(file="bestSeparation.pdf")
sc_minp<-ggsurvplot(fit, linetype = "strata", conf.int = F, pval = TRUE,palette = c("#D95F02","#1B9E77"),legend.title="",legend=c(0.7,0.9),legend.labs=c("High-expression","low-expression"))
sc_minp
dev.off()
sc_minp