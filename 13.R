## 输入数据
第一列是生存时间、第二列是追踪情况、第三列以后是基因表达水平

```{r,message=FALSE,warning=FALSE}
svdata <- read.csv("easy_input.csv",header = T,row.names = 1)
dim(svdata) #一共18个基因
#查看前3个sample的前4个基因的表达量
head(svdata[1:3,1:6])
```

## 开始画图

找best separation用的是survminer的函数，结果与FigureYa4bestSeparation一致

```{r}
##加载必需的包
library(survival)
library(survminer)

##对数据集的基因进行bestSeparation统计
res.cut <- surv_cutpoint(svdata, time = "futime", 
                         event = "fustat", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.25) #默认组内sample不能低于30%

##按照bestSeparation分高低表达
res.cat <- surv_categorize(res.cut)

##统计作图
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  
  ##计算HR以及95%CI
  ##修改分组参照
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  
  #只画出p value<=0.05的基因，如果不想筛选，就删掉下面这行
  #if (p.val>0.05) next
  
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  
  #按照基因表达量从低到高排序，便于取出分界表达量
  svsort <- svdata[order(svdata[,i]),]
  
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
           #ggtheme = theme_bw(), #想要网格就运行这行
           conf.int = F, #不画置信区间，想画置信区间就把F改成T
           #conf.int.style = "step",#置信区间的类型，还可改为ribbon
           censor = F, #不显示观察值所在的位置
           palette = c("#D95F02","#1B9E77"), #线的颜色对应高、低
           
           legend.title = i,#基因名写在图例题目的位置
           font.legend = 12,#图例的字体大小
           font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小

           #在图例上标出高低分界点的表达量，和组内sample数量
           legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                           paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
           
           #在左下角标出pvalue、HR、95% CI
           #太小的p value标为p < 0.001
           pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                             paste("p = ",round(p.val,3), sep = "")),
                        HR, CI, sep = "\n"))
    
  #如果想要一个图保存为一个pdf文件，就把下面这行前面的“#”删掉
  ggsave(paste0(i,".pdf"),width = 4,height = 4)
}

length(pl)
```

只保留p<=0.05的10个基因，筛掉了p>0.05的8个基因

### 批量出图

用survminer包自带的函数组图

```{r,fig.width=16,fig.height=20}
res <- arrange_ggsurvplots(pl, 
                           print = T,
                           ncol = 4, nrow = 5)#每页纸画几列几行

#保存到pdf文件
ggsave("bestSurvPlot.pdf",res,width = 16,height = 20)
```

```{r}
#source("rcode.R")
sessionInfo()