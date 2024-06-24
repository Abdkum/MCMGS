rm(list = ls())
options(stringsAsFactors = F)

setwd("C:\\Users\\Administrator\\Desktop\\20230225.LUAD.ICD\\05.Drug sensitive")

#这里有两种方法，第一种：用pRRophetic；第二种：oncoPredict供选择

#方法一：pRRophetic
#05文件中有一个压缩包，用pRRophetic.gz进行本地安装，较大，建议本地安装，压缩包放在你的指定目录即可，安装时在压缩包名称前加上路径信息，如下：
#install.packages("C:\\Users\\Administrator\\Desktop\\20230225.LUAD.ICD\\04.Drug sensitive\\pRRophetic.gz", dependencies = TRUE,repos = NULL,type="source")
#药物敏感性预测
library(ggpubr)
library(pRRophetic)
library(ggplot2)
allDrugs=c("A.443654", "A.770041", "ABT.263", "ABT.888", "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "Z.LLNle.CHO", "ZM.447439")

#读取表达矩阵
data <- read.delim('C:\\Users\\Administrator\\Desktop\\20230225.LUAD.ICD\\03.Data.pre\\TCGA\\results\\tcga_dat_T.txt',
                      stringsAsFactors = F, 
                      check.names = F, 
                      row.names = 1)

data=as.matrix(data)
drug_ic50.dat<-data.frame()
for (drug in allDrugs) {
  print(drug)
  set.seed(12345)
  #预测药物敏感性
  senstivity=pRRopheticPredict(data, drug, selection=1)
  senstivity=senstivity[senstivity!="NaN"]
  senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
  tmpic50=as.data.frame(senstivity)
  colnames(tmpic50) <- drug
  if(nrow(drug_ic50.dat)==0){
    drug_ic50.dat=tmpic50
  }
  drug_ic50.dat <- cbind.data.frame(drug_ic50.dat, tmpic50)
}
write.table(drug_ic50.dat,file = 'results/drug_ic50.dat.txt',quote = F,sep='\t',row.names = T)
#drug_ic50.dat=read.delim("results/drug_ic50.dat.txt",stringsAsFactors = F, check.names = F, row.names = 1)
# 读取分组信息
dat_group <- read.csv('C:\\Users\\Administrator\\Desktop\\20230225.LUAD.ICD\\04.ConsensusClusterPlus\\results\\cluster.csv',
                      stringsAsFactors = F, 
                      check.names = F, 
                      row.names = 1)

#合并数据
drug_ic50.dat.res<-cbind.data.frame(dat_group,drug_ic50.dat)
drug_ic50.dat.res$Cluster=factor(drug_ic50.dat.res$Cluster,levels = c('clust1','clust2'))




drug.res<-data.frame()
for (drug in allDrugs){
  gr=levels(drug_ic50.dat.res$Cluster)
  dat1=as.numeric(na.omit(drug_ic50.dat.res[drug_ic50.dat.res$Cluster==gr[1],drug]))
  dat2=as.numeric(na.omit(drug_ic50.dat.res[drug_ic50.dat.res$Cluster==gr[2],drug]))
  fc=median(dat1)/median(dat2)-1
  #fc=mean(dat1)/mean(dat2)-1
  p=wilcox.test(x=dat1,y = dat2)$p.value
  if(nrow(drug.res)==0){
    drug.res=data.frame(x=drug,y=fc,p=p)
  }else{
    drug.res=rbind.data.frame(drug.res,data.frame(x=drug,y=fc,p=p))
  }
}
drug.res$FDR<-p.adjust(drug.res$p,method = 'BH')
head(drug.res)
drug.res$group<-ifelse(drug.res$p<0.05 & drug.res$FDR<0.05,'P < 0.05 and FDR < 0.05',ifelse(drug.res$p<0.05 & drug.res$FDR>0.05 ,'P < 0.05 and FDR > 0.05','P > 0.05 and FDR > 0.05'))
drug.res=drug.res[order(drug.res$y,decreasing = T),]
drug.res$x=factor(drug.res$x,levels = drug.res$x)
library(ggplot2)
#可以去除p>0.05,FDR>0.05的数据
del <- which(drug.res$group=="P > 0.05 and FDR > 0.05")
drug.res <- drug.res[-del,]
range(drug.res$y)
p1=ggplot(drug.res, aes(x, y, fill = group)) + 
  geom_bar(stat = 'identity') +
  theme_classic()+ylim(-2,2)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1),
        legend.position = c(0.8,0.8))+
  xlab('')+ylab('Median Estimated IC50 (C1) / Median Estimated IC50 (C2)')
p1
ggsave(p1,filename = "./results/drug.pdf",he=6,wi=20)
#选择p值和FDR有意义的药物绘制箱线图
table(drug.res$group)
drug.res.fit=drug.res[which(drug.res$group=='P < 0.05 and FDR < 0.05'),]

sig_point<-function(dat,leg,ylab,palette=ggsci::pal_lancet()(10)[3:4]){
  library(ggpubr)
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggplot(dat,aes(x=group,y = gene))+
    #geom_jitter(aes(colour=group),size=1,position = position_jitter(seed=1))+
    geom_boxplot(aes(fill=group),colour='black',alpha=1)+
    theme_classic()+
    xlab('')+ylab(ylab)+labs(fill=leg)+#scale_color_manual(values = palette)+
    scale_fill_manual(values = palette)+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")
  return(pp)
}
drug_ic50.dat.res[1:4,1:4]
dir.create('results/boxplot')
for (i in unique(as.character(drug.res.fit$x))){
  dd=sig_point(dat = drug_ic50.dat.res[,c("Cluster",i)],
               leg = 'Cluster',ylab = paste0('Drugs(IC50)'),
               palette = ggsci::pal_aaas()(10))
  ggsave(paste0('results/boxplot/',i,'_boxplot.pdf'),dd,height = 5,width = 5)
}



#方法二：oncoPredict
#此方法速度更快，绘制箱线图，自己选择药物，并修改boxplot药物名称
rm(list = ls())
options(stringsAsFactors = F)
# install.packages("oncoPredict") 
#报错缺哪个包，逐步安装就好
setwd("C:\\Users\\Administrator\\Desktop\\20230225.LUAD.ICD\\05.Drug sensitive")
library(oncoPredict)

# 直接运行函数，不用管
get_oncoPredict_res <- function(data = data,
                                traData = c('GDSC2', 'GDSC1', 'CTRP2')[1],
                                minNumSamples = 10) {
  library(oncoPredict)
  
  data <- as.matrix(data)
  if (traData == 'GDSC2') {
    traDataExp <- readRDS("DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
    traDataRes <- readRDS(file = "DataFiles/Training Data/GDSC2_Res.rds")
  } else if (traData == 'GDSC1') {
    traDataExp <- readRDS("DataFiles/Training Data/GDSC1_Expr (RMA Normalized and Log Transformed).rds")
    traDataRes <- readRDS(file = "DataFiles/Training Data/GDSC1_Res.rds")
  } else if (traData == 'CTRP2') {
    traDataExp <- readRDS("DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
    traDataExp <- log2(traDataExp + 1)
    traDataRes <- readRDS(file = "DataFiles/Training Data/CTRP2_Res.rds")
  } else {
    stop('请检查 traData 是否正确, 仅支持GDSC2, GDSC1, CTRP2 三种药物信息库')
  }
  
  calcPhenotype(trainingExprData = traDataExp,
                trainingPtype = traDataRes,
                testExprData = data,
                batchCorrect = 'eb',  #   "eb" for ComBat  
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = minNumSamples, 
                printOutput = TRUE, 
                removeLowVaringGenesFrom = 'rawData')
}


# 读取Data.pre准备好的肿瘤数据和分组
dat_exp <- read.delim('C:\\Users\\Administrator\\Desktop\\20230225.LUAD.ICD\\03.Data.pre\\TCGA\\results\\tcga_dat_T.txt',
                      stringsAsFactors = F, 
                      check.names = F, 
                      row.names = 1)
boxplot(dat_exp[, 1:5])
dat_exp[1:5, 1:5]


# 读取分组信息
dat_group <- read.csv('C:\\Users\\Administrator\\Desktop\\20230225.LUAD.ICD\\04.ConsensusClusterPlus\\results\\cluster.csv',
                      stringsAsFactors = F, 
                      check.names = F, 
                      row.names = 1)
head(dat_group)
table(dat_group$Cluster)

# 药物敏感性预测 ###############
# 默认GDSC2
get_oncoPredict_res(data = dat_exp,
                    traData = 'GDSC2')

# # GDSC1
# get_oncoPredict_res(data = dat_exp,
#                     traData = 'GDSC1')
# 
# # CTRP2
# get_oncoPredict_res(data = dat_exp,
#                     traData = 'CTRP2')


library(stringr)
GDSC2_res <- read.csv('calcPhenotype_Output/DrugPredictions.csv',
                      row.names = 1, check.names = F, stringsAsFactors = F)
GDSC2_res[1:5, 1:5]



# 展示画图 ##################
dat_GDSC2_res <- cbind(dat_group,
                       GDSC2_res)
dat_GDSC2_res[1:5, 1:6]
class(dat_GDSC2_res)
str(dat_GDSC2_res)

# unique(colnames(dat_GDSC2_res))

# 绘图 
library(ggsci)
library(ggplot2)
library(ggpubr)

#查看allDrugs中药物，选择感兴趣药物，如Cisplatin_1005
Cisplatin_1005 <- ggplot(dat_GDSC2_res, aes(x=Cluster, y=Cisplatin_1005)) +
  geom_boxplot(aes(fill = Cluster), position=position_dodge(0.9), outlier.colour = NA, notch = T) +
  stat_compare_means(aes(group=Cluster)) +
  scale_fill_d3()+
  theme_bw()+ 
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5,  hjust = 0.5,colour="black"),
        legend.position = 'top') + xlab('Cluster') + 
  ylab('Estimated IC50 (Cisplatin_1005)')#这里也要改
Cisplatin_1005
ggsave(p1,filename = "./results/Cisplatin_1005.pdf",he=6,wi=6)





       