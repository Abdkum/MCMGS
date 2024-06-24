# 安装所需的R包
# install.packages("survival")
# install.packages("survminer")
# install.packages("timeROC")

# 加载R包
library(survival)
library(survminer)
library(timeROC)

# 设置工作目录
setwd("D:\\WeChat Files\\wxid_7n5sysorplir21\\FileStorage\\File\\2024-06")

# 定义绘制ROC曲线的函数并计算AUC的置信区间
bioROC <- function(inputFile = NULL, rocFile = NULL) {
  # 读取输入文件
  rt <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
  
  # 计算ROC曲线
  ROC_rt <- timeROC(T = rt$futime, delta = rt$fustat,
                    marker = rt$riskScore, cause = 1,
                    weighting = 'marginal',
                    times = c(1, 3, 5), ROC = TRUE, iid = TRUE)
  
  # 计算AUC的置信区间
  ci_1 <- confint(ROC_rt, times = 1)$CI_AUC
  ci_3 <- confint(ROC_rt, times = 3)$CI_AUC
  ci_5 <- confint(ROC_rt, times = 5)$CI_AUC
  
  # 提取置信区间的上下限
  ci_1_lower <- ci_1[1, 1] / 100
  ci_1_upper <- ci_1[1, 2] / 100
  ci_3_lower <- ci_3[2, 1] / 100
  ci_3_upper <- ci_3[2, 2] / 100
  ci_5_lower <- ci_5[3, 1] / 100
  ci_5_upper <- ci_5[3, 2] / 100
  
  # 打开PDF设备
  pdf(file = rocFile, width = 5, height = 5)
  
  # 绘制ROC曲线
  plot(ROC_rt, time = 1, col = 'green', title = FALSE, lwd = 2)
  plot(ROC_rt, time = 3, col = 'blue', add = TRUE, title = FALSE, lwd = 2)
  plot(ROC_rt, time = 5, col = 'red', add = TRUE, title = FALSE, lwd = 2)
  
  # 设置图例字体大小
  par(cex = 0.7)
  
  # 添加图例，显示AUC值和置信区间
  legend('bottomright',
         c(paste0('1 year AUC: ', sprintf("%.03f", ROC_rt$AUC[1]),
                  " (95% CI: ", sprintf("%.03f", ci_1_lower), "-", sprintf("%.03f", ci_1_upper), ")"),
           paste0('3 year AUC: ', sprintf("%.03f", ROC_rt$AUC[2]),
                  " (95% CI: ", sprintf("%.03f", ci_3_lower), "-", sprintf("%.03f", ci_3_upper), ")"),
           paste0('5 year AUC: ', sprintf("%.03f", ROC_rt$AUC[3]),
                  " (95% CI: ", sprintf("%.03f", ci_5_lower), "-", sprintf("%.03f", ci_5_upper), ")")),
         col = c("green", 'blue', 'red'), lwd = 2, bty = 'n')
  
  # 关闭PDF设备
  dev.off()
}

# 调用函数，绘制ROC曲线并计算AUC的置信区间
bioROC(inputFile = "tcga.txt", rocFile = "tcga.pdf")
