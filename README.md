# 下载所需要的R包 ----------------------------------------------------------------
install.packages("VIM")
install.packages("rms")
install.packages("nomogramFormula")
install.packages("pROC")
install.packages("rmda")

# 加载需要的R包 -----------------------------------------------------------------
library(VIM)#缺失值可视化
library(rms)#拟合模型
library(nomogramFormula)#计算列线图得分
library(pROC)#绘制ROC曲线，计算AUC和95%置信区间
library(rmda)#临床决策曲线和临床影响曲线

# 载入示例数据 ------------------------------------------------------------------
mydata <- read.csv("训练集-csv.csv",row.names = 1,sep = ",",header = T)

mydata_copy <- mydata#数据备份


#生存变量的处理
table(mydata$Outcome)#检查生存状态
mydata$Outcome <- factor(mydata$Outcome,level=c(0,1),labels=c("DB without AF","DB with AF"))

#二分类变量的转换（计算评分时不要运行）
mydata$CRRT <- factor(mydata$CRRT,level=c(0,1),labels=c("No","Yes"))
mydata$CVC <- factor(mydata$CVC,level=c(0,1),labels=c("No","Yes"))
mydata$Broad.spectrum.antibiotics <- factor(mydata$Broad.spectrum.antibiotics,level=c(0,1),labels=c("No","Yes"))
mydata$BDG <- factor(mydata$BDG,level=c(0,1),labels=c("No","Yes"))


# 开始构建logistic回归模型 ------------------------------------------------------------
#数据打包 
dd = datadist(mydata)
option <- options(datadist = "dd")

#拟合模型
colnames(mydata)#查看列名，选择你要构建模型的变量
formula <- as.formula(Outcome ~ Hypertension + Non.Hispanic.White  + Male + Age+ Coronary.Heart.Disease+Chronic.Bronchitis+NYHA.III.IV+Hyperthyroidism+Infection+Renal.Insufficiency )
model <- lrm(formula, # 回归模型的公式,指定自变量和因变量
           data = mydata, # 包含所有变量的数据框
           x=TRUE, # logistic回归也称为"广义线性模型",这个参数指定响应变量的二分类
           y=TRUE) # 参数y也是指定因变量是二分类的
model#查看模型具体情况

OR <- exp(model$coefficients)#计算Logistic回归模型中每个自变量的比值比(Odds Ratio, OR)
OR
# 计算每个变量的OR值
OR <- exp(model$coefficients)

# 计算每个变量的标准误差
SE <- sqrt(diag(model$var))  # diag提取模型的方差协方差矩阵

# 计算95%置信区间 (Wald 置信区间)
lower_bound <- exp(model$coefficients - 1.96 * SE)  # 下限
upper_bound <- exp(model$coefficients + 1.96 * SE)  # 上限

# 结果输出
OR_table <- data.frame(OR = OR, Lower_CI = lower_bound, Upper_CI = upper_bound)
print(OR_table)
# 结果可视化 -------------------------------------------------------------------
Nomogram_1 <- nomogram(model,
                       fun = function(x)1/(1+exp(-x)),
                       lp=F,
                       fun.at = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
                       funlabel = "Risk")
plot(Nomogram_1)
# 答疑：
# 在Nomogram图的绘制中,线段长度同时取决于两个因素:
# 1) 变量的回归系数:系数绝对值越大的变量,其线段越长,其他条件相同。
# 2) 变量的取值范围:范围越大的变量,其线段也越长,其他条件相同。

#调整一些参数
plot(Nomogram_1,
     #模型的名称
     xfrac = .35,
     #变量与图形的占比（调整变量与坐标抽距离）
     cex.var = 1.6,
     #变量字体加粗
     cex.axis = 1.4,
     #数轴：字体的大小
     tcl = -0.5,
     #数轴：刻度的长度
     lmgp = 0.3,
     #数轴：文字与刻度的距离
     label.every = 1,
     #数轴：刻度下的文字，1=连续显示，2=隔一个显示一个
     col.grid = gray(c(0.8,0.95)))

#计算列线图评分(要计算评分，数据中不能含有因子变量)
options(option) 
#设置options()函数相关参数
results <- formula_rd(nomogram = Nomogram_1)
#使用formula_rd()函数基于nomogram对象Nomogram_1生成公式formula和相应的数据框rd
mydata$points<-points_cal(formula = results$formula,rd=mydata)
#使用points_cal()函数基于formula公式和原始数据框mydata计算points列,并添加到mydata
head(mydata) 
#查看mydata的数据框的前6行

#可以根据列线图评分（中位数、均值、四分位数）对所有样本进行分组

#可以根据列线图评分（中位数、均值、四分位数）对所有样本进行分组

# 模型评估 -----------------------------------------------------------------
# ROC曲线 -------------------------------------------------------------------
#AUC计算
model_ROC <- glm(formula,data = mydata,family = binomial())
#lrm()函数属于Design包,专门用于拟合Logistic回归模型。
#glm()函数属于stats包,是用于拟合广义线性模型(Generalized Linear Models)的泛用函数,可以拟合Logistic回归、Poisson回归等多种模型。

#type=“response”给出具体的预测概率，而type=“class”按规定的阈值给出分类
mydata$predvalue <- predict(model_ROC,type="response")
# 计算AUC和95%CI
ROC <- roc(mydata$Outcome,mydata$predvalue)
auc(ROC)
ci(auc(ROC))
# 确保ROC曲线数据按照specificity升序排列
ROC_data <- data.frame(
  specificity = ROC$specificities, 
  sensitivity = ROC$sensitivities
)

# 按照specificity升序排列
ROC_data <- ROC_data[order(ROC_data$specificity), ]

# 重新绘制ROC曲线，确保X轴和Y轴的标签、坐标轴范围正确
plot(1 - ROC_data$specificity,  
     ROC_data$sensitivity, 
     type = "l", 
     col = "#00008B",  # 深蓝色, 
     lty = 1,
     xlab = "1 - Specificity",
     ylab = "Sensitivity",
     lwd = 2,
     xlim = c(0, 1),  # X轴范围从0到1
     ylim = c(0, 1))  # Y轴范围从0到1
# 添加对角线（随机分类器的线）
abline(0, 1, col = "black", lty = 2)

# 修改图例部分
legend("bottomright", 
       legend = c("AUC: 0.83 (95% CI 0.85 - 0.88)"),  # 图例内容
       lty = 1,  # 曲线类型：实线
       lwd = 2,  # 曲线宽度
       col = "#00008B",  # 曲线颜色：深蓝色
       bty = "n",  # 无边框
       cex = 0.8)  # 字体大小
# 计算AUC和最大约登指数
ROC <- roc(mydata$Outcome, mydata$predvalue)

# 计算AUC
AUC_value <- auc(ROC)
AUC_ci <- ci(auc(ROC))  # 计算95%置信区间

# 输出AUC和95%置信区间
cat("AUC of Nomogram:", AUC_value, "\n")
cat("95% CI for AUC:", AUC_ci, "\n")

# 计算最大约登指数 (Youden's Index)
coords_result <- coords(ROC, "best", ret = c("threshold", "sensitivity", "specificity", "youden"))

# 输出最大约登指数的阈值、敏感性、特异性以及约登指数
cat("Best Threshold:", coords_result$threshold, "\n")
cat("Sensitivity at Best Threshold:", coords_result$sensitivity, "\n")
cat("Specificity at Best Threshold:", coords_result$specificity, "\n")
cat("Youden's Index at Best Threshold:", coords_result$youden, "\n")
# 载入测试集数据
mydata2 <- read.csv("同济医院-测试集-1226-2257.csv", row.names = 1, sep = ",", header = T)




# 计算测试集数据的预测值
# 使用训练集模型（model）对测试集进行预测
# 使用训练集模型（model）对测试集进行预测
mydata2$predvalue <- predict(model, newdata = mydata2, type = "fitted")


# 计算新的ROC曲线
ROC_test <- roc(mydata2$Outcome, mydata2$predvalue)

# 输出AUC和95%置信区间
AUC_test <- auc(ROC_test)
AUC_test_ci <- ci(auc(ROC_test))  # 计算95%置信区间
cat("AUC of Test Set:", AUC_test, "\n")
cat("95% CI for AUC of Test Set:", AUC_test_ci, "\n")

# 绘制测试集的ROC曲线
plot(1 - ROC_test$specificities,  
     ROC_test$sensitivities, type = "l", 
     col = "blue", 
     lty = 1,
     xlab = "1-Specificity",
     ylab = "Sensitivity",
     lwd = 2)

# 添加对角线
abline(0, 1)

# 添加图例
legend(0.45, 0.05, 
       c(paste("AUC of Test Set:", round(AUC_test, 3), "(95% CI", round(AUC_test_ci[1], 3), "-", round(AUC_test_ci[2], 3), ")")),
       lty = c(1),
       lwd = c(2),
       col = c("blue"),
       bty = "0")
# 计算灵敏度、特异度和约登指数
coords_result <- coords(ROC_test, "best", ret = c("threshold", "sensitivity", "specificity", "youden"))

# 输出灵敏度、特异度和约登指数
cat("Best Threshold:", coords_result$threshold, "\n")
cat("Sensitivity at Best Threshold:", coords_result$sensitivity, "\n")
cat("Specificity at Best Threshold:", coords_result$specificity, "\n")
cat("Youden's Index at Best Threshold:", coords_result$youden, "\n")
# 输出训练集数据（mydata）为CSV文件
write.csv(mydata, file = "训练集-csv-mydata.csv", row.names = FALSE)

# 输出测试集数据（mydata2）为CSV文件
write.csv(mydata2, file = "验证集-csv-mydata2.csv", row.names = FALSE)

