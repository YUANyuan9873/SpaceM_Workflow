## 设定工作目录
setwd("./")
## 数据路径
path <- "./Data/"
## R包加载
## library(Cardinal)

## neg 数据
file_neg <- list.files(path = paste0(path,"1 imzml/"), pattern = "-neg.imzML")
## pos 数据
file_pos <- list.files(path = paste0(path,"1 imzml/"), pattern = "-pos.imzML")
## sample 名字
sample_name <- c("liver-1","liver-2","liver-3","liver-4")

## 读入imzML数据
neg_1 <- readMSIData(file = paste0(path, "1 imzml/",file_neg[1]),
                     resolution = 5, 
                     units = "ppm",
                     mass.range = c(70,1000))



## SSCC聚类：空间质心收缩聚类分析
neg_1_SSCC <- spatialShrunkenCentroids(neg_1, 
                                       r = 1, # 邻域平滑半径`r`：应根据你的数据集中空间区域的大小和颗粒度来选择
                                       k = 8, # 最大聚类数
                                       s = 3, # 缩减或稀疏参数`s`：这个数字越大，用于确定最终分割的峰值就越少
                                       iter.max = 30)
## SSCC结果存储为rds文件
## 结果输出路径设置
out_path <- paste0(path,"2 sscc-rds/sscc/")
if (dir.exists(out_path)==FALSE) {
  dir.create(out_path, recursive = TRUE)
} else {
  print("it has been already exited")}
saveRDS(neg_1_SSCC, file = paste0(path,"2 sscc-rds/sscc/",sample_name[1],"-neg-sscc.rds"))
neg_1_SSCC <- readRDS(file = paste0(path,"2 sscc-rds/sscc/",sample_name[1],"-neg-sscc.rds"))
## SSCC可视化结果展示
image(neg_1_SSCC)
## SSCC聚类T统计量棒状图
plot(neg_1_SSCC, values="statistic", lwd=2)

## 质谱图
plot(neg_1)


## 平均质谱峰
mean_feature <- summarizeFeatures(neg_1, "mean")
plot(mean_feature)

## 多离子成像图
image(neg_1, mz = mz(neg_1)[1:5])
## 多颜色通道成像图
image(neg_1, mz = mz(neg_1)[1:3],
      col = c("#d32f2f", "#303f9f","#fbc02d"),
      superpose = TRUE)

## 离子共定位分析
coloc <- colocalized(neg_1, mz = mz(neg_1)[1])
coloc
image(neg_1, mz = coloc$mz[1:3], layout = c(1,3))

