## 数据路径
path <- "./Data/"

## 结果输出路径设置
out_path <- paste0(path,"3 Diffresult/")
if (dir.exists(out_path)==FALSE) {
  dir.create(out_path, recursive = TRUE)
} else {
  print("it has been already exited")}


## neg 数据
file_neg <- list.files(path = paste0(path,"1 imzml/"), pattern = "-neg.imzML")

## sample 名字
sample_name <- unlist(strsplit(file_neg,"-neg.imzML", fixed = TRUE))

## 读取选区信息和分组信息
ROI_select <- read.xlsx(paste0(path,"Registration.xlsx"),sheet = "聚类选区")
names(ROI_select)[2:4] <- c("AreaName", "Mode", "Cluster")
Group_data <- read.xlsx(paste0(path,"Registration.xlsx"),sheet = "分组信息")
names(Group_data)[1:3] <- c("Group", "Sample", "AreaName")

## 把选区信息和分组信息整合在一起
select_data <- merge(ROI_select,
                     Group_data,
                     by = "AreaName")

## 读入SSCC结果
ssc_neg <- readRDS(file = paste0(path,"2 sscc-rds/sscc/",sample_name[1],"-neg-sscc.rds"))
ssccluster_neg <- data.frame(pixelData(ssc_neg)@listData,
                             class = ssc_neg$class, 
                             select = T)

# 选择要用到的cluster
ssc_neg_select <- subset(ssccluster_neg, class == ROI_select$Cluster[3] |
                           class == ROI_select$Cluster[4]) 
ssc_neg_select$Pixel <- paste0(sample_name[1], "-neg-",
                               ssc_neg_select$X3DPositionX, "-",
                               ssc_neg_select$X3DPositionY)


## 开始匹配数据
## 先读入imzML数据
neg_1 <- readMSIData(file = paste0(path, "1 imzml/",file_neg[1]),
                     resolution = 5, 
                     units = "ppm",
                     mass.range = c(70,1000))

## 提取丰度矩阵
spectradata <- as.matrix(spectra(neg_1, "intensity"))
spectradata <- as.data.frame(spectradata)

## 命名：样本名-mode-x坐标-y坐标
colnames(spectradata) <- paste(sample_name[1],"neg",
                               coord(neg_1)$x, coord(neg_1)$y,
                               sep = "-")

## mz保留5位小数
row.names(spectradata) <- format(mz(neg_1), nsmall = 5, trim = T)
spectradata_t <- as.data.frame(t(spectradata))
spectradata_t$Pixel <- rownames(spectradata_t)

## 将cluster和丰度矩阵整合在一起
neg_data <- merge(ssc_neg_select,
                  spectradata_t, 
                  by = "Pixel")
rownames(neg_data) <- neg_data$Pixel

## 选取第一列和第四列
meta_data_neg <- neg_data[,c(1,4)]
names(meta_data_neg)[2] <- "Cluster"
select_area_neg <- subset(select_data, Mode == "neg")
meta_data_neg <- merge(select_area_neg[,c(4,6)],
                         meta_data_neg,
                         by = "Cluster")
rownames(meta_data_neg) <- meta_data_neg$Pixel

## 形成用于分析的数据矩阵
data_for_ana_neg <- neg_data[,6:length(neg_data)]

## 开始分析
## PCA

# 要求行为样本，列为变量
pca_data <- prcomp(data_for_ana_neg, scale. = TRUE) 
# 提取score数据
data_pca_data <- as.data.frame(pca_data$x) 
data_pca_data$Pixel <- rownames(data_pca_data)
data_pca_data <- merge(meta_data_neg[,2:3],
                       data_pca_data,
                       by = "Pixel")
pca_score_data <- data_pca_data[,c(1:4)]
write.csv(pca_score_data, file = paste(out_path,'pca_score_data.csv',sep = "/"))

## 计算loading
pca_loading_data <- as.data.frame(summary(pca_data)[2])[,c(1:2)]
# 这里只选取了前两个PC的loading
names(pca_loading_data)[1:2] = c("PC1","PC2")
pca_loading_data$var <- rownames(pca_loading_data)
write.csv(pca_loading_data, file = paste(out_path,'pca_loading_data.csv',sep = "/"))

## 计算每个主成分对方差的解释度
pca_var_data <- pca_data$sdev^2 %>% as.data.frame()
# 计算各主成分所占百分比
pca_var_data$var <- round(pca_var_data$. / sum(pca_var_data) * 100, 2) 
pca_var_data$pc <- colnames(data_pca_data)[3:(ncol(data_pca_data))]

## 绘制碎石图看每个主成分的解释量
ggplot(pca_var_data[c(1:5),], aes(pc, var, fill = pc)) +
  geom_bar(stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw() +
  labs(x = '主成分',
       y = '主成分解释量（%）')

## PCA可视化
p_score_data <- ggplot(pca_score_data, 
                       aes(x = PC1,
                           y = PC2,
                           color = Group))+ 
  # 选择X轴Y轴并映射颜色
  geom_point(size = 4)+ # 画散点图并设置大小
  geom_hline(yintercept = 0,linetype="dashed") + 
  geom_vline(xintercept = 0,linetype="dashed") + 
  theme_bw() + 
  stat_ellipse(level = 0.95)+ 
  # 提取主成分解释度进行绘图
  labs(x = paste('PC1 (', pca_var_data$var[1],'%)', sep = ''),
       y = paste('PC2 (', pca_var_data$var[2],'%)', sep = ''),
       title = "PCA Score Plot") +
  theme(plot.title = element_text(size = 16, face =  "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12,face = "bold", colour = "black",
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12,face = "bold", colour = "black",
                                   vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 12,face = "bold",
                                   colour = "black"),
        legend.title = element_text(size = 14,face = "bold",
                                    colour = "black")) 
p_score_data
ggsave(filename = paste0(out_path,'pca_score.png'),plot = p_score_data)


p_loading_data <- ggplot(pca_loading_data, aes(PC1, PC2))+ 
  # 选择X轴Y轴并映射颜色
  geom_point(size = 4, alpha = 0.4, color = "#3a8db0")+
  geom_hline(yintercept = 0,linetype="dashed") + 
  geom_vline(xintercept = 0,linetype="dashed") + 
  theme_bw() + 
  # 提取主成分解释度进行绘图
  labs(x = "Loading 1",
       y = "Loading 2",
       title = "PCA loading Plot") +
  theme(plot.title = element_text(size = 16, face =  "bold",
                                  hjust = 0.5),
        axis.text.x = element_text(size = 12,face = "bold", colour = "black",
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12,face = "bold", colour = "black",
                                   vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 12,face = "bold",
                                   colour = "black"),
        legend.title = element_text(size = 14,face = "bold",
                                    colour = "black")) 
p_loading_data
ggsave(filename = paste0(out_path,'pca_loading.png'),plot = p_loading_data)

## OPLS-DA
data_plsda <- opls(x = data_for_ana_neg, 
                   y = meta_data_neg[, 'Group'], 
                   predI = NA, 
                   orthoI = NA,  
                   crossvalI = 5,
                   scaleC = "pareto") 
                   
## Score
data_score <- data_plsda@scoreMN %>%  
  as.data.frame() %>%
  mutate(Group = meta_data_neg[,'Group'],
         o1 = data_plsda@orthoScoreMN[,1]) 
data_score$Pixel <- rownames(data_score)
write.csv(data_score, file = paste(out_path,"OPLS-DA_score.csv",sep = "/"))

## Loading
data_loading <- data.frame("t1" = data_plsda@loadingMN,
                           "o1" = data_plsda@orthoLoadingMN[,1])
data_loading$mz <- rownames(data_loading)
write.csv(data_loading, file = paste(out_path,"OPLS-DA_loading.csv",sep = "/"))

## 得分图
data_score_plot <- ggplot(data_score, aes(p1, o1, color = Group)) +
  # 选择X轴Y轴并映射颜色
  geom_point(size = 4)+ 
  geom_hline(yintercept = 0,linetype="dashed") + 
  geom_vline(xintercept = 0,linetype="dashed") + 
  theme_bw() + 
  stat_ellipse(level = 0.95)+ 
  # 提取主成分解释度进行绘图
  labs(x = paste('p1 (', data_plsda@modelDF[1, "R2X"] * 100,'%)', sep = ''),
       y = paste('o1 (', data_plsda@modelDF[2, "R2X"] * 100,'%)', sep = ''),
       title = "OPLS-DA Score Plot") +
  theme(plot.title = element_text(size = 16, face =  "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12,face = "bold", colour = "black",
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12,face = "bold", colour = "black",
                                   vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 12,face = "bold",
                                   colour = "black"),
        legend.title = element_text(size = 14,face = "bold",
                                    colour = "black"))
data_score_plot
ggsave(filename = paste0(out_path,'OPLS-DA_score.png'),plot = data_score_plot)

## 载荷图
data_loading_plot <- ggplot(data_loading, aes(p1, o1)) +
  geom_point(size = 5, alpha = 0.4, color = "#3a8db0")+ 
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_vline(xintercept = 0, linetype="dashed") + 
  theme_bw() + 
  labs(x = 'Loading 1',
       y = 'Loading 2',
       title = "OPLS-DA Loading Plot") +
  theme(plot.title = element_text(size = 16, face =  "bold", hjust = 0.5),
        axis.text.x = element_text(size = 12,face = "bold", colour = "black",
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12,face = "bold", colour = "black",
                                   vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 12,face = "bold",
                                   colour = "black"),
        legend.title = element_text(size = 14,face = "bold",
                                    colour = "black")) 
data_loading_plot
ggsave(filename = paste0(out_path,'OPLS-DA_loading.png'),plot = data_loading_plot)

## VIP
data_VIP <- as.data.frame(getVipVn(data_plsda))
data_VIP$mz <- rownames(data_VIP)
names(data_VIP) [1]<- "VIP"
# 按 VIP 降序排序
data_VIP <- data_VIP[order(data_VIP$VIP, decreasing = TRUE),
                     c(2,1)]
write.csv(data_VIP, file = paste(out_path,"OPLS-DA_VIP.csv",sep = "/"))


## log2FC与T test
data_for_ana_neg$Pixel <- rownames(data_for_ana_neg)
data1 <- merge(meta_data_neg[,2:3],
               data_for_ana_neg,
               by = "Pixel") # 将数据与meta数据整合起来
data2 <- melt(data1,
              variable.name = "mz",
              value.name = "Abundance")
data2$mz <- as.numeric(as.character(data2$mz)) 
Mz <- unique(data2$mz)

## 先建立一个基准的均值数据框
Mz_1 <- Mz[1]
data_1 <- subset(data2, mz == Mz_1)
mean_data <- aggregate(data_1$Abundance,
                       by = list(Group = data_1$Group),
                       mean)
names(mean_data)[2] <- Mz_1 
rownames(mean_data) <- mean_data$Group

## 开始循环计算
for (i in 2:length(Mz)) {
  Mz_i <- Mz[i]
  data_i <- subset(data2, mz == Mz_i)
  mean_data_i <- aggregate(data_i$Abundance,
                           by = list(Group = data_i$Group),
                           mean)
  names(mean_data_i)[2] <- Mz_i
  rownames(mean_data_i) <- mean_data_i$Group
  mean_data <- merge(mean_data, mean_data_i,
                     by = "Group")
}
rownames(mean_data) <- mean_data$Group
mean_data <- mean_data[,-1]

## 需要转置
mean_data_t <- as.data.frame(t(mean_data))
names(mean_data_t) <- c("Average(case)","Average(control)")
## 计算FC
mean_data_t$FoldChange <- mean_data_t$`Average(case)`/mean_data_t$`Average(control)`
## 计算log2FC
mean_data_t$log2FoldChange <- log2(mean_data_t$FoldChange)
mean_data_t$mz <- rownames(mean_data_t)

## 计算p值，T test方法
## 同样先建立一个基准的数据框
data_1_case <- subset(data_1, Group == "case")
data_1_con <- subset(data_1, Group == "control")
p_value <- t.test(data_1_case$Abundance, 
                  data_1_con$Abundance,
                  paired = FALSE,
                  alternative= "two.sided")$p.value
p_data <- data.frame(mz = Mz_1,
                     p_value = p_value)

# 开始循环计算
for (i in 2:length(Mz)) {
  Mz_i <- Mz[i]
  data_i <- subset(data2, mz == Mz_i)
  data_i_case <- subset(data_i, Group == "case")
  data_i_con <- subset(data_i, Group == "control")
  p_value_i <- t.test(data_i_case$Abundance, data_i_con$Abundance,
                    paired = FALSE,
                    alternative= "two.sided")$p.value
  p_data_i <- data.frame(mz = Mz_i,p_value = p_value_i)
  p_data <- rbind(p_data, p_data_i)
}


## p值校正，方法是BH
p_data$q_value <- p.adjust(p_data$p_value, method = "BH")
## 将之前算的FC、log2FC等与p值和校正p值整合在一起
diff_data_neg <- merge(mean_data_t, p_data,by = "mz")
## 读取定性结果
qua <- read.xlsx(paste0(path,"Qualitative.xlsx"),sheet = "neg")
## 匹配定性结果信息
diff_data_neg <- merge(qua, diff_data_neg,by = "mz")
diff_data_neg <- diff_data_neg[,c(2,1,3:length(diff_data_neg))]
diff_data_neg_select <- subset(diff_data_neg, abs(log2FoldChange)>=0 &p_value < 0.05)
## 输出表格
data_out <- list("差异表达矩阵" = diff_data_neg,"差异表达矩阵(未筛选)" =diff_data_neg_select)
write.xlsx(data_out, file = paste0(out_path,"case-vs-control.xlsx"))

