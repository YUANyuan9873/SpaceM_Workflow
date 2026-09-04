## 数据路径
path <- "./Data/"
## 结果输出路径设置
out_path <- paste0(path,"4 Diffvisualization/")
if (dir.exists(out_path)==FALSE) {
  dir.create(out_path, recursive = TRUE)
} else {
  print("it has been already exited")}


## 读入数据
# 要用未筛选的
data <- read.xlsx(xlsxFile = paste0(path,"case-vs-control-select.xlsx"),sheet = 2) 
# 只选择需要的差异数据
data_volcano <- data[,c(2,24,27:28,30:31)] 
rownames(data_volcano) <- data_volcano$mz
data_volcano <- data_volcano[,-1]

# 用筛选的
data2 <- read.xlsx(xlsxFile = paste0(path,"case-vs-control-select.xlsx"),sheet = 1) 
# 除了各离子在样本像素点中的丰度外，留下p值供之后绘图选择离子
data_heatmap <- data2[,c(2,30,32:length(data2))]
rownames(data_heatmap) <- data_heatmap$mz
data_heatmap <- data_heatmap[,-1]
## 创建一个分组和像素对应关系的数据框meta_data
data_heatmap_all <- data_heatmap[,-1]
meta_data <- data.frame(Pixel = colnames(data_heatmap_all),
                        Group = rep(c("Case", "Control"), each = 100))
rownames(meta_data) <- meta_data$Pixel

## 火山图，差异标准为|log2FC|≥0，p＜0.05
names(data_volcano)[4] <- "p"
## 设定差异筛选阈值
log2FC_cutoff <- 0
p_value <- 0.05

## 根据差异标准对离子进行标注“up”和“down”
data_volcano[which(abs(data_volcano$log2FoldChange) < log2FC_cutoff), "sig"] <- "No diff"
data_volcano[which(data_volcano$p > p_value), "sig"] <- "No diff"
data_volcano[which(data_volcano$log2FoldChange >= log2FC_cutoff & 
                     data_volcano$p < p_value), "sig"] <- "Up"
data_volcano[which(data_volcano$log2FoldChange <= -log2FC_cutoff & 
                     data_volcano$p < p_value), "sig"] <- "Down"
head(data_volcano)

## 设置颜色
colors <- data.frame("Up" = "#ff4757", "No diff" = "gray20",
                     "Down" = "#1e90ff",
                     check.names = FALSE)
## 绘制火山图
volcano_plot <- ggplot(data_volcano, aes(log2FoldChange, -log(p, 10))) +
  geom_point(aes(color = sig), 
             alpha = 0.6, 
             size = 2 ) +
  scale_color_manual(values = colors)+
  geom_vline(xintercept = 0, 
             color = "gray", 
             linewidth = 1) +
  geom_hline(yintercept = -log(0.05, 10), 
             color = "gray", 
             linewidth = 1) +
  labs(x = "log2 Fold Change", 
       y = "-log10 p-value",
       color = NA) +
  ggtitle("Case vs Ctrl") +
  theme(title = element_text(size = 16,
                             color = "black", face = "bold",
                             vjust = 0.5, hjust = 0.5, angle = 0),
        panel.grid = element_blank(),
        panel.background = element_rect(color = "black", 
                                        fill = "transparent"),
        legend.title = element_text(face = "bold",
                                    size = 14),
        legend.text = element_text(face = "bold",
                                   size = 14),
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title = "Case vs Ctrl"))
volcano_plot
ggsave(filename = paste0(out_path,'volcano.png'),
       plot = volcano_plot)

## 热图
## 先准备列注释
annotation_col <- as.data.frame(meta_data$Group)
rownames(annotation_col) <- row.names(meta_data)
names(annotation_col) <- "Group"
## 设置下列注释也就是分组的颜色
ann_colors = list(Group = c(Case = "#1abc9c", Control = "#e67e22"))

## 全部差异离子的热图
pheatmap(data_heatmap_all,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         filename = paste0(out_path, "heatmap_all.png"),
         width = 10,
         height = 30) # 看不全的

# p值最小的前50个差异离子热图
data_heatmap <- data_heatmap[order(data_heatmap$`p-value`),]
data_heatmap_top50 <- data_heatmap[1:50,-1]
pheatmap(data_heatmap_top50,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_colors = ann_colors,
         annotation_col = annotation_col,
         filename = paste0(out_path, "heatmap_top50.png"),
         width = 10,
         height = 10)

## 箱型图
data_box <- data2[,c(2,32:length(data2))] 
rownames(data_box) <- data_box$mz
data_box_m <- melt(data_box,
                   id.vars = "mz", 
                   variable.name = "Pixel", 
                   value.name = "abundance") 
## 单离子绘制
data_plot <- subset(data_box_m, mz == 214.05205)
data_plot <- merge(data_plot, meta_data, by = "Pixel")
data_plot$Group <- factor(data_plot$Group,levels = c("Case","Control"))

ggplot(data_plot, aes(x = Group,
                      y = abundance,
                      fill = Group,
                      color = Group))+ 
  geom_boxplot()+
  scale_fill_manual(values = c("#1abc9c", "#e67e22"))+ 
  scale_color_manual(values = c("#1abc9c", "#e67e22")) + 
  stat_compare_means(method = "t.test",
                     paired = FALSE, 
                     method.args = list(alternative = "two.sided"),
                     label = "p.signif", 
                     hide.ns = TRUE) + 
  xlab("Group")+ 
  ylab("abundance")+ 
  ggtitle("mz = 214.05205") +
  theme(
    panel.border = element_blank(), 
    axis.line.x = element_line(color = "black", 
                               linewidth = 0.5), 
    axis.line.y = element_line(color = "black", 
                               linewidth = 0.5), 
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12,
                              hjust = 0.5),
    axis.text.x = element_text(angle = 45, 
                               size = 10, 
                               hjust = 1, 
                               color = "black"),                                                   
    axis.text.y = element_text(size = 10))
ggsave(filename = paste0(out_path, 214.05205,"-neg-箱型图.png"))

## 循环出图
path1 <- paste0(out_path, "循环出图/")
if (!dir.exists(path1)){
  dir.create(path1,recursive = TRUE)
} else {
  print("Dir already exists!")
}
Mz <- unique(data_box_m$mz)
for (i in 1:10) {
  data_plot <- subset(data_box_m, mz == Mz[i])
  data_plot <- merge(data_plot, meta_data, by = "Pixel")
  data_plot$Group <- factor(data_plot$Group,
                            levels = c("Case",
                                       "Control"))
  ggplot(data_plot, aes(x = Group,
                        y = abundance,
                        fill = Group,
                        color = Group))+ 
    geom_boxplot()+
    scale_fill_manual(values = c("#1abc9c", "#e67e22"))+ 
    scale_color_manual(values = c("#1abc9c", "#e67e22")) + 
    stat_compare_means(method = "t.test",
                       paired = FALSE, 
                       method.args = list(alternative = "two.sided"),
                       label = "p.signif", 
                       hide.ns = TRUE) + 
    xlab("Group")+ 
    ylab("abundance")+ 
    ggtitle(Mz[i]) +
    theme(
      panel.border = element_blank(), 
      axis.line.x = element_line(color = "black", 
                                 linewidth = 0.5), 
      axis.line.y = element_line(color = "black", 
                                 linewidth = 0.5), 
      panel.background = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12,
                                hjust = 0.5), 
      axis.text.x = element_text(angle = 45, 
                                 size = 10, 
                                 hjust = 1, 
                                 color = "black"),                                                   
      axis.text.y = element_text(size = 10))
  ggsave(filename = paste0(path1, Mz[i],"-neg-箱型图.png")
  )
}


