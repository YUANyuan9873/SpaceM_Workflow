## 数据路径
path <- "./Data/"
## 结果输出路径设置
out_path <- paste0(path,"5 Pathwayenrich/")
if (dir.exists(out_path)==FALSE) {
  dir.create(out_path, recursive = TRUE)
} else {
  print("it has been already exited")}


# 用筛选的
data2 <- read.xlsx(xlsxFile = paste0(path,"case-vs-control-select.xlsx"),sheet = 1)
# 只用mz，log2FC、p值等
data <- data2[,c(2,12,28,30)] 
rownames(data) <- data$mz
data <- data[,-1]

## 删掉无KEGG ID的物质
enrich_data <- data %>% filter(!is.na(data$KEGG))
head(enrich_data)
## 删掉KEGG列只含有“;”以及“; ;”的行
enrich_data1 <- enrich_data %>% filter(!trimws(KEGG) %in% ";")
enrich_data1 <- enrich_data1 %>% filter(!trimws(KEGG) %in% "; ;")

## 拆分
enrich_data1 <-  enrich_data1 %>% separate_rows(KEGG, sep = "; ")
enrich_data1$ano <- enrich_data1$KEGG # 这个ano实际上就是KEGG ID，只不过后面要用到，为了与“KEGG”区分，所以新建了一列

## 读入通路与描述的对应关系
KEGG_des <- read.xlsx(paste0(path, "background.xlsx"))[,1:2]
names(KEGG_des)[1] <- "Pathway"
KEGG_des <- KEGG_des[!duplicated(KEGG_des), ]

## 读入背景数据并整理
bg_KEGG <- read.xlsx(paste0(path, "gene_kegg_backgroud.xlsx"))
bg_KEGG <- bg_KEGG[,-3]
# 分别处理通路和功能列
bg_KEGG <-  bg_KEGG %>% separate_rows(Pathway, sep = ",")

# 富集得分
# enrichscore = (m/n)/(M/N)
# N为具有该数据库注释的背景代谢物数目；
# n为具有该数据库注释的前景代谢物数目；
# M为某特定条目的背景代谢物数目；
# m为注释为某特定条目的前景代谢物数目。
N_kegg <- length(unique(bg_KEGG$KEGG))
n_KEGG <- length(unique(enrich_data1$KEGG))

enrich_df1 <- merge(bg_KEGG, enrich_data1,
                    by = "KEGG",
                    all = TRUE) 
head(enrich_df1)

## 有的pathway是na，需要删除
enrich_df1 <- enrich_df1 %>% filter(!is.na(enrich_df1$Pathway))
Pathway <- unique(enrich_df1$Pathway)

enrich_result <- data.frame()
# 通过循环来算每一个通路的M和m
for (i in 1:length(Pathway)) { 
  g <- Pathway[i]
  # 挑选特定通路的数据
  a <- subset(enrich_df1, Pathway == g) 
  # 计算M值，某特定条目的背景代谢物数目
  M_i <- length(unique(a$KEGG)) 
  enrich_df2 <- a %>% filter(!is.na(a$ano))
  # 计算m值，注释为某特定条目的前景代谢物数目
  m_i <- length(unique(enrich_df2$KEGG)) 
  # 计算特定通路的富集得分
  enrich_score_i <- (m_i/n_KEGG)/(M_i/N_kegg) 
  # 把特定通路富集到的代谢物拼起来
  DEM <- paste(sort(unique(enrich_df2$ano)), collapse = ";") 
  # 生成富集结果数据框
  enrich_result_i <- data.frame(Pathway = g,
                                m = m_i,
                                M = M_i,
                                n_KEGG = n_KEGG,
                                N_KEGG = N_kegg,
                                enrichment_score = enrich_score_i,
                                Metabolites_ID = DEM) 
  enrich_result <- rbind(enrich_result,enrich_result_i)
}

## 富集p值：超几何检验
enrich_result <- enrich_result %>%
  mutate(
    p_value = phyper(
      q = m - 1, 
      m = M, 
      n = N_KEGG - M, 
      k = n_KEGG, 
      lower.tail = FALSE 
    )
  ) %>%
  mutate(
    adj_p = p.adjust(p_value, method = "BH"))

## 删除enrich_result为0的，也就是没富集到代谢物的通路
enrich_result <- enrich_result %>% filter(enrich_result$m != 0) 
enrich_result <- merge(KEGG_des[,1:2],
                       enrich_result,
                       by = "Pathway")
## 排个序，p从小到大
enrich_result <- enrich_result[order(enrich_result$p_value, decreasing = FALSE),]
head(enrich_result)
write.xlsx(enrich_result,file = paste0(out_path, "enrichment_result.xlsx"))

## 可视化
## 气泡图
bubble_plot_top10 <- ggplot(enrich_result[1:10,], 
                            aes(x = enrichment_score,
                                y = Annotation, 
                                size = m,
                                color = p_value))+
  geom_point(alpha = 0.8)+
  scale_size(range = c(4, 10))+
  scale_colour_gradient(low = "red", high = "green")+
  theme_bw()+
  labs(title = "Pathway Enrichment")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = 16, face =  "bold",hjust = 0.5),
        axis.text.x=element_text(size = 12, colour = "black",
                                 vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14,face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold", colour = "black",
                                    vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 12, face = "bold",
                                   colour = "black"),
        legend.title = element_text(size = 14, face = "bold",
                                    colour = "black"))
bubble_plot_top10
ggsave(filename = paste0(out_path,"富集气泡图_top10.png"),
       plot = bubble_plot_top10)

## 棒棒糖图
lolipop_plot_top10 <- ggplot(enrich_result[1:10,],
                             aes(x = enrichment_score,
                                 y = Annotation)) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 1) + 
  geom_point(aes(color = p_value), size = 4) +         
  geom_bar(aes(fill = "#e67e22"), stat = "identity", width = 0.1) + 
  scale_colour_gradient(low = "red", high = "green")+
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),      
        axis.text.x = element_text(angle = 90),
        legend.position = "None",
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))   
ggsave(filename = paste0(out_path,"富集棒棒糖图_top10.png"),
       plot = lolipop_plot_top10)




