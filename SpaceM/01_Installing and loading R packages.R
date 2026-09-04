## 定义R包安装函数
metanr_packages <- function(){
  metr_pkgs <- c("Cardinal","openxlsx","ggplot2","ropls","tidyverse","reshape2","pheatmap","ggpubr")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
  
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
  #安装Cardinal指定版本
  remotes::install_version("Cardinal", 
                           version = "2.14.0", 
                           repos = BiocManager::repositories())
  
}
metanr_packages()

## 批量导入R包
metr_pkgs <- c("Cardinal","openxlsx","ggplot2","ropls","tidyverse","reshape2","pheatmap","ggpubr")
sapply(metr_pkgs, library, character.only = T)
