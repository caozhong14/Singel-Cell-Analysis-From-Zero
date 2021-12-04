#!/usr/bin/env Rscript

# . 依赖关系检查、安装和加载
# See whether these packages exist on comp. If not, install.
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("MAST")
package_list <- c("optparse", "Seurat","ggplot2","patchwork", "dplyr", "showtext", "ggprism", "rstatix", "clusterProfiler", "ggpubr", "ggrepel", "DOSE", "ggnewscale")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    # install.packages(p, repos="http://cran.r-project.org")
    install.packages(p, repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

if  (!requireNamespace("remotes",  quietly  =  TRUE))  {
 install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

# 另两种常见R包安装方法
# if (FALSE){
#  # Bioconductor安装
#  source("https://bioconductor.org/biocLite.R")
#  biocLite(c("reshape2"))
  # Github安装
#  install.packages("devtools", repo="http://cran.us.r-project.org")
#  library(devtools)
#  install_github("kassambara/ggpubr")
#}
