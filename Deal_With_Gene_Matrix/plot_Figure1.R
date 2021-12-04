#!/usr/bin/env Rscript

# 1. 程序功能描述和主要步骤

# 程序功能 聚类结果与marker基因展示
# Script functions: show umap cluster and plot markers 
# Main steps: 
# - Reads data input.rds
# - Draw dim feature, vln feature and save in output.pdf

# 程序使用示例
# USAGE
# Default
# plot_Fugure1.R   -i input pbmc.rds
#                  -c input cell type
#                  -o otuput filename prefix for output directory name 

# 参数说明
# Options
# -i/--input    输入数据 pbmc.rds
# -o/--output   输出结果文件名前缀 output_prefix, 通常会有统计表txt和矢量图pdf
options(warn = -1)


# 2. 依赖关系检查、安装和加载
# See whether these packages exist on comp. If not, install.
package_list <- c("optparse", "Seurat","ggplot2","patchwork", "dplyr", "showtext", "ggprism", "rstatix")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 另两种常见R包安装方法
if (FALSE){
  # Bioconductor安装
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("reshape2"))
  # Github安装
  install.packages("devtools", repo="http://cran.us.r-project.org")
  library(devtools)
  install_github("kassambara/ggpubr")
}

# 清理工作环境 clean enviroment object
rm(list=ls()) 
set.seed(123)
WORKDIR = '~/Singel-Cell-Analysis-From-Zero/Deal_With_Gene_Matrix/'
setwd(WORKDIR)

# 加载依赖关系 Load essential packages
library(optparse)
library(Seurat)
# library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(showtext)
library(ggprism)
library(rstatix)

# 解析命令行
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="../output/pbmc_ASS_ILD_HD_cluster.rds",
                help="Input pbmc file to read [default %default]"),
    make_option(c("-c", "--celltype"), type="character", default="celltypeI",
                help="Input pbmc file to read [default %default]"),
    make_option(c("-o", "--output"), type="character", default="../Figure1_Cluster_I/",
                help="output directory or prefix [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  
  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("The input celltype is ", opts$celltype,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
}

# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 从给定的文件读取
if (FALSE){
  # pbmc <- readRDS("../output/pbmc_ASS_ILD_HD_orig.rds")
  pbmc <- readRDS("../output/pbmc_ASS_ILD_HD_cluster.rds")
}

# 从文件中读取
if (TRUE){
  pbmc <- readRDS(opts$input)  
}

# 弹出窗口选择文件
if (FALSE){
  pbmc <- read.table(file.choose())
}

DefaultAssay(pbmc) <- 'RNA'
pbmc <- NormalizeData(pbmc, verbose = FALSE)
Idents(pbmc) <- opts$celltype

treat = stringr::str_split(colnames(pbmc),"1|2|3|4|5|6|7|8|9",simplify = T)[,1]
# treat <- gsub('ASS-ILD', 'ASS', treat)
# treat <- gsub('HD', 'HD', treat)
table(treat)
pbmc@meta.data$treat <- treat

PATH = opts$output
if (!dir.exists(PATH)){
  dir.create(path = PATH, recursive = TRUE)
} 

print ("Load success")



# 4. 统计与绘图
## Plot setting
my16colors <- c("#D93424", "#5186B5", "#60B345", "#895095", "#EA7E1C", 
                "#E7E086", "#9F5425", "#A0CBE1", "#959597", "#5BB89D", 
                "#E74641", "#8F9CC6", "#DA84B4", "#9FCA5A", "#EBAC6A", "#84B0CA")
r16colors <- c("#FF5AA3", "#FF4ED5", "#FF5DF4", "#D674FF", "#7F95FC", 
               "#008DBE", "#00BCE6", "#00C2C7", "#00C598", "#00C25D", 
               "#00B000", "#6DAF00", "#AAA400", "#D29400", "#EF8300", "#FF6F6A", "#EF8300")
colors = c(my16colors, r16colors)

theme_me <- theme(text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  axis.text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  plot.caption = element_text(family="Times New Roman", size=10.5, face = "italic"),
                  plot.title = element_text(family="Times New Roman", size=12, face = "bold"),
                  plot.subtitle = element_text(family="Times New Roman", size=11, face = "bold"),
                  legend.title = element_text(family="Times New Roman", size=11, face = "plain"),
                  strip.text = element_text(family="Times New Roman", size=10.5, face = "italic"))

theme_set(theme_me)

# levels(pbmc) <- cluster.names
if (TRUE) {
  # vln.features <- c("CD3D", "CCR7", "CD8A", "TRDV2", "FOXP3", "MKI67",
  #                   "KLRF1", "NCAM1",
  #                   "CD79A", "MS4A1", "IGHA1", 
  #                   "CD14", "FCGR3A", "MS4A7", 
  #                   "CD1C", "FCER1A", "LILRA4", "CLEC4C", "PPBP", "GP9")
  # vln.features <- c("CD3D", "IL7R", "MAL", "CCR7", "CD8A", "CD8B", "TRDV2", "FOXP3", "MKI67",
  #                   "KLRF1", "KLRD1", "NCAM1",
  #                   "CD79A", "CD79B", "MS4A1", "MZB1", "IGHA1", "IGHA2", "IGLC1", 
  #                   "CD14", "FCGR3A", "MS4A7", "LST1", 
  #                   "CD1C", "FCER1A", "CLEC10A", "LILRA4", "CLEC4C", "PLD4",
  #                   "PPBP", "GP9", "TUBB", "PF4")
  # vln.features <- c("CD3D", "CD3G", "CCR7", "LEF1", "CD4", "CD8A", "CD8B", 
  #                   "GZMA", "GZMK", "GNLY", "PRF1", "TRDC", "MKI67", "RRM2", 
  #                   "KLRF1", "NKG7", "FCGR3A", 
  #                   "CD14", "MS4A7", "LYZ", "S100A8", 
  #                   "CD79A", "MS4A1", "CD1C", "CLEC10A", "LILRA4", "CLEC4C", "PPBP", "GP9")
  # feat.features <- c("CD3D", "CD3G", "CD4", "CD8A", "CD8B", "CCR7", "SELL", "LEF1", "GZMA", "GZMK", "GNLY", "PRF1", "TRDC", "TRDV2", "TRGV9", "TRAV1-2", "SLC4A10", "MKI67", "KLRF", "NCAM1", "FCGR3A", "CD14", "LYZ", "S100A8", "MS4A7", "CD1C", "LILRA4", "CLEC4C", "CD79A", "MZB1", "MS4A1", "PPBP", "GP9")
  vln.features <- c("CD3D", #"CD4", "CD8A", "CD8B", 
                      "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "LTB", "S100A4", "S100A11", "MAL", # Naive/TCM
                      "GPR183", "CCL5",  # CD4 TCM
                      "LINC02446", # CD8 Naive
                      "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "GNLY", "FGFBP2", "TRBV20-1", "TNFRSF4", # Effector
                      "FCGR3A", "NKG7", # NK
                      "TYROBP", "GATA3", "TBX21", "RORC", #
                      "TRGV9", "TRDV2", "TRDC", "CMC1", # gdT
                      "FOXP3", "CTLA4", # Treg
                      "TRAV1-2", "SLC4A10", # MAIT
                      "MKI67", "RRM2", # ProT
                      "CD14", "MS4A7", "LYZ", "S100A8", "S100A9", "LGALS2", "FCN1", "MS4A6A", "IFI30", # Mono
                      "CD19", "MS4A1", "IGHD", "IGHM", "FCER2", "CD19", "CD27", "MZB1", "TCL1A", # B
                      "FCER1A", "CD1C", "CLEC10A", "PLD4", "LILRA4", "CLEC4C", "PPBP", "GP9" # DC, pDC, Platelet
                      )

  feat.features <- vln.features
  if (opts$celltype == "celltypeI"){
    # vln.features <- c("CD3D", "CCR7", "LEF1", "CD8A", "TRDC", 
    #                   "CD14", "MS4A7", "S100A8", 
    #                   "SPON2", "FCGR3A", 
    #                   "CD79A", "MS4A1", "MZB1",
    #                   "CD1C", "CLEC10A", 
    #                   "LILRA4", "CLEC4C", 
    #                   "PPBP", "PF4")
    # feat.features <- c("CD3D", "CCR7", "LEF1", "CD8A", "TRDC", 
    #                   "CD14", "MS4A7", "S100A8", 
    #                   "SPON2", "FCGR3A", 
    #                   "CD79A", "MS4A1", "MZB1",
    #                   "CD1C", "CLEC10A", 
    #                   "LILRA4", "CLEC4C", 
    #                   "PPBP", "PF4")
    vln.features <- c("CD3D", #"CD4", "CD8A", "CD8B", 
                        "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "LTB", "S100A4",  # Naive/TCM
                        "GPR183",
                        "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "GNLY", # Effector
                        "MKI67", "RRM2", # ProT
                        "FCGR3A", "NKG7", # NK
                        "CD14", "MS4A7", "LYZ", "S100A8", "S100A9", # Mono
                        "CD19", "MS4A1",  # B
                        "CD1C", "CLEC10A", 
                        "LILRA4", "CLEC4C", 
                        "PPBP", "GP9" # DC, pDC, Platelet
                        )
    feat.features <- vln.features
    print("To be finished!")
  }
  if (opts$celltype == "celltypeII_T"){
    ## vln.features <- c("CD3D", "CD3G", "CD4", "CD8A",
    ##                   "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "LTB",
    ##                   "TSHZ2", "GNLY", "CD40LG", "FOXP3",
    ##                   "LINC02446", "LRRN3", "CCL5", "GZMK",
    ##                   "GZMH", "GZMA", "TOP2A", "MKI67", "TRGV9",
    ##                   "TRDV2", "FXYD2", "AC004585.1", "CD1D", "KLRF1", "NCAM1")
    ## feat.features <- c("CD3D", "CD3G", "CD4", "CD8A",
    ##                   "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "LTB",
    ##                   "TSHZ2", "GNLY", "CD40LG", "FOXP3",
    ##                   "LINC02446", "LRRN3", "CCL5", "GZMK",
    ##                   "GZMH", "GZMA", "TOP2A", "MKI67", "TRGV9",
    ##                   "TRDV2", "FXYD2", "AC004585.1", "CD1D", "KLRF1", "NCAM1")
    vln.features <- c("CD3D", "CD4", "CD8A",
                      "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "LTB", "S100A4", "S100A11", "MAL",
                      "GPR183",
                      "LINC02446", 
                      "GZMA", "GZMB", "GZMH", "GZMK", "GNLY", "FGFBP2", "TRBV20-1", "TNFRSF4", 
                      "FCGR3A", "NKG7", "TYROBP", 
                      "TRGV9", "TRDV2", # gdT
                      "FOXP3", "CTLA4", # Treg
                      "TRAV1-2", "SLC4A10", # MAIT
                      "MKI67", "RRM2" # ProT
                      )
    feat.features <- vln.features
  }
  # if (opts$celltype == "celltypeII_B"){
  #   vln.features <- c("CD79A", "MS4A1", "CXCR4", "IGHD",
  #                     "AIM2", "CD19", "CD27", "CD38", "MZB1", "CD74", "CD1D")
  #   feat.features <- c("CD79A", "MS4A1", "CXCR4", "IGHD",
  #                     "AIM2", "CD19", "CD27", "CD38", "MZB1", "CD74", "CD1D")
  # }
  if (opts$celltype == "celltypeII_B") {
    vln.features <- c("CD79A", "CD19", "CD27", "CD38", 
                      "TCL1A", "FCER2", "FCRL1", "CCR7", 
                      "IGHM", "IGHD", "IGHA1", "IGHG2", "IGKC", "IGLC1", "AIM2", 
                      "MZB1", "NEIL1") #, "FGR", "TNFRSF1B", "FCRL5", "IFIT3", "IFI6")
    if (opts$input == "../output/ASSHD_pbmc_subcluster_B.rds") {
      vln.features <- c("CD79A", "CD19", "CD27", "CD38", 
                        "TCL1A", "FCER2", "FCRL1", "CCR7", 
                        "IGHM", "IGHD", "IGHA1", "IGHG2", "IGKC", "IGLC1", "AIM2", 
                        "MZB1", "NEIL1", "FGR", "TNFRSF1B", "FCRL5", "IFIT3", "IFI6")      
    }
    feat.features <- vln.features
  }
} 



if (TRUE){
  count_table <- table(Idents(pbmc))
  print(count_table)
  PATH = opts$output
  write.table(count_table, paste(PATH, "cellcountbefore.csv", sep = ""))
  # cluster.names <-  c("CD4+ T", "CD8+ T", "gdT", "Treg", "Pro T", "uc T", "NK", "NKb", "B", "Plasma", "CD14+ Mono", "CD16+ Mono", "DC", "pDC", "Platelet", "doublet")
  cluster.names <- levels(x = pbmc) 
  print ("count and subset")
  if ("doublet" %in% levels(pbmc)){
    pbmc <- subset(pbmc, idents =  cluster.names[cluster.names != "doublet"])
    count_new <- table(Idents(pbmc))
    print(count_new)
    write.table(count_new, paste(PATH, "cellcount.csv", sep = ""))
  } 
  print ("count and subset")
}

vln.features <- unique(vln.features)
feat.features <- unique(feat.features)

if (TRUE){
  plot.cluster.primary <- DimPlot(pbmc, reduction = "umap", label = "FALSE", label.size = 3, repel = TRUE, cols = colors) + theme_me
  plot.cluster.primary.patient <- DimPlot(pbmc, reduction = "umap", label = "FALSE", label.size = 3, repel = TRUE, cols = colors, split.by = "orig.ident") + theme_me
  
  plot.cluster.markers.vln <- VlnPlot(pbmc, vln.features, stack = TRUE, fill.by = "ident", assay = 'RNA', cols = colors, flip = TRUE) + NoLegend() + theme_me #, adjust default or =1
  
  levels(pbmc) <- rev(cluster.names)
  plot.cluster.markers.vln.rev <- VlnPlot(pbmc, vln.features, stack = TRUE, fill.by = "ident", assay = 'RNA', cols = rev(colors[1:length(levels(pbmc))]), flip = TRUE) + NoLegend() + theme_me # adjust default or =1
  plot.cluster.markers.dot.rev <- DotPlot(pbmc, features = vln.features, split.by = "treat", assay = 'RNA', cols = c("blue", "lightgrey"), dot.scale = 5) + RotatedAxis() + NoLegend() + theme_me 
  levels(pbmc) <- cluster.names
  
  plots <- FeaturePlot(pbmc, reduction = "umap", feat.features, combine = FALSE)
  plots <- lapply(X = plots, FUN = function(x) {
    p <- x + theme_me
  })
  plot.cluster.markers.feat <- CombinePlots(plots = plots)
}

if (TRUE){
  showtext_auto(TRUE)
  font_add("Times New Roman", regular = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf", 
           italic = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Italic.ttf",
           bold = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold_Italic.ttf",
           bolditalic = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold.ttf")
  
  
  PATH = opts$output
  print(paste("The output PATH is ", PATH, sep = ""))
  ggsave("cluster_primary.pdf", plot.cluster.primary, path = PATH,  width = 7, height = 5)
  print("The output figure is cluster.primary.pdf",  sep = "")
  ggsave("cluster_primary_patient.pdf", plot.cluster.primary.patient, path = PATH,  width = 16, height = 5)
  print("The output figure is cluster.primary.patient.pdf",  sep = "")
  ggsave("cluster_markers_vln.pdf", plot.cluster.markers.vln, path = PATH, width = 7, height = 7)
  print("The output figure is cluster.markers.vln.pdf",  sep = "")
  ggsave("cluster_markers_vln_rev.pdf", plot.cluster.markers.vln.rev, path = PATH, width = 7, height = 7)
  print("The output figure is cluster.markers.vln.rev.pdf",  sep = "")
  ggsave("cluster_markers_dot.rev.pdf", plot.cluster.markers.dot.rev, path = PATH, width = 7, height = 7)
  print("The output figure is cluster.markers.dot.rev.pdf",  sep = "")  
  # ggsave("cluster_markers_feat.pdf", plot.cluster.markers.feat, path = PATH, width = 16, height = 10)
  # print("The output figure is cluster.markers.feat.pdf",  sep = "")  
  
  filename = paste(PATH, "/cluster_markers_feat.png",sep="") # cluster_markers_feat.pdf is large for publication, this png is for check
  png(filename, width = 3000, height = 2000) 
  print(plot.cluster.markers.feat)
  dev.off()
}

if (TRUE) {
  theme_plain <- theme(text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                       axis.text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                       plot.caption = element_text(family="Times New Roman", size=10.5, face = "plain"),
                       plot.title = element_text(family="Times New Roman", size=12, face = "bold"),
                       plot.subtitle = element_text(family="Times New Roman", size=11, face = "italic"),
                       legend.title = element_text(family="Times New Roman", size=11, face = "plain"),
                       strip.text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                       strip.background = element_blank(),
                       panel.background = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"))
  
  theme_set(theme_plain)
  count_each <- table(Idents(pbmc), pbmc$orig.ident)
  write.table(count_each, paste(PATH, "cellcounteach.csv", sep = ""))
  cell.patient.prop <- as.data.frame(prop.table(table(Idents(pbmc), pbmc$orig.ident), 2))
  colnames(cell.patient.prop) <- c('cluster', 'patient', 'proportion')
  plot.patient.prop <- ggplot(cell.patient.prop, aes(patient, proportion, fill=cluster)) + geom_bar(stat="identity", position="fill") + ggtitle("") + theme_plain + scale_fill_manual(values = colors) + guides(fill=guide_legend(title=NULL)) + theme(axis.text.x = element_text(angle = 45))
  
  
  cell.treat.prop <- as.data.frame(prop.table(table(Idents(pbmc), pbmc$treat), 2))
  colnames(cell.treat.prop) <- c('cluster', 'treat', 'proportion')
  plot.treat.prop <- ggplot(cell.treat.prop, aes(treat, proportion, fill=cluster)) + geom_bar(stat="identity", position="fill") + ggtitle("") + theme_plain + scale_fill_manual(values = colors) + guides(fill=guide_legend(title=NULL)) + NoLegend()
  
  cell.treat.box <- cell.patient.prop
  cell.treat.box$patient <- gsub ("1|2|3|4|5|6|7|8|9", "", cell.treat.box$patient)
  cell.treat.box$percentage <- 100 * cell.treat.box$proportion
  write.table(cell.treat.box, paste(PATH, "cellcountprop.csv", sep = ""))
  df_p_val <- cell.treat.box %>%
    rstatix::group_by(cluster) %>%
    rstatix::t_test(percentage ~ patient) %>%
    rstatix::add_xy_position()
  plot.treat.box <- ggplot(cell.treat.box, aes(patient, percentage, fill=cluster)) + geom_boxplot() + geom_point(position = position_jitterdodge()) + facet_wrap(~cluster, scales = "free", nrow = 2) + add_pvalue(df_p_val, label = "p = {p}", remove.bracket = TRUE, tip.length = 0, y.position = 0) +  ggtitle("") + theme_plain + scale_fill_manual(values = colors) + guides(fill=guide_legend(title=NULL)) + NoLegend() 
}


if (TRUE){
  showtext_auto(TRUE)
  font_add("Times New Roman", regular = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")
  
  PATH = opts$output
  print(paste("The output PATH is ", PATH, sep = ""))
  ggsave("count_patient_proportion.pdf", plot.patient.prop, path = PATH,  width = 7, height = 7)
  print("The output figure is count_patient_proportion.pdf")
  ggsave("count_treat_proportion.pdf", plot.treat.prop, path = PATH,  width = 3, height = 7)
  print("The output figure is count_treat_proportion.pdf")  
  ggsave("count_treat_box.pdf", plot.treat.box, path = PATH,  width = 14, height = 5)
  print("The output figure is count_treat_box.pdf")    
}

print ("Process success")
print ("***************")
print ("")

