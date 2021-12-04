#!/usr/bin/env Rscript


# Intro: Gene set

## Clear workspace
#########{r packages, cache=FALSE}
rm(list=ls())
set.seed(123)
WORKDIR = '~/Singel-Cell-Analysis-From-Zero/Deal_With_Gene_Matrix/'
setwd(WORKDIR)


## BiocManager::install("clusterProfiler")
## BiocManager::install("org.Hs.eg.db")
## BiocManager::install("MAST")

library(patchwork)
library(dplyr)
library(showtext)
library(ggprism)
library(rstatix)
library(optparse)
library(Seurat)
library(gplots)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db) ##加载人类
library(Seurat)

library(fgsea)
library(dplyr)
library(tibble)
library(Seurat)
library(enrichplot)
library(ggpubr)
library(ggrepel)
library(DOSE)
library(ggnewscale)
library(showtext)
#########


if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="../output/pbmc_ASS_ILD_HD_cluster.rds",
                help="Input pbmc file to read [default %default]"),
    make_option(c("-o", "--output"), type="character", default="../Figure2_DEG_0.3/",
                help="output directory or prefix [default %default]"),
    make_option(c("-c", "--celltype"), type="character", default="celltypeI",
                help="Input pbmc file to read [default %default]"),
    make_option(c("-l", "--log2fc"), type="double", default=0.3,
                help="Input pbmc file to read [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))
  

}

# 解析命令行
TEMPDIR <- paste(opts$output, stringr::str_split(tempdir(), '/', simplify=T)[1,3], sep="/")
if (!dir.exists(TEMPDIR)){
  dir.create(path = TEMPDIR, recursive = TRUE)
} 

# 显示输入输出确认是否正确
print(paste("The input file is ", opts$input,  sep = ""))
print(paste("The input log2fC is ", opts$log2fc,  sep = ""))
print(TEMPDIR)

# 3. 读取输入文件
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 从给定的文件读取
if (FALSE){
  # pbmc <- readRDS("../output/pbmc_ASS_ILD_HD_orig.rds")
  pbmc <- readRDS("../output/ASSHD_pbmc_cluster_com.rds")
  DefaultAssay(pbmc) <- "RNA"
  summary(pbmc)
  Idents(pbmc) <- opts$celltype
  pbmc <- NormalizeData(pbmc, verbose = FALSE)
}

# 从文件中读取
if (TRUE){
  pbmc <- readRDS(opts$input)  
  DefaultAssay(pbmc) <- "RNA"
  Idents(pbmc) <- opts$celltype
  pbmc <- NormalizeData(pbmc, verbose = FALSE)
}

# 弹出窗口选择文件
if (FALSE){
  pbmc <- read.table(file.choose())
  DefaultAssay(pbmc) <- "RNA"
  Idents(pbmc) <- opts$celltype
  pbmc <- NormalizeData(pbmc, verbose = FALSE)  
}



print ("Load success")



## Section 1 Loading data ----
#########{r dataset}

showtext_auto(TRUE)
font_add("Times New Roman", regular = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf", 
         italic = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Italic.ttf",
         bold = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold_Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold.ttf")

treat = stringr::str_split(colnames(pbmc),"1|2|3|4|5|6|7|8|9",simplify = T)[,1]
# treat <- gsub('T', 'TREAT', treat)
# treat <- gsub('C', 'CTRL', treat)
table(treat)
pbmc@meta.data$treat <- treat

#########
#########{r theme}
theme_me <- theme(text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  axis.text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  plot.caption = element_text(family="Times New Roman", size=10.5, face = "italic"),
                  plot.title = element_text(family="Times New Roman", size=12, face = "bold"),
                  plot.subtitle = element_text(family="Times New Roman", size=11, face = "bold"),
                  legend.title = element_text(family="Times New Roman", size=11, face = "plain"),
                  strip.text = element_text(family="Times New Roman", size=10.5, face = "italic"))

theme_plain <- theme(text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  axis.text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  plot.caption = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  plot.title = element_text(family="Times New Roman", size=12, face = "bold"),
                  plot.subtitle = element_text(family="Times New Roman", size=11, face = "bold"),
                  legend.title = element_text(family="Times New Roman", size=11, face = "plain"),
                  strip.text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  strip.background = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"))

theme_set(theme_plain)
# theme_set(theme_me)
#########

## Find markers
#########{r function for markers}
funMarkers <- function(celltype){
  pbmc_cell <- subset(pbmc, idents = celltype)
  DefaultAssay(pbmc_cell) <- "RNA"
  Idents(pbmc_cell) <- "treat"
  # markers <- FindMarkers(pbmc_cell, ident.1 = "TREAT", ident.2 = "CTRL", min.pct = 0.1)
  markers <- FindMarkers(pbmc_cell, ident.1 = "ASS-ILD", ident.2 = "HD", min.pct = 0.1)
  
  PATH = paste(TEMPDIR, celltype, sep="/")
  if (!dir.exists(PATH)){
    dir.create(path = PATH, recursive = TRUE)
  }
  filename = paste(PATH, "markers.csv", sep="/")
  write.table(markers, file = filename)  
  return(markers)
}
#########



## Volcano plots
#########{r function for geneset}
volcano_plot <- function(markers){
  # add a column of NAs
  P_VALUE = 0.05
  LOG2FC = opts$log2fc
  de <- markers
  de$gene_symbol <- rownames(markers)
  de$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de$diffexpressed[de$avg_log2FC > LOG2FC & de$p_val_adj < P_VALUE] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$diffexpressed[de$avg_log2FC < -LOG2FC & de$p_val_adj < P_VALUE] <- "DOWN"
  mycolors <- c("blue", "red", "grey")
  names(mycolors) <- c("DOWN", "UP", "NO")
  
  de <- de[order(de$p_val_adj), ]
  up10 <- head(de$gene_symbol[which(de$diffexpressed == "UP")], 10)
  down10 <- head(de$gene_symbol[which(de$diffexpressed == "DOWN")], 10)
  top10.gene <- c(as.character(up10), as.character(down10))
  
  de$delabel = NA
  de$delabel[match(top10.gene, de$gene_symbol)] <- top10.gene  
  
  volplot <- ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
    geom_point() + 
    theme_plain +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "grey", "red")) +
    geom_vline(xintercept=c(- LOG2FC,  LOG2FC), col="black") +
    geom_hline(yintercept=-log10(P_VALUE), col="black")

  PATH = paste(TEMPDIR, celltype, sep="/")
  if (!dir.exists(PATH)){
    dir.create(path = PATH, recursive = TRUE)
  }
  ggsave("volcano.pdf", volplot, path = PATH, width = 7, height = 7)  
  return (volplot)
}  
#########



## GO ontology Analysis
#########{r function for Go ontology Analysis}
## script for go analysis
scriptGo <- function(celltype, markers){
  print ("GO processing  ...")
  genedic <- funFindGene(markers)
  ego_up <- funGoEnrich(genedic$up)
  ego_down <- funGoEnrich(genedic$down)
  ego_all <- funGoEnrich(genedic$all)
  ## Save csv files
  PATH = paste(TEMPDIR, celltype, "Figure2_Geneset", "GO", sep="/")
  if (!dir.exists(PATH)){
    dir.create(path = PATH, recursive = TRUE)
  }
  write.csv(ego_up, file =paste(PATH, "/ego_up.csv", sep="")) 
  write.csv(ego_down, file =paste(PATH, "/ego_down.csv", sep="")) 
  write.csv(ego_all, file =paste(PATH, "/ego_all.csv", sep=""))   
  
  # plot
  p1 <- barplot(ego_up, showCategory = 20)
  p2 <- barplot(ego_down, showCategory = 20)
  p3 <- barplot(ego_all, showCategory = 20)
  p4 <- dotplot(ego_up, showCategory = 20)
  p5 <- dotplot(ego_down, showCategory = 20)
  p6 <- dotplot(ego_all, showCategory = 20)

  ## Save figures
  ggsave("barplot_up.pdf", p1, path = PATH, width = 15, height = 7)
  ggsave("barplot_down.pdf", p2, path = PATH, width = 15, height = 7)
  ggsave("barplot_all.pdf", p3, path = PATH, width = 15, height = 7)
  ggsave("dotplot_up.pdf", p4, path = PATH, width = 15, height = 7)
  ggsave("dotplot_down.pdf", p5, path = PATH, width = 15, height = 7)
  ggsave("dotplot_all.pdf", p6, path = PATH, width = 15, height = 7)

  ## GSEA
  PATH = paste(TEMPDIR, celltype, "Figure3_GSEA", "GO", sep="/")
  if (!dir.exists(PATH)){
    dir.create(path = PATH, recursive = TRUE)
  }  
  try ({
    geneList <- funRenameGeneList(markers)
    gsea <- funGoGSEA(geneList)
    write.csv(gsea, file =paste(PATH, "/go_gsea.csv", sep=""))   

    p7 <- dotplot(gsea, showCategory = 20)
    p8 <- cnetplot(gsea, foldChange = geneList)
    p9 <- heatplot(gsea, foldChange = geneList)
    p10 <- gseaplot2(gsea, geneSetID = 1:3)  

    ggsave("dotplot_gsea.pdf", p7, path = PATH, width = 15, height = 7)
    ggsave("cnetplot_.pdf", p8, path = PATH, width = 15, height = 7)
    ggsave("headplot.pdf", p9, path = PATH, width = 15, height = 7)
    ggsave("gseaplot.pdf", p10, path = PATH, width = 15, height = 7)    
  })

  
}


scriptKEGG <- function(celltype, markers){
  print ("KEGG processing ...")
  genedic <- funFindGene(markers)
  kk_up <- funKEGGEnrich(genedic$up)
  kk_down <- funKEGGEnrich(genedic$down)
  kk_all <- funKEGGEnrich(genedic$all)
  ## Save csv files
  PATH = paste(TEMPDIR, celltype, "Figure2_Geneset", "KEGG", sep="/")
  if (!dir.exists(PATH)){
    dir.create(path = PATH, recursive = TRUE)
  }  
  write.csv(kk_up, file =paste(PATH, "/kk_up.csv", sep="")) 
  write.csv(kk_down, file =paste(PATH, "/kk_down.csv", sep="")) 
  write.csv(kk_all, file =paste(PATH, "/kk_all.csv", sep="")) 

  # plot
  p1 <- barplot(kk_up, showCategory = 20)
  p2 <- barplot(kk_down, showCategory = 20)
  p3 <- barplot(kk_all, showCategory = 20)
  p4 <- dotplot(kk_up, showCategory = 20)
  p5 <- dotplot(kk_down, showCategory = 20)
  p6 <- dotplot(kk_all, showCategory = 20)
  
  ## Save figures
  ggsave("barplot_up.pdf", p1, path = PATH, width = 15, height = 7)
  ggsave("barplot_down.pdf", p2, path = PATH, width = 15, height = 7)
  ggsave("barplot_all.pdf", p3, path = PATH, width = 15, height = 7)
  ggsave("dotplot_up.pdf", p4, path = PATH, width = 15, height = 7)
  ggsave("dotplot_down.pdf", p5, path = PATH, width = 15, height = 7)
  ggsave("dotplot_all.pdf", p6, path = PATH, width = 15, height = 7)


  # GSEA
  PATH = paste(TEMPDIR, celltype, "Figure3_GSEA", "KEGG", sep="/")
  if (!dir.exists(PATH)){
    dir.create(path = PATH, recursive = TRUE)
  }  

  try ({
    geneList <- funRenameGeneList(markers)
    gsea <- funKEGGGSEA(geneList)
    write.csv(gsea, file =paste(PATH, "/go_gsea.csv", sep=""))   

    p7 <- dotplot(gsea, showCategory = 20)
    p8 <- cnetplot(gsea, foldChange = geneList)
    p9 <- heatplot(gsea, foldChange = geneList)
    p10 <- gseaplot2(gsea, geneSetID = 1:3)    
    
    ggsave("dotplot_gsea.pdf", p7, path = PATH, width = 15, height = 7)
    ggsave("cnetplot_.pdf", p8, path = PATH, width = 15, height = 7)
    ggsave("headplot.pdf", p9, path = PATH, width = 15, height = 7)
    ggsave("gseaplot.pdf", p10, path = PATH, width = 15, height = 7)    
  })

  
}

### Find Gene 
funFindGene <- function(markers){
  P_VALUE = 0.05
  LOG2FC = opts$log2fc
  up <-rownames(markers[intersect(which(markers [,1] < P_VALUE),which(markers [,2]>= LOG2FC)),]) # p < 0.05, avg_log2FC > 1
  down <-rownames(markers[intersect(which(markers [,1] < P_VALUE),which(markers [,2]<=(- LOG2FC))),]) # p < 0.05, avg_log2FC < 1
  all <- rownames(markers[intersect(which(markers [,1] < P_VALUE),which(abs(markers [,2]) >= LOG2FC)),]) # p < 0.05, |avg_log2FC| > 1
  out <- list(up=up, down=down, all=all)  
  return(out)
}


### Enrichment analysis
funGoEnrich <- function(gene){
  gene.df <- bitr(gene, fromType="SYMBOL", 
          toType="ENTREZID", 
          OrgDb="org.Hs.eg.db")
  head(gene.df)
  
  ggo <- groupGO(gene     = gene.df$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

  head(ggo)
  
  ## GO over-representation test
  ego <- enrichGO(gene          = gene.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                keyType       = "ENTREZID",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
  head(ego)
  return (ego)
}

funKEGGEnrich <- function(gene){
  gene.df <- bitr(gene, fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Hs.eg.db")
  head(gene.df)
  ## GO over-representation test
  kk <- enrichKEGG(gene     = gene.df$ENTREZID, 
                  organism  = 'hsa',
                  pvalueCutoff = 0.05)
  head(kk)
  return (kk)
}

### Rename gene list
funRenameGeneList <- function(markers){
  markers <- markers[which(markers [,2]>-Inf), ]
  markers <- markers[which(markers [,2]<Inf), ]
  # print (markers)
  nrDEG = markers[, c('avg_log2FC', 'p_val')]
  colnames(nrDEG)=c('log2FoldChange','pvalue') ##更改列名
  gene <- bitr(rownames(nrDEG),     
               fromType = "SYMBOL",     
               toType =  "ENTREZID",    
               OrgDb = org.Hs.eg.db)  
  ## 基因名、ENTREZID、logFC一一对应起来
  gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))]
  ## 构建 gene list
  geneList=gene$logFC
  names(geneList)=gene$ENTREZID 
  geneList=sort(geneList,decreasing = T) # 降序，按照logFC的值来排序
  geneList <- geneList[!duplicated(names(geneList))]
  # print(geneList)
  return (geneList)
}

## GSEA
funGoGSEA <- function(geneList){
  ego <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                nPermSimple  = 1000, 
                minGSSize    = 10,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                verbose      = FALSE,
                by = 'fgsea',
                seed = 123)  
  gsea = setReadable(ego, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  return (gsea)
}

funKEGGGSEA <- function(geneList){
  kk <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPermSimple  = 1000, 
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose      = FALSE,
               by = 'fgsea',
               seed = 123)  
  gsea = setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  return (gsea)
}

## 
#########

## Test
#########{r cellnames}
cellnames <- levels(pbmc)
# cellnames <- c("DC", "Platelet")
for (celltype in cellnames) {
  if (celltype == "doublet" ) next
  if (celltype == "Platelet" ) next
  print ("deal with celltype:")
  print(celltype)
  markers <- funMarkers(celltype)
  try(volcano_plot(markers))
  try(scriptGo(celltype, markers))
  try(scriptKEGG(celltype, markers))
}  
#########