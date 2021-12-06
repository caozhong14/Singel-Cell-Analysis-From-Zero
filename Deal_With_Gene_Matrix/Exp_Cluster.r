#!/usr/bin/env Rscript

# ---
# title: "process_cluster"
# author: "zhong"
# date: "8/12/2021"
# output: html_document
# ---

## Clear workspace
##########{r packages, cache=FALSE}
rm(list=ls())
set.seed(123)
WORKDIR = '~/Singel-Cell-Analysis-From-Zero/Deal_With_Gene_Matrix/' # nolint
setwd(WORKDIR)


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(showtext)  
library(ggprism)
library(rstatix)
##########


showtext_auto(TRUE)
font_add("Times New Roman", regular = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")
PATH = "../Figure1_cluster"
if (!dir.exists(PATH)){
  dir.create(path = PATH, recursive = TRUE)
} 
if (!dir.exists("../output/")){
  dir.create("../output/", recursive = TRUE)
} 

theme_me <- theme(text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  axis.text = element_text(family="Times New Roman", size=10.5, face = "plain"),
                  plot.caption = element_text(family="Times New Roman", size=10.5, face = "italic"),
                  plot.title = element_text(family="Times New Roman", size=12, face = "bold"),
                  plot.subtitle = element_text(family="Times New Roman", size=11, face = "bold"),
                  legend.title = element_text(family="Times New Roman", size=11, face = "plain"),
                  strip.text = element_text(family="Times New Roman", size=10.5, face = "italic"))

theme_set(theme_me)


## Data loading
# To demonstrate mapping to this multimodal reference, we will use a dataset of 2,700 PBMCs generated by 10x Genomics and available via `SeuratData`. 
# Noted that sct transformed can also be replaced with NormalizeData(), ScaleData() and FindVariableFeatures()



######{r tc.load echo=FALSE}
if (!file.exists("../output/pbmc_ASSHD_pbmc3k.orig.rds")){
  DATADIR = '~/Singel-Cell-Analysis-From-Zero/data_filtered_feature_bc_matrix/ASS5HD3'
  WORKDIR = '~/Singel-Cell-Analysis-From-Zero/Deal_With_Gene_Matrix/'
  folders = list.files(DATADIR)
  library(stringr)
  samples=str_split(folders,'PBMC',simplify = T)[,1]

  setwd(DATADIR)
  scrna.list = lapply(X = folders, FUN=function(x){ 
    x <- CreateSeuratObject(counts = Read10X(x), project=x, min.cells = 3, min.features = 200)
  })
  setwd(WORKDIR)

  scrna.all <- merge(scrna.list[[1]], 
                  y = c(scrna.list[[2]], scrna.list[[3]], scrna.list[[4]], scrna.list[[5]],
                        scrna.list[[6]], scrna.list[[7]], scrna.list[[8]]), 
                  add.cell.ids = samples, 
                  project = "LungProj")


  scrna.list <- SplitObject(scrna.all, split.by = "orig.ident")
  rm(scrna.all)
  scrna.list <- lapply(X=scrna.list, FUN = SCTransform)
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = scrna.list, nfeatures=3000)
  scrna.list <- PrepSCTIntegration(object.list = scrna.list, anchor.features = features)
  # Perform integration
  scrna.list <- lapply(X=scrna.list, FUN = RunPCA, verbose = FALSE, features = features)
  pbmc3k.anchors <- FindIntegrationAnchors(object.list = scrna.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca", dim = 1:15)
  # this command creates an 'integrated' data assay
  rm(scrna.list)
  pbmc3k <- IntegrateData(anchorset = pbmc3k.anchors, normalization.method = "SCT", dim = 1:15)

  ######
  #### We then save primary integration of pbmc3k and save it to files.
  ######{r tc.save.orig}
  WORKDIR = '~/Singel-Cell-Analysis-From-Zero/Deal_With_Gene_Matrix/'
  setwd(WORKDIR)
  # saveRDS(pbmc3k, file = "../output/pbmc_TC_orig.rds")
  # saveRDS(pbmc3k, file = "../output/pbmc_ASS_ILD_HD_orig.rds")
  saveRDS(pbmc3k, file = "../output/pbmc_ASSHD_pbmc3k.orig.rds")
  print("Saved data!")
  ######  
}

## Quality Control
##########{r tc.quality.contrl fig.width=10}
if (!file.exists("../output/ASSHD_pbmc_orig.rds")) {
  
  print("Loading data: pbmc_ASSHD_pbmc3k.orig.rds")
  pbmc3k <- readRDS("../output/pbmc_ASSHD_pbmc3k.orig.rds")
  print("Loaded")
  
  DefaultAssay(pbmc3k) <- 'RNA'
  pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")
  plot.qc.vlnplot <- VlnPlot(pbmc3k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot1 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend() + theme_me
  plot2 <- FeatureScatter(pbmc3k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_me
  plot.qc.scatter <- plot1 + plot2

  ##########

  ###### We then choose subset of pbmc3k and save it to files.
  ##########{r tc.save.qc}
  tab1 <- table(pbmc3k$orig.ident)
  write.table(tab1, paste(PATH, "cell_count_wb.csv", sep='/'))
  ##########


  ##########{r with.batchcorr}
  DefaultAssay(pbmc3k) <- 'RNA'
  pbmc3k[["percent.mt"]] <- PercentageFeatureSet(pbmc3k, pattern = "^MT-")
  pbmc <- subset(pbmc3k, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & nCount_RNA < h & percent.mt < 10)

  DefaultAssay(pbmc) <- "integrated"
  pbmc <- ScaleData(object = pbmc, verbose = FALSE)
  pbmc <- RunPCA(object = pbmc, verbose = FALSE)
  pbmc <- RunUMAP(object = pbmc, reduction = "pca", dims = 1:15)
  pbmc <- FindNeighbors(object = pbmc, reduction = "pca", dims = 1:15)
  pbmc <- FindClusters(object = pbmc, resolution = c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4))

  tab2 <- table(pbmc$orig.ident)
  write.table(tab2, paste(PATH, "cell_count_wb_qc.csv", sep='/'))
  plot.pca.elbow <- ElbowPlot(pbmc)

  ######################################

  ######{r tc.save.cluster}
  WORKDIR = '~/Singel-Cell-Analysis-From-Zero/Deal_With_Gene_Matrix/'
  setwd(WORKDIR)

  Idents(pbmc) <- "integrated_snn_res.0.4"
  p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
  Idents(pbmc) <- "integrated_snn_res.0.6"
  p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
  Idents(pbmc) <- "integrated_snn_res.0.8"
  p3 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
  Idents(pbmc) <- "integrated_snn_res.1"
  p4 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
  Idents(pbmc) <- "integrated_snn_res.1.2"
  p5 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
  Idents(pbmc) <- "integrated_snn_res.1.4"
  p6 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
  plot.cluster.resolutions <- p1 + p2 + p3 + p4 + p5 + p6
  ggsave("gene_refmap_integrated.pdf", plot.cluster.resolutions, path = PATH, width = 15, height = 10)
  ggsave("quality_control.pdf", plot.qc.scatter, path = PATH, width = 15, height = 7)
  ggsave("pca_elbow.pdf", plot.pca.elbow, path = PATH, width = 15, height = 7)
  saveRDS(pbmc, file = "../output/ASSHD_pbmc_orig.rds")
} else {
  print("Loading data: ASSHD_pbmc_orig.rds")
  pbmc <- readRDS("../output/ASSHD_pbmc_orig.rds")
  print("Loaded")  
}

##########

#####{r tc.rename.idents.ASSILDHD}
if (TRUE){
  ######{r tc.rename.idents.ASSILDHD}
  WORKDIR = '~/Singel-Cell-Analysis-From-Zero/Deal_With_Gene_Matrix/'
  setwd(WORKDIR)
  Idents(pbmc) <- "integrated_snn_res.0.4"

  pbmc <- RenameIdents(pbmc, 
              "0" = "CD14 Mono", 
              "1" = "TEM",
              "2" = "Naive T",
              "3" = "NK",
              "4" = "TCM",
              "5" = "B", 
              "6" = "CD16 Mono", 
              "7" = "doublet",
              "8" = "gdT",
              "9" = "Platelet",
              "10" = "DC",
              "11" = "Pro T", 
              "12" = "pDC",
              "13" = "doublet")

  cluster.names <-  c("Naive T", "TCM", "TEM", "gdT", "Pro T", "NK", "CD14 Mono", "CD16 Mono", "B", "DC", "pDC", "Platelet", "doublet")
  levels(x = pbmc) <- cluster.names

  plot.cluster.primary <- DimPlot(pbmc, reduction = "umap", label = "TRUE", label.size = 3, repel = TRUE)
  pbmc[["celltypeI"]] <- Idents(pbmc)
  saveRDS(pbmc, file = "../output/ASSHD_pbmc_cluster.rds")
}
