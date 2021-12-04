library(gplots)
library(ggplot2)

input_list <- c("/Users/caozhong/tmp/DEG_gene_analysis_choose/T GO.csv",
                "/Users/caozhong/tmp/DEG_gene_analysis_choose/B GO.csv",
                "/Users/caozhong/tmp/DEG_gene_analysis_choose/CD14 GO.csv",
                "/Users/caozhong/tmp/DEG_gene_analysis_choose/CD16 GO.csv",
                "/Users/caozhong/tmp/DEG_gene_analysis_choose/NK GO.csv",
                "/Users/caozhong/tmp/DEG_gene_analysis_choose/T KEGG.csv",
                "/Users/caozhong/tmp/DEG_gene_analysis_choose/B KEGG.csv")

output_list <- c("T GO ",
                 "B GO ",
                 "CD14 GO ",
                 "CD16 GO ",
                 "NK GO",
                 "T KEGG ",
                 "B KEGG ")

PATH = "TMP"


  
funppplot <- function(input_name, output_name) {
  filename <- input_name
  titlename <- paste(output_name, "Enrichment",  sep = "")
  savename1 <- paste(output_name, "dotplot_up.pdf",  sep = "")
  savename2 <- paste(output_name, "barplot_up.pdf",  sep = "")
  
  ego_up <- read.csv(filename, header = TRUE)
  x <- - log(ego_up$p.adjust)
  y <- factor(ego_up$Description, levels=rev(unique(ego_up$Description)))
  c <- factor(ego_up$Count)
  
  p1 <- ggplot(ego_up, aes(x,y)) +
    geom_point(aes(size=Count,color=p.adjust)) +
    scale_color_gradient(high = "red", low = "blue", trans = 'reverse') +
    labs(color=expression(p.adjust), 
         size="GeneCount",
         x="-LogP",
         y="",
         title=titlename) +
    theme_bw() + 
    scale_size_continuous(range=c(4,8))
  
  
  p2 <- ggplot(ego_up, aes(y, c, fill=p.adjust)) + 
    geom_bar(stat = "identity")  + coord_flip() +
    scale_fill_gradient(high = "red", low = "blue", trans = 'reverse') +
    labs(color=expression(p.adjust), 
         x="",
         y="GeneCount",
         title=titlename) +
    theme_bw()
  
  ggsave(savename1, p1, path = PATH, width = 10, height = 7)
  ggsave(savename2, p2, path = PATH, width = 10, height = 7)  
}

## MAIN
for (i in 1:7){
  print (i)
  funppplot(input_list[i], output_list[i])
}


