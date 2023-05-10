
LXcor <- function (gene_file,meta_file,group1,group2,row_n){
  
#-----------------------------------------------------------------------
R_packs_install <- function() {
  
  R_packs <- c("psych", "pheatmap", "ggplot2", "openxlsx", 
               "ggrepel", "dplyr", "magrittr", "ggplotify","limma")
  list_installed <- installed.packages()
  new_R_packs <- subset(R_packs, !(R_packs %in% list_installed[, 
                                                               "Package"]))
  if (length(new_R_packs) != 0) {
    install.packages(new_R_packs, force = TRUE, quietly = TRUE)
    print(c(new_R_packs, " packages added..."))
  }
  if ((length(new_R_packs) < 1)) {
    print("No new dependency packages added...")
  }
}
R_packs_install()
Bio_packages <- function() {
  Bio_pkgs <- c("WGCNA", "GO.db")
  list_installed <- installed.packages()
  new_pkgs <- subset(Bio_pkgs, !(Bio_pkgs %in% list_installed[, 
                                                              "Package"]))
  if (length(new_pkgs) != 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) 
      install.packages("BiocManager")
    library(BiocManager)
    BiocManager::install(new_pkgs, force = TRUE, quietly = TRUE)
    print(c(new_pkgs, " packages added..."))
  }
  if ((length(new_pkgs) < 1)) {
    print("No new BiocManager packages added...")
  }
}
Bio_packages()
packages <- c("psych", "pheatmap", "ggplot2", "openxlsx", 
              "ggrepel", "dplyr", "magrittr", "WGCNA", "GO.db", "ggplotify")
for (i in packages) {
  library(i, character.only = T)
}
rm(i)


group1 <- trimws(group1)
group2 <- trimws(group2)
group <- data.frame(group1, group2)

genes0 <- read.xlsx(gene_file)

#colnames(genes0)[1] <- "gene_id"

genes_num <- sapply(genes0[,-1],as.numeric,simplify = T)  # 转换成数字
genes0[,c(2:ncol(genes0))] <- genes_num # 替换原表

genes <- limma::avereps(genes0[,-1],genes0[,1]) %>%  # 去重，同时取平均值
        t() # 列、行转置 

#genes <- distinct(genes0, gene_id, .keep_all = TRUE)
#rownames(genes) <- genes$gene_id
#genes <- genes[, -1] %>% t()
gene_n <- ncol(genes)

#------------------------------------------------
meta0 <- read.xlsx(meta_file)

meta_num <- sapply(meta0[,-1],as.numeric,simplify = T)  # 转换成数字
meta0[,c(2:ncol(meta0))] <- meta_num # 替换原表

meta <- limma::avereps(meta0[,-1],meta0[,1]) %>%  # 去重，同时取平均值
  t() # 列、行转置 

#colnames(meta0)[1] <- "meta_id"
#meta <- distinct(meta0, meta_id, .keep_all = TRUE)
#rownames(meta) <- meta$meta_id
#meta <- meta[, -1] %>% t()
meta_n <- ncol(meta)

#------------------------------------------------
r_p <- WGCNA::corAndPvalue(genes, meta, method = "spearman") # 计算相关系数r和p值
r_data <- r_p[["cor"]]
p_data <- r_p[["p"]]
titile_text <- paste("Heatmap graphics", "(", group1, "VS", group2, ")")

if (gene_n > 40 | meta_n > 40) {
  y_name <- rep("", gene_n)
  y_n = round(gene_n/2, 0)
  y_name[y_n] <- "Genes"
  x_name <- rep("", meta_n)
  x_n = round(meta_n/2, 0)
  x_name[x_n] <- "Metabolites"
  p1 <- pheatmap(r_data, 
                 main = titile_text, fontsize = 10, 
                 scale = "row", drop_levels = TRUE, 
                 cluster_rows = TRUE, cluster_cols = TRUE, 
                 show_rownames = T, show_colnames = T, 
                 labels_row = y_name, labels_col = x_name, fontsize_row = 10, 
                 fontsize_col = 10, angle_col = "0",
                 treeheight_col = 50, treeheight_row = 50, border_color = "#fff8dc", 
                 color = colorRampPalette(c("deepskyblue","white", "red"))(100))
} else {
  p1 <- pheatmap(r_data, main = titile_text, fontsize = 10, 
                 scale = "row", drop_levels = TRUE, 
                 cluster_rows = TRUE, cluster_cols = TRUE, 
                 show_rownames = T, show_colnames = T, 
                 fontsize_row = 9, fontsize_col = 9, angle_col = "45", 
                 treeheight_col = 50, treeheight_row = 50, border_color = "#fff8dc", 
                 color = colorRampPalette(c("deepskyblue", "white", "red"))(100))
}

#------------------------------------------------
if (dir.exists("analysis result") == FALSE) 
  dir.create("analysis result")
p1_plot <- as.ggplot(p1)
p1_heat <- p1_plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
p1_heat
ggsave("analysis result/Heatmap graphics 01.png", p1_heat, 
       width = 1200, height = 1000, dpi = 180, units = "px")

if (!is.null(p_data)) {
  sig_01 <- p_data < 0.01
  p_data[sig_01] <- "**"
  sig_05 <- p_data > 0.01 & p_data < 0.05
  p_data[sig_05] <- "*"
  p_data[!sig_01 & !sig_05] <- ""
  }else {
  p_data <- F
  }


if (gene_n > 40 | meta_n > 40) {
  y_name <- rep("", gene_n)
  y_n = round(gene_n/2, 0)
  y_name[y_n] <- "Genes"
  x_name <- rep("", meta_n)
  x_n = round(meta_n/2, 0)
  x_name[x_n] <- "Metabolites"
  p2 <- pheatmap(r_data, main = titile_text, fontsize = 10, 
                 scale = "row", drop_levels = TRUE, cluster_rows = TRUE, 
                 cluster_cols = TRUE, show_rownames = T, show_colnames = T, 
                 labels_row = y_name, labels_col = x_name, fontsize_row = 10, 
                 fontsize_col = 10, angle_col = "0", treeheight_col = 50, 
                 treeheight_row = 50, cutree_rows = row_n, cutree_cols = 1, 
                 border_color = "#fff8dc", 
                 display_numbers = p_data,
                 color = colorRampPalette(c("deepskyblue","white", "red"))(100))
}else {
  p2 <- pheatmap(r_data, main = titile_text, fontsize = 10, 
                 scale = "row", drop_levels = TRUE, cluster_rows = TRUE, 
                 cluster_cols = TRUE, show_rownames = T, show_colnames = T, 
                 treeheight_col = 50, treeheight_row = 50, fontsize_row = 8, 
                 fontsize_col = 8, cutree_rows = row_n, cutree_cols = 1, 
                 border_color = "#fff8dc", angle_col = 45, 
                 display_numbers = p_data, 
                 fontsize_number = 12, color = colorRampPalette(c("deepskyblue","white", "red"))(100))
}

p2_plot <- as.ggplot(p2)
p2_heat <- p2_plot + theme(plot.margin = unit(c(0.5, 0.5, 
                                                0.5, 0.5), "cm"))
p2_heat
ggsave("analysis result/Heatmap graphics 02.png", p2_heat, 
       width = 1200, height = 1000, dpi = 180, units = "px")
print("---------------------------------------------------------------------")
print("The analysis results can be found in the folder of <analysis result>")
row_cluster <- cutree(p2$tree_row, k = row_n)
row_cluster
neworder <- p2$tree_row$order
neworder_df <- r_data[neworder, ] %>% data.frame()
neworder_df$Cluster <- row_cluster[match(rownames(neworder_df), 
                                         names(row_cluster))]
write.xlsx(neworder_df, "analysis result/Cluster_data.xlsx", 
           rowNames = T, colNames = T)

}
