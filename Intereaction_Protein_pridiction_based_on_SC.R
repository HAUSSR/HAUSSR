# ========== 1. 加载必要的包 ==========
library(Seurat)
library(ggplot2)
library(cowplot)  # 提供更好的图排版（可选）

# ========== 2. 设置工作目录（可选）==========
# 可以不设置，直接用完整路径
setwd("/media/desk16/tly9108/riceomics/")

# ========== 3. 加载 RDS 文件 ==========
obj <- readRDS("omics.Rds")

# ========== 4. 检查对象是否包含 UMAP ==========
# 查看有哪些降维方法
names(obj@reductions)
# 正常输出应包含： "pca", "umap" 或类似



p12 <- DimPlot(obj, reduction = "wsnnumap", label = TRUE, pt.size = 0.8)
print(p12)
ggsave(p12, filename = "wsnnumap.pdf", width = 8, height = 7, dpi = 300)

p_gene <- FeaturePlot(
  object = obj,
  features = "Os09g0401200",
  reduction = "wsnnumap",
  pt.size = 0.8,
  min.cutoff = "q10",
)
print(p_gene)

# 保存为 PDF（矢量格式）
ggsave(p_gene, filename = "Os05g0204600_wsnnumap.pdf", width = 8, height = 7, dpi = 300)


# ========== 1. 设置基因名 ==========
gene_of_interest <- "Os05g0204600"

# ========== 2. 检查基因是否在对象中 ==========
if (!gene_of_interest %in% rownames(obj)) {
  stop(paste("Gene", gene_of_interest, "not found in the object."))
}

# ========== 3. 获取 log-normalized 表达数据（'data' 层）==========
expr_matrix <- GetAssayData(obj, assay = "RNA", layer = "data")

# ========== ✅ 新增：提取目标基因的表达向量 ==========
expr_target <- expr_matrix[gene_of_interest, ]

# ========== 4. 筛选高表达细胞（例如 top 25%）==========
high_cells <- names(expr_target)[expr_target >= quantile(expr_target, 0.75)]

# 输出数量
cat("Number of high-expression cells:", length(high_cells), "\n")

# ========== 5. 提取这些细胞的子集表达矩阵 ==========
expr_sub <- expr_matrix[, high_cells, drop = FALSE]

# ========== 6. 计算与其他基因的皮尔逊相关性 ==========
correlations <- apply(expr_sub, 1, function(x) {
  cor(expr_target[high_cells], x, method = "pearson")
})

# ========== 7. 整理结果 ==========
coexp_df <- data.frame(
  gene = names(correlations),
  pearson_cor = correlations,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(gene != gene_of_interest) %>%           # 去掉自己
  dplyr::arrange(desc(pearson_cor))

# ========== 8. 查看和保存结果 ==========
head(coexp_df, 20)

write.csv(coexp_df, "Os05g0204600_coexp_in_high_cells.csv", row.names = FALSE)
