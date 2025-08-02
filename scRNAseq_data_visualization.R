# ========== 1. 加载包 ==========
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

# 设置路径
setwd("/media/desk16/tly9108/riceomics/GSE232863/")
ana.dir <- "./analysis_scRNA_only/"
dir.create(ana.dir, showWarnings = FALSE)
setwd(ana.dir)

# ========== 2. 正确读取 H5 文件，并提取 Gene Expression ==========
file_path <- "../GSM8865409_E06_1_filtered_feature_bc_matrix.h5"

# 读取整个 H5 文件（返回 list）
data_list <- Read10X_h5(filename = file_path)

# 查看内容
names(data_list)  # 应该输出: "Gene Expression" "Peaks"

# ✅ 提取 RNA 矩阵（这才是我们要的 scRNA-seq 数据）
rna_matrix <- data_list[["Gene Expression"]]

# 检查维度
dim(rna_matrix)
# 输出示例：37960 x 20000 → 37960 个基因，20000 个细胞（或核）

# ========== 3. 创建 Seurat 对象（仅基于 RNA）==========
seurat_obj <- CreateSeuratObject(
  counts = rna_matrix,
  project = "rice_E06_1_RNA",
  assay = "RNA",
  min.cells = 3,        # 基因至少在 3 个细胞中表达
  min.features = 200    # 每个细胞至少检测到 200 个基因
)

print(seurat_obj)
# 应该输出：An object of class Seurat 
# 37960 features across 20000 samples

# ========== 4. 添加元数据（QC 指标）==========
# 添加 nCount 和 nFeature
seurat_obj$nCount_RNA <- Matrix::colSums(seurat_obj@assays$RNA@counts)
seurat_obj$nFeature_RNA <- seurat_obj@assays$RNA@data@Dimnames[[1]] |> length()

# （可选）计算线粒体基因比例（水稻需确认命名）
# 常见水稻线粒体基因前缀：如 "mt:" 或 "LOC_OsMt..." 或 "COX", "ATP6" 等
# 先查看一些基因名：
head(x = rownames(seurat_obj), 20)

# 示例：如果线粒体基因包含 "mt:" 或以 "COX" 开头
# 假设你发现有基因如 "mt:COX1", "mt:ATP6" 等
# 则可以：
seurat_obj[["percent_mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt:")

# 如果没有明确标记，可跳过或用其他方式估计

# ========== 5. 初步过滤（根据 QC 图调整）==========
seurat_obj <- subset(seurat_obj, subset = 
                       nFeature_RNA > 200 & 
                       nFeature_RNA < 8000 & 
                       nCount_RNA > 500 &
                       percent_mt < 20
)

print(paste("过滤后细胞数：", ncol(seurat_obj)))

# ========== 6. QC 可视化 ==========
p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)
ggsave(p1, filename = "QC_violin.pdf", width = 12, height = 5)

p2 <- FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA") + ggtitle("nCount vs nFeature")
ggsave(p2, filename = "QC_scatter.pdf", width = 7, height = 5)

# ========== 7. 标准化、高变基因、PCA ==========
seurat_obj <- seurat_obj %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(.), npcs = 30)

# PCA 图
p3 <- DimPlot(seurat_obj, reduction = "pca", label = TRUE)
ggsave(p3, filename = "PCA.pdf", width = 8, height = 6)

# Elbow 图
p4 <- ElbowPlot(seurat_obj, ndims = 30)
ggsave(p4, filename = "elbow_plot.pdf", width = 8, height = 6)

# ========== 8. UMAP + 聚类（仅 RNA）==========
seurat_obj <- seurat_obj %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(dims = 1:20, n.neighbors = 30, min.dist = 0.3)

# UMAP 可视化
p5 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.8) +
  ggtitle("UMAP: RNA-only Clusters")
ggsave(p5, filename = "UMAP_clusters.pdf", width = 9, height = 7)

# 按基因数着色
p6 <- FeaturePlot(seurat_obj, feature = "nFeature_RNA", min.cutoff = "q10")
ggsave(p6, filename = "UMAP_nFeature.pdf", width = 8, height = 7)

# ========== 9. 保存对象 ==========
saveRDS(seurat_obj, file = "seurat_obj_RNA_only.rds")