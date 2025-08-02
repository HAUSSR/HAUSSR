library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(motifmatchr)
library(TFBSTools)

load("Alfalfa.sub.Rdata")
pfm <-
  readRDS("Alfalfa.motif.pfmList2.Rds")
motif.family <-
  lapply(pfm, function(x) {
    return(x@tags$family)
  }) %>%
  unlist() %>%
  data.frame(family = ., motif = names(.))
DefaultAssay(Alfalfa.sub) <- "ATAC"
Alfalfa.sub <- AddMotifs(
  object = Alfalfa.sub,
  genome = ZM.1,
  pfm = pfm
)
Alfalfa.sub <- RunChromVAR(
  object = rice.sub,
  genome = ZM.1
)

Alfalfa.sub.atac.markers <-
  parallel::mclapply(unique(Alfalfa.sub$tissue_cluster), function(x) {
    print(x)
    xx <-
      FindMarkers(
        Alfalfa.sub,
        ident.1 = x,
        only.pos = T,
        test.use = "LR",
        logfc.threshold = 0.05,
        max.cells.per.ident = 300,
        latent.vars = "nCount_ATAC"
      )
    
    if (nrow(xx) > 0) {
      return(data.frame(xx, gene = rownames(xx), cluster = x))
    }
  }, mc.cores = 10)

names(Alfalfa.sub.atac.markers) <- unique(Alfalfa.sub$tissue_cluster)
Alfalfa.sub.atac.markers <-
  Alfalfa.sub.atac.markers[!sapply(Alfalfa.sub.atac.markers, is.null)]

get_motif <- function(Alfalfa.sub, cluster = "") {
  top.da.peak <-
    Alfalfa.sub.atac.markers[[cluster]] %>%
    dplyr::filter(p_val < 0.05, avg_log2FC >
                    0.5) %>%
    pull(gene)
  if (length(top.da.peak) < 10) {
    top.da.peak <-
      Alfalfa.sub.atac.markers[[cluster]] %>%
      top_n(300, avg_log2FC) %>%
      pull(gene)
  }
  open.peaks <- AccessiblePeaks(Alfalfae.sub, idents = cluster)
  # match the overall GC content in the peak set
  meta.feature <-
    GetAssayData(Alfalfa.sub, assay = "ATAC", slot = "meta.features")
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[top.da.peak, ],
    n = 50000
  )
  enriched.motifs <- FindMotifs(
    object = rice.sub,
    features = top.da.peak,
    background = peaks.matched
  ) %>% mutate(family = gsub(":.+", "", motif))
  enriched.motifs$cluster <- cluster
  return(enriched.motifs)
}

allmotifs <-
  lapply(names(Alfalfa.sub.atac.markers), function(x) {
    get_motif(Alfalfa.sub, cluster = x)
  })
allmotifs <- do.call(rbind, allmotifs)
DefaultAssay(Alfalfa.sub) <- "chromvar"

DefaultAssay(Alfalfa.sub) <- "ATAC"

pdf("PRR7.CoveragePlot.pdf", 8, 6)
for (xx in unique(Alfalfa.sub$tissue)) {
  p <- CoveragePlot(
    object = Alfalfa.sub,
    region = "MsG0480020869.01",
    # PRR7:MsG0480020869.01
    features = "MsG0480020869.01",
    expression.assay = "RNA",
    extend.upstream = 5000,
    extend.downstream = 5000,
    idents = grep(xx, unique(Alfalfa.sub$tissue_cluster), value = T),
    annotation = T
  )
  print(cowplot::plot_grid(p, labels = xx))
}
dev.off()
