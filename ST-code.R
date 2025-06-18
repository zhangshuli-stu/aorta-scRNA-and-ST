rm(list=ls())
options(stringsAsFactors = F) 

Sys.setenv(GITHUB_PAT = "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")
Sys.getenv("GITHUB_PAT")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glmGamPoi)
library(viridis)
library(harmony)

options(timeout = 600000000)
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
devtools::install_github("dmcable/RCTD", build_vignettes = F)
BiocManager::install("spacexr")
library(spacexr)

getwd()

mycol = c("#A0C2E7","#6894B9","#8798A6","#E0D9E0","#EDBAA7","#FADB7F","#F3B646",
          "#EF9749","#B27466","#646F3F","#899678","#C2BC9A","#868A63","#C4C3BE",
          "#DFA0A6","#98B3D9","#E4BE92","#CB6B7A","#D5CBDA","#f1707d","#f15536",
          "#ef5767","#ae716e","#cb8e85","#cf8878","#c86f67","#f1ccb8","#f2debd",
          "#b8d38f","#ddff95","#ff9b6a","#f1b8f1","#d9b8f1","#f1ccb8","#f1f1b8",
          "#b8f1ed","#e7dbca","#e26538","#f3d751","#fd803a","#fe997b")

spatial_counts1 <- Read10X("./GSM7917897_NC1/filtered_feature_bc_matrix/")
spatial_image1 <- Read10X_Image(
  image.dir = file.path("./GSM7917897_NC1", "spatial"), filter.matrix = TRUE)
spatial_obj_nc1 <- CreateSeuratObject(counts = spatial_counts1, 
                                      project = "GSM7917897_NC1", assay = "Spatial")
# 37487 X 4992
spatial_image1 <- spatial_image1[Cells(x = spatial_counts1)]
DefaultAssay(spatial_counts1 = spatial_image1) <- "Spatial"
spatial_obj_nc1[["slice1"]] <- spatial_image1
SpatialDimPlot(spatial_obj_nc1)
SpatialFeaturePlot(spatial_obj_nc1, features = "nFeature_Spatial")
#各样本分别进行SCT 
spatial_obj_nc1 <- SCTransform(spatial_obj_nc1, assay = "Spatial", verbose = FALSE)
# transformsample_objects<- lapply(sample_objects,
#                                   SCTransform, assay="Spatial", method="poisson")
saveRDS(spatial_obj_nc1, "./spatial_obj_nc1.rds")

spatial_counts2 <- Read10X("./GSM7917898_NC2/filtered_feature_bc_matrix/")
spatial_image2 <- Read10X_Image(
  image.dir = file.path("./GSM7917898_NC2", "spatial"), filter.matrix = TRUE)
spatial_obj_nc2 <- CreateSeuratObject(counts = spatial_counts2, 
                                      project = "GSM7917898_NC2", assay = "Spatial")
# 36601 X 4992
spatial_image2 <- spatial_image2[Cells(x = spatial_counts2)]
DefaultAssay(spatial_counts2 = spatial_image2) <- "Spatial"
spatial_obj_nc2[["slice1"]] <- spatial_image2
SpatialDimPlot(spatial_obj_nc2)
SpatialFeaturePlot(spatial_obj_nc2, features = "nFeature_Spatial")
                  # brewer.pal(9, name = "RdPu"))
spatial_obj_nc2 <- SCTransform(spatial_obj_nc2, assay = "Spatial", verbose = FALSE)
saveRDS(spatial_obj_nc2, "./spatial_obj_nc2.rds")

spatial_counts3 <- Read10X("./GSM7917899_NC3/filtered_feature_bc_matrix/")
spatial_image3 <- Read10X_Image(
  image.dir = file.path("./GSM7917899_NC3", "spatial"), filter.matrix = TRUE)
spatial_obj_nc3 <- CreateSeuratObject(counts = spatial_counts3, 
                                  project = "GSM7917899_NC3", assay = "Spatial")
# 36601 X 4992
spatial_image3 <- spatial_image3[Cells(x = spatial_counts3)]
DefaultAssay(spatial_counts3 = spatial_image3) <- "Spatial"
spatial_obj_nc3[["slice1"]] <- spatial_image3
SpatialDimPlot(spatial_obj_nc3)
SpatialFeaturePlot(spatial_obj_nc3, features = "nFeature_Spatial")
spatial_obj_nc3 <- SCTransform(spatial_obj_nc3, assay = "Spatial", verbose = FALSE)
saveRDS(spatial_obj_nc3, "./spatial_obj_nc3.rds")

spatial_obj_nc1 <- readRDS("./spatial_obj_nc1.rds")
spatial_obj_nc2 <- readRDS("./spatial_obj_nc2.rds")
spatial_obj_nc3 <- readRDS("./spatial_obj_nc3.rds")
sample_objects <- list(spatial_obj_nc1, spatial_obj_nc2, spatial_obj_nc3)
combined_object <- merge(
  x = spatial_obj_nc1,
  y = c(spatial_obj_nc2, spatial_obj_nc3),
  add.cell.ids = c("nc1", "nc2", "nc3"),  # 每个样本加前缀，避免细胞名重复
  project = "Spatial_GSE248608_NC"
)
saveRDS(combined_object, "./combined_nc_object.rds")
#合并所有样本ST<-merge(sample_objects[[1]], y=sample_objects[2:4])

ST <- readRDS("./combined_nc_object.rds")
ST@images[["slice1.GSM7917897_NC1"]] <- ST@images[["slice1"]]
ST@images[["slice1"]] <- NULL
ST@images <- ST@images[c("slice1.GSM7917897_NC1", 
                         "slice1.GSM7917898_NC2", 
                         "slice1.GSM7917899_NC3")]
SpatialDimPlot(ST)
SpatialFeaturePlot(ST, features = "nFeature_Spatial")
VariableFeatures(ST) <- c(VariableFeatures(sample_objects[[1]]),
                          VariableFeatures(sample_objects[[2]]),
                          VariableFeatures(sample_objects[[3]]))
ST <- RunPCA(ST, assay="SCT", verbose=FALSE)

ST <- RunHarmony(ST, reduction = "pca", group.by.vars = "orig.ident",
                 reduction.save = "harmony")
ST <- FindNeighbors(ST, reduction = "harmony", dims = 1:30) 
ST <- FindClusters(ST, resolution = 0.4)
# ST <- RunUMAP(ST, dims = 1:30)
ST <- RunUMAP(ST, reduction = "harmony", dims = 1:30, 
                reduction.name = 'umap', reduction.key = 'UMAP_')

saveRDS(ST, "./ST_nc.rds")
ST <- readRDS("./ST_nc.rds")

display.brewer.all()
colors <- colorRampPalette(brewer.pal(12, 'Paired'))(10)
st_cols = c("#E41A1C", "#4A72A6", "#48A462", "#7E6E85", "#D16948", "#FFB716", "#E1C62F", "#B75F49",
            "#EC83BA", "#999999")
DimPlot(ST, reduction = "umap", cols = st_cols)
ggsave("./ST_UMAP_clusters.pdf", width = 5.5, height = 5)

Idents(ST) <- ST$seurat_clusters
names(st_cols) <- Idents(ST) %>% levels()
SpatialDimPlot(ST, cols = st_cols, pt.size.factor = 2, label.size = 5)
ggsave("./ST_clusters.pdf", width = 15, height = 5)
SpatialDimPlot(ST, label = TRUE, label.size = 3)
ST <- PrepSCTFindMarkers(ST)
ST.markers <- FindAllMarkers(ST, only.pos = F, min.pct = 0.25, 
                             logfc.threshold = 0.25, verbose = T, assay = 'SCT')
st_markers_top100 <- ST.markers %>% group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC)
st_markers_top100 <- unstack(st_markers_top100, gene ~ cluster)
names(st_markers_top100) <- gsub("X", "cluster", names(st_markers_top100))

gene_lists <- split(st_markers_top100$gene, st_markers_top100$cluster)
max_len <- max(sapply(gene_lists, length))
gene_lists_padded <- lapply(gene_lists, function(x) {
  length(x) <- max_len
  return(x)
})
gene_df <- as.data.frame(gene_lists_padded)
colnames(gene_df) <- paste0("cluster_", colnames(gene_df))
write.csv(file="./st_markers_top100.csv", gene_df, row.names=F)

#### COSG: 使用COSG寻找差异基因
pro = 'cosg_celltype_'
library(COSG)
marker_cosg <- cosg(ST, groups='all', assay='SCT',
                    slot='data', mu=1, n_genes_user=100)
save(marker_cosg, file = paste0(pro,'marker_cosg.Rdata'))

### 
### GO: ALL 
x_ALL <- compareCluster(marker_cosg[[1]], fun='enrichGO', OrgDb = 'org.Hs.eg.db', 
                        keyType = 'SYMBOL', ont="BP", pvalueCutoff = 0.05)
cluster_gene_list <- split(ST.markers$gene, ST.markers$cluster)
x_ALL_1 <- compareCluster(cluster_gene_list, fun='enrichGO', OrgDb = 'org.Hs.eg.db', 
                        keyType = 'SYMBOL', ont="BP", pvalueCutoff = 0.05)
save(x_ALL, file = "enrichGO.Rdata")
p <- dotplot(x_ALL, showCategory=5, label_format=55) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 20),
        axis.text.y = element_text(size=18),
        panel.spacing = unit(7, "mm")) + 
  scale_fill_gradientn(colours = c("#980043","#E7E1EF")) 
# "#330066", "#FFCC33","#980043","#C994C7","#E7E1EF"
p
ggsave(p, filename = './GO_cosg_go_all_top5.pdf', width = 12, height = 19)

## gene名字的转换
symbols_list <- marker_cosg[[1]]
symbols_list
x <- lapply(symbols_list, function(y){
  # print(y)
  out <- bitr(y,fromType = 'SYMBOL', 
              toType = 'ENTREZID', OrgDb = org.Hs.eg.db)[, 2]
})
# 找到最长向量长度 
max_len = max(sapply(x, length))
# 用NA补齐其他向量，因为不一定能够所有的ID都能完成转换，会有缺失，如果有缺失就没法
# 将list直接转为data.frame
my_list <- lapply(x, function(x) {
  if(length(x) < max_len){
    c(x, rep(NA, max_len - length(x)))
  }else{
    x
  }
})
# 转换为数据框
gene <- as.data.frame(my_list)
colnames(gene) <- names(x)
# KEGG通路富集
x2 <- compareCluster(gene, fun ="enrichKEGG", organism ='hsa', pvalueCutoff = 0.05)
save(x2, file = "enrichKEGG.Rdata")
KEGG_plot <- dotplot(x2, showCategory=7, label_format=45) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 22),
        axis.text.y = element_text(size=20),
        panel.spacing = unit(7, "mm")) + 
  scale_fill_gradientn(colours = brewer.pal(11, "RdGy"))
KEGG_plot                  
ggsave(KEGG_plot, filename = './cosg_KEGG_top5.pdf', width = 11, height = 15)

res <- x2@compareClusterResult
## 将富集结果中的 ENTREZID 重新转为 SYMBOL， geneID 转换 method2
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res[i,"geneID"] = gene_names
}
# geneID 转换 method2
compare_results <- x2@compareClusterResult
compare_results$geneID_list <- str_split(compare_results$geneID, "/")
# 创建一个函数，将 ENTREZID 转换为 SYMBOL
convert_to_symbol <- function(entrez_ids) {
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,                     # 数据库（替换为适合的物种）
    keys = entrez_ids,                # ENTREZID 列表
    column = "SYMBOL",                # 转换为 SYMBOL
    keytype = "ENTREZID",             # 输入类型是 ENTREZID
    multiVals = "first"               # 如果有多个映射结果，取第一个
  )
  return(symbols)
}
compare_results$geneSYMBOL <- lapply(compare_results$geneID_list, convert_to_symbol)
compare_results$geneSYMBOL <- sapply(compare_results$geneSYMBOL, function(x) paste(na.omit(x), collapse = ","))
compare_results$geneID_list <- NULL
class(compare_results)

sample_names <- unique(ST$orig.ident)
sample_names

st.all <- readRDS("./ST_nc.rds")
st1 <- readRDS("./st_nc1.rds")
st1_1 <- readRDS("./spatial_obj_nc1.rds")
table(st1_1@images[["slice1"]]@boundaries[["centroids"]]@coords == 
        st1@images[["slice1.GSM7917897_NC1"]]@boundaries[["centroids"]]@coords)
ST_split <- list()
for (sample in sample_names) {
  subset_obj <- subset(ST, subset = orig.ident == sample)
  ST_split[[sample]] <- subset_obj
}

saveRDS(ST_split[["GSM7917897_NC1"]], "st_nc1.rds")                                                                                                                                                                                                                                                                                                                                                                                                                                             
saveRDS(ST_split[["GSM7917898_NC2"]], "st_nc2.rds") 
saveRDS(ST_split[["GSM7917899_NC3"]], "st_nc3.rds") 
