
library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
setwd("~/Desktop/model")
# 01-1----
rawdata_path <- './gse124310'
filename <- list.files(rawdata_path)   
rawdata_path <- paste(rawdata_path,filename,sep = '/')
rawdata_path

sceList <- lapply(rawdata_path, function(x){
  obj <- CreateSeuratObject(counts = Read10X(x),
                            project = str_split(x,'/')[[1]][3])
})

names(sceList) <- filename   

sce <- merge(sceList[[1]],sceList[-1],add.cell.ids = names(sceList),project = '###')
sce@meta.data$group <- str_split(sce@meta.data$orig.ident,'_',simplify = T)[,2]     

# 01-2----
# 
grep('^MT',x=rownames(sce@assays$RNA@data),value = T)
grep('^RP[SL]',x=rownames(sce@assays$RNA@data),value = T)
grep('^HB[^(P)]',x=rownames(sce@assays$RNA@data),value = T)
sce <- PercentageFeatureSet(sce,'^MT',col.name = 'percent_MT')      
sce <- PercentageFeatureSet(sce,'^RP[SL]',col.name = 'percent_RP')
sce <- PercentageFeatureSet(sce,'^HB[^(P)]',col.name = 'percent_HB')

VlnPlot(sce,features = "nCount_RNA",pt.size = 0,y.max = 10000)
VlnPlot(sce,features = "nFeature_RNA",pt.size = 0,y.max = 2500)
VlnPlot(sce,features = "percent_MT",pt.size = 0)
VlnPlot(sce,features = "percent_RP",pt.size = 0)
VlnPlot(sce,features = "percent_HB",pt.size = 0,y.max = 0.1)
VlnPlot(sce,features = c("nCount_RNA","nFeature_RNA","percent_MT"),pt.size = 0,group.by = 'orig.ident')

dim(sce) 

# 01-3----
# 
sce <- subset(sce,subset = nCount_RNA>1000 & nFeature_RNA>300 & percent_MT<25)
# 
sce <- sce[rowSums(sce@assays$RNA@counts>0)>3,]
# 
sce <- subset(sce,subset = percent_MT<25 & percent_RP<30 & percent_HB<0.1)  

# ----
s_feature <- cc.genes.updated.2019$s.genes
g2m_feature <- cc.genes.updated.2019$g2m.genes

sce <- CellCycleScoring(sce,
                        s.features = s_feature,
                        g2m.features = g2m_feature,
                        set.ident = T)
VlnPlot(sce,features = c('S.Score','G2M.Score'),group.by = 'orig.ident',pt.size = 0)
saveRDS(sce,'sce_qc.rds')

#----
library(harmony)
library(dplyr)
library(Seurat)
library(clustree)
# 02-1Normalize----
sce <- readRDS('sce_qc.rds')
sce <- NormalizeData(sce,
                     normalization.method = 'LogNormalize',
                     scale.factor = 10000)
sce <- FindVariableFeatures(sce,
                            selection.method = "vst",
                            nfeatures = 2000)


sce <- ScaleData(sce)
sce <- RunPCA(sce,features = VariableFeatures(sce))
DimPlot(sce,reduction = 'pca',group.by = 'group')  

# 02-2harmony----
sce <- RunHarmony(sce,group.by.vars = 'orig.ident')  

ElbowPlot(sce,reduction = 'harmony')      
sce <- RunUMAP(sce,dims = 1:10,reduction = 'harmony')  
sce <- RunTSNE(sce,dims = 1:10,reduction = 'harmony')
DimPlot(sce,reduction = 'umap',label = T,group.by = 'group')
DimPlot(sce,reduction = 'tsne',label = T,group.by = 'group')

sce <- FindNeighbors(sce,reduction = 'harmony',dims = 1:10)  
# 
sce_res <- sce
for (i in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)){   
  sce_res <- FindClusters(sce_res,resolution = i)
}
clustree(sce_res,prefix = 'RNA_snn_res.')   

sce <- FindClusters(sce,resolution = 0.5)  
saveRDS(sce,file = 'step2_harmony.rds')
DimPlot(sce,reduction = 'umap',group.by = 'orig.ident')
DimPlot(sce,reduction = 'umap',group.by = 'seurat_clusters')


# 06-1-----
sparse_data <- as(as.matrix(sce@assays$RNA@counts),'sparseMatrix')  
mdata <- new('AnnotatedDataFrame',data=sce@meta.data)              
fData <- data.frame(gene_short_name=row.names(sparse_data),row.names = row.names(sparse_data))
fd <- new('AnnotatedDataFrame',data=fData)         

monocle_cds <- newCellDataSet(cellData = sparse_data,
                              phenoData = mdata,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

# 06-2----
monocle_cds <- estimateSizeFactors(monocle_cds)  
monocle_cds <- estimateDispersions(monocle_cds)  

monocle_cds <- detectGenes(monocle_cds,min_expr = 0.1)  



sce_var_gene <- VariableFeatures(sce)
monocle_cds <- setOrderingFilter(monocle_cds,sce_var_gene)
monocle_cds <- reduceDimension(monocle_cds,num_dim=10,norm_method = 'log',reduction_method = 'tSNE')
monocle_cds <- clusterCells(monocle_cds,num_clusters = 10)
plot_cell_clusters(monocle_cds,color_by = 'singleR_label')

diff_test_gene <- differentialGeneTest(monocle_cds[sce_var_gene,],fullModelFormulaStr = '~singleR_label')
diff_gene <- row.names(subset(diff_test_gene,qval<0.01))
monocle_cds <- setOrderingFilter(monocle_cds,diff_gene)   

monocle_cds <- reduceDimension(monocle_cds,reduction_method = 'DDRTree')

# 06-4----
source('order_cells.R')       
library('igraph')
my_ordercell <- function(cds){
  root_state = NULL
  num_paths = NULL
  reverse = NULL
  root_cell <- select_root_cell(cds, root_state, reverse)
  cds@auxOrderingData <- new.env(hash = TRUE)
  
  if (cds@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering[row.names(pData(cds)), ]$pseudo_time
    K_old <- reducedDimK(cds)
    old_dp <- cellPairwiseDistances(cds)
    old_mst <- minSpanningTree(cds)
    old_A <- reducedDimA(cds)
    old_W <- reducedDimW(cds)
    cds <- project2MST(cds, project_point_to_line_segment)
    minSpanningTree(cds) <- cds@auxOrderingData[[cds@dim_reduce_type]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(V(old_mst)$name == root_cell, arr.ind = T)
    cells_mapped_to_graph_root <- which(cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- V(minSpanningTree(cds))[cells_mapped_to_graph_root]$name
    tip_leaves <- names(which(degree(minSpanningTree(cds)) == 1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- select_root_cell(cds, root_state, reverse)
    }
    cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- extract_ddrtree_ordering(cds, root_cell)
    pData(cds)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(cds)), ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      pData(cds)$State <- cc_ordering[closest_vertex[, 1], ]$cell_state
    }
  }
  cds
}
monocle_cds <- my_ordercell(monocle_cds)   

# 06-5----
plot_cell_trajectory(monocle_cds,color_by = 'Pseudotime')
plot_cell_trajectory(monocle_cds,color_by = 'State')
plot_cell_trajectory(monocle_cds,color_by = 'singleR_label')

plot_cell_trajectory(monocle_cds,color_by = 'singleR_label')+facet_wrap(~singleR_label,nrow = 3)  

saveRDS(monocle_cds,file = 'monocle.rds')

