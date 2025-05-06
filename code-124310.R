

setwd("~/MM_TME/raw/GSE124310_RAW_1")
library(readr)
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(scales)
library(cowplot)
library(RCurl)
#非必需包
library(xlsx)
library(DT)

filelist=list.files('./',pattern='_')   #30

sceList = lapply(filelist,function(x){ 
  CreateSeuratObject(counts = Read10X(data.dir = x), 
                     project = x )
})


merged_seurat <- merge(x=sceList[[1]], 
                       y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                              sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]],
                             sceList[[10]],sceList[[11]],sceList[[12]],sceList[[13]],
                             sceList[[14]],sceList[[15]],sceList[[16]],sceList[[17]],sceList[[18]],sceList[[19]],sceList[[20]]
                             ,sceList[[21]],sceList[[22]],sceList[[23]],sceList[[24]],sceList[[25]]
                             ,sceList[[26]],sceList[[27]],sceList[[28]],sceList[[29]],sceList[[30]]), 
                       add.cell.ids = filelist, 
                       project = "GRCh38")

table(merged_seurat@meta.data$orig.ident)


#check
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
save(merged_seurat, file="data/raw_merged_seurat.RData")


# QC ----------------------------------------------------------------------
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA)/log10(merged_seurat$nCount_RNA)

merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata <- merged_seurat@meta.data
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# Create sample column
metadata$sample <- str_split(metadata$cells,"_",simplify=T)[,1]
table(metadata$sample)
#mgus   mm  nbm  smm
#3475 9708 6726 7126 

merged_seurat@meta.data <- metadata

save(merged_seurat, file="data/merged_marked_filtered_seurat.RData")

load("data/merged_marked_filtered_seurat.RData")
metadata <- merged_seurat@meta.data
#  ------------------------------------------------------------------

pdf("QC_fig/NCells.pdf",pointsize = 10)
metadata %>%
  ggplot(aes(x=sample, fill=sample)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

pdf("QC_fig/nUMIs.pdf",pointsize = 10)
metadata %>%
  ggplot(aes(color=sample, x=nUMI, fill= sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = c(500,30000))
dev.off()

pdf("QC_fig/nGene.pdf",pointsize = 10)
metadata %>%
  ggplot(aes(color=sample, x=nGene, fill= sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_x_log10() +
  geom_vline(xintercept = c(200,5000))
dev.off()

pdf("QC_fig/NCells_vs_NGenes.pdf",pointsize = 10)
metadata %>%
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
dev.off()

#Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
pdf("QC_fig/UMIs_vs_genes.pdf",pointsize = 10)
metadata %>%
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) +
  geom_point() +
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  geom_vline(xintercept = 3000) +
  geom_hline(yintercept = 1500) +
  facet_wrap(~sample)
dev.off()

pdf("QC_fig/mitoRatio.pdf",pointsize = 10)
metadata %>%
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  geom_vline(xintercept = 0.15)
dev.off()

pdf("QC_fig/log10GenesPerUMI.pdf",pointsize = 10)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()


# --------------------------------------------------------------------
#
filtered_seurat <- subset(x = merged_seurat,
                          subset= (nUMI >= 500) &
                            (nUMI <= 30000) &
                            (nGene >= 500) &
                            (nGene <= 7000) &
                            (log10GenesPerUMI > 0.80) &
                            (mitoRatio < 0.10))
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

saveRDS(filtered_seurat, file="data/merged_filtered_seurat.rds")


filtered_seurat <- readRDS("data/merged_filtered_seurat.rds")

# harmony流---------------------------------------------------------------

#Cell cycle scoring
# Normalize the counts
harmony_seurat <- NormalizeData(filtered_seurat)

#download.file("https://www.dropbox.com/s/hus4mrkueh1tfpr/cycle.rda?dl=1","cycle.rda")
# Load cell cycle markers
load("cycle.rda")
# Score cells for cell cycle
harmony_seurat <- CellCycleScoring(harmony_seurat,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes) #这里没跑

# Identify the most variable genes
harmony_seurat <- FindVariableFeatures(harmony_seurat,
                                       selection.method = "vst",
                                       nfeatures = 2000,
                                       verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(harmony_seurat)
# Perform PCA
seurat_phase <- RunPCA(seurat_phase)
# Plot the PCA colored by cell cycle phase
pdf("Cell_cycle/Cell_cycle_scoring.pdf",pointsize = 10)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "sample")
dev.off()

pdf("Cell_cycle/sample.pdf",pointsize = 10)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "sample")
dev.off()


harmony_seurat_1 <- ScaleData(harmony_seurat,vars.to.regress = c("mitoRatio"))
harmony_seurat_1 <- RunPCA(harmony_seurat_1)
pdf("QC_fig/before_harmony_pca_1.pdf",height = 5,width = 10)
#options(repr.plot.height = 10, repr.plot.width = 12)
p1 <- DimPlot(object = harmony_seurat_1, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = harmony_seurat_1, features = "PC_1", group.by = "sample", pt.size = .1)
print(plot_grid(p1,p2))
dev.off()


harmony_seurat_2 <- ScaleData(harmony_seurat,vars.to.regress = c("mitoRatio","S.Score","G2M.Score"))
harmony_seurat_2 <- RunPCA(harmony_seurat_2)
pdf("QC_fig/before_harmony_pca_2.pdf",height = 10,width = 12)
#options(repr.plot.height = 10, repr.plot.width = 12)
p1 <- DimPlot(object = harmony_seurat_2, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = harmony_seurat_2, features = "PC_1", group.by = "sample", pt.size = .1)
p3 <- DimPlot(object = harmony_seurat_2, reduction = "pca", pt.size = .1, group.by = "Phase",split.by = "sample")
print(plot_grid(p1,p2,p3))
dev.off()


#harmony
library(harmony)
#without cell cycle
harmony_embeddings <- HarmonyMatrix(
  data_mat  = as.matrix(Embeddings(harmony_seurat_1,"pca")),
  meta_data = harmony_seurat_1@meta.data,
  vars_use  = "sample",
  do_pca = FALSE
)
rownames(harmony_embeddings) <- rownames(Embeddings(harmony_seurat_1,"pca"))
harmony_seurat_1[["harmony"]] <- CreateDimReducObject(embeddings = harmony_embeddings, 
                                                    key = "harmony_", 
                                                    assay = DefaultAssay(harmony_seurat_1))

pdf("QC_fig/after_harmony_pca_1.pdf",height = 5,width = 10)
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = harmony_seurat_1, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = harmony_seurat_1, features = "harmony_1", group.by = "sample", pt.size = .1)
print(plot_grid(p1,p2))
dev.off()

saveRDS(harmony_seurat_1, "data/harmony_seurat_1.rds")  

######
harmony_seurat <- readRDS("data/harmony_seurat_1.rds")

seurat <- RunUMAP(harmony_seurat,reduction = "harmony", dims = 1:20)

# Plot UMAP
pdf("figures/UMAP_harmony_spe_1.pdf",height = 10,width = 12)
DimPlot(seurat,
        pt.size = 1,
        reduction = "umap",
        group.by = "sample",
        split.by = "sample")
dev.off()

pdf("figures/UMAP_harmony_inter_1.pdf",height = 10,width = 12)
DimPlot(seurat,
        pt.size = 1,
        reduction = "umap",
        shuffle=T,
        group.by = "sample")
dev.off()



# Determine the K-nearest neighbor graph
seurat <- FindNeighbors(object = seurat, reduction = "harmony", dims = 1:20)
# -------------------------------------------------------------------------

# Determine the clusters for various resolutions                                
seurat <- FindClusters(object = seurat,
                                  resolution = c(0.2,0.4, 0.6, 0.8, 1.0, 1.2))

library(clustree)
pdf("figures/clustree_2.pdf",width = 20,height = 10)
clustree(seurat,prefix = "RNA_snn_res.") +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
dev.off()



# 
res=0.2

{Idents(object = seurat) <- paste("RNA_snn_res.",res,sep = "")
  # Plot the UMAP
  pdf(paste("figures/cluster_spe_res_",res,".pdf",sep = ""),height = 6,width = 16)
  print(DimPlot(seurat,
                pt.size = 2,
                reduction = "umap",
                shuffle=T,
                label = TRUE,
                label.size = 6,
                split.by = "sample"))
  dev.off()
  pdf(paste("figures/cluster_inter_res_",res,".pdf",sep = ""),height = 6,width = 8)
  print(DimPlot(seurat,
                pt.size = 1.5,
                reduction = "umap",
                shuffle=T,
                label = TRUE,
                label.size = 6))
  dev.off()}
Idents(object = seurat) <- paste("RNA_snn_res.",res,sep = "")

seurat <- RenameIdents(object = seurat, 
                       "0" = "CD4 cell",
                       "1" = "CD14 monocyte",
                       "2" = "CD8 cell",
                       "3" = "NK cell",
                       "4" = "Plasma",
                       "5" = "B cell",
                       "6" = "CD4 cell",
                       "7" = "Hematopoietic stem cell",
                       "8" = "Hematopoietic stem cell",
                       "9" = "CD16 monocyte",
                       "10" = "pDC",
                       "11" = "mDC",
                       "12" = "pre-B cell",
                       "13" = "Hematopoietic progenitor cell",
                       "14" = "Plasmablast")

seurat@meta.data$seurat_clusters <- Idents(object = seurat)
# 
custom_colors <- c("salmon1", "#33a02c", "#e31a1c", "orange", "#6a3d9a",
                   "#a6cee3", "#b2df8a", "#fb9a99", "purple2", "#cab2d6", "deepskyblue4",
                   "#a65628", "royalblue3")

pdf(paste("figures/rename",res,".pdf",sep = ""),height = 6,width = 10)
DimPlot(seurat,
        pt.size = 1.5,
        reduction = "umap",
        shuffle=T,
        label = F,
        label.size = 6,
        cols=custom_colors,
        group.by = "seurat_clusters")
dev.off()


# SCTransform -------------------------------------------------------------
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_seurat, split.by = "sample")

split_seurat <- split_seurat[c("ctrl", "arac")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio","S.Score","G2M.Score"))
}

#
integ_features <- SelectIntegrationFeatures(object.list = split_seurat,
                                            nfeatures = 3000)
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat,
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors,
                                   normalization.method = "SCT")

# Save integrated seurat object
saveRDS(seurat_integrated, "data/integrated_seurat.rds")

seurat_integrated <- readRDS("data/integrated_seurat.rds")

seurat_integrated <- RunPCA(object = seurat_integrated)

pdf("QC_fig/after_SCT_pca.pdf",height = 5,width = 12)
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = seurat_integrated, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = seurat_integrated, features = "PC_1", group.by = "sample", pt.size = .1)
print(plot_grid(p1,p2))
dev.off()

# Run TSNE
seurat_integrated_TSNE <- RunTSNE(seurat_integrated,
                                  dims = 1:20,
                                  reduction = "pca")
# Plot TSNE
pdf("figures/TSNE_SCT_spe.pdf",pointsize = 10)
DimPlot(seurat_integrated_TSNE,
        reduction = "tsne",
        split.by = "sample")
dev.off()

pdf("figures/TSNE_SCT_inter.pdf",pointsize = 10)
DimPlot(seurat_integrated_TSNE,
        reduction = "tsne",)
dev.off()

# Run UMAP
seurat <- RunUMAP(seurat_integrated,
                  dims = 1:20,
                  reduction = "pca")

# Plot UMAP
pdf("figures/UMAP_SCT_spe.pdf",pointsize = 10)
DimPlot(seurat,
        pt.size = 1,
        reduction = "umap",
        group.by = "sample",
        split.by = "sample")
dev.off()

pdf("figures/UMAP_SCT_inter.pdf",pointsize = 10)
DimPlot(seurat,
        reduction = "umap",)
dev.off()

# Determine the K-nearest neighbor graph
seurat <- FindNeighbors(object = seurat, reduction = "pca", dims = 1:20)

# HKU workflow ------------------------------------------------------------
#Normalizing the data
seurat <- NormalizeData(filtered_seurat)
#Identification of highly variable features (feature selection)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat) + 
  theme(legend.position="top")
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + 
  theme(legend.position="none")
pdf("HKU/variable_features.pdf",pointsize = 10)
print(plot1 + plot2)
dev.off()
#Scaling the data
seurat <- ScaleData(seurat, vars.to.regress = "mitoRatio")



# Downstream analysis ----------------------------------------------------
# Explore heatmap of PCs
pdf("figures/heatmap_PCs.pdf",pointsize = 10)
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()
sink(file = "results/5 most variable genes driving PCs.txt")
print(x = seurat_integrated[["pca"]], 
      dims = 1:9, 
      nfeatures = 5)
sink()


#  --------------------------------------------------------------------
# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:20)

# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.2,0.4, 0.6, 0.8, 1.0, 1.2))

library(clustree)
pdf("figures/clustree.pdf",width = 20,height = 10)
clustree(seurat_integrated,prefix = "integrated_snn_res.") +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
dev.off()


res=0.2

{Idents(object = seurat) <- paste("RNA_snn_res.",res,sep = "")
  # Plot the UMAP
  pdf(paste("figures/cluster_spe_res_",res,".pdf",sep = ""),height = 6,width = 16)
  print(DimPlot(seurat,
                pt.size = 2,
                reduction = "umap",
                shuffle=T,
                label = TRUE,
                label.size = 6,
                split.by = "sample"))
  dev.off()
  pdf(paste("figures/cluster_inter_res_",res,".pdf",sep = ""),height = 6,width = 8)
  print(DimPlot(seurat,
                pt.size = 1.5,
                reduction = "umap",
                shuffle=T,
                label = TRUE,
                label.size = 6))
  dev.off()}
Idents(object = seurat) <- paste("RNA_snn_res.",res,sep = "")

seurat <- RenameIdents(object = seurat, 
                       "0" = "CD4 cell",
                       "1" = "CD14 monocyte",
                       "2" = "CD8 cell",
                       "3" = "NK cell",
                       "4" = "Plasma",
                       "5" = "B cell",
                       "6" = "CD4 cell",
                       "7" = "Hematopoietic stem cell",
                       "8" = "Hematopoietic stem cell",
                       "9" = "CD16 monocyte",
                       "10" = "pDC",
                       "11" = "mDC",
                       "12" = "pre-B cell",
                       "13" = "Hematopoietic progenitor cell",
                       "14" = "Plasmablast")

seurat@meta.data$seurat_clusters <- Idents(object = seurat)
# 
custom_colors <- c("salmon1", "#33a02c", "#e31a1c", "orange", "#6a3d9a",
                   "#a6cee3", "#b2df8a", "#fb9a99", "purple2", "#cab2d6", "deepskyblue4",
                   "#a65628", "royalblue3")

pdf(paste("figures/rename",res,".pdf",sep = ""),height = 6,width = 10)
DimPlot(seurat,
        pt.size = 1.5,
        reduction = "umap",
        shuffle=T,
        label = F,
        label.size = 6,
        cols=custom_colors,
        group.by = "seurat_clusters")
dev.off()


