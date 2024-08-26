library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(spacexr)

InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")

#Data preprocessing
png(filename="Violin_And_Spatial_Feature_Plot.png")
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
dev.off()

options(future.globals.maxSize = 1024 * 1024 * 1024) 
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE) #Normalization
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)


png(filename="Dimplot_And_Spatial_Dimplot.png")
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2
dev.off()

png(filename="Spatial_Dimplot_Cells_By_Identities.png")
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c(1, 6, 13)), facet.highlight = TRUE)
dev.off()


#Integration with a scRNA-seq reference
ref <- readRDS("mouse_hippocampus_reference.rds")
ref <- UpdateSeuratObject(ref)

##Predicting major cell types
anchors <- FindTransferAnchors(reference = ref, query = slide.seq, normalization.method = "SCT", npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE, weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay

#Visualization of the prediction scores of the major expected classes
DefaultAssay(slide.seq) <- "predictions"
png(filename="Prediction_scores_major_expected_classes.png")
SpatialFeaturePlot(slide.seq, features = c("Dentate Principal cells", "CA3 Principal cells", "Entorhinal cortex", "Endothelial tip", "Ependymal", "Oligodendrocyte"), alpha = c(0.1, 1))
dev.off()

slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
png(filename="SpatialDimPlot_Cells_By_Identities_Integrated_scRNA.png")
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c("CA3 Principal cells", "Dentate Principal cells", "Endothelial tip")), facet.highlight = TRUE)
dev.off()

#Identification of Spatially Variable Features
#DefaultAssay(slide.seq) <- "SCT"
#slide.seq <- FindSpatiallyVariableFeatures(slide.seq, assay = "SCT", slot = "scale.data", features = VariableFeatures(slide.seq)[1:1000], selection.method = "moransi", x.cuts = 100, y.cuts = 100)
#png(filename="SpatialFeaturePlot_Top6features.png")
#SpatialFeaturePlot(slide.seq, features = head(SpatiallyVariableFeatures(slide.seq, selection.method = "moransi"), 6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")
#dev.off()

#Spatial deconvolution using RCTD
Idents(ref) <- "celltype"

# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$celltype)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA
slide.seq <- SeuratData::LoadData("ssHippo")
counts <- slide.seq[["Spatial"]]$counts
coords <- GetTissueCoordinates(slide.seq)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 60)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
slide.seq <- AddMetaData(slide.seq, metadata = RCTD@results$results_df)

#Plotting RCTD Annotations
png(filename="RCTD_Annotation_Plot.png")
p1 <- SpatialDimPlot(slide.seq, group.by = "first_type")
p2 <- SpatialDimPlot(slide.seq, group.by = "second_type")
p1 | p2
dev.off()

