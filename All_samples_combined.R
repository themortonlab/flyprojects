rm(list = ls()) 


library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(ggplot2)
library(stringr)
library(sctransform)
#BiocManager::install("glmGamPoi")
library(glmGamPoi)



setwd("/Users/laurynhigginson/Desktop/sn_raw_data/G11A1/")
G11A1 <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
mydata_final <- CreateSeuratObject(counts = G11A1, project = "G11A1")

setwd("/Users/laurynhigginson/Desktop/sn_raw_data/G11A2/")
G11A2 <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
mydata_final <- CreateSeuratObject(counts = G11A2, project = "G11A2")

setwd("/Users/laurynhigginson/Desktop/sn_raw_data/WT1/")
WT1 <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
mydata_final <- CreateSeuratObject(counts = WT1, project = "WT1")

setwd("/Users/laurynhigginson/Desktop/sn_raw_data/WT2/")
WT2 <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
mydata_final <- CreateSeuratObject(counts = WT2, project = "WT2")

setwd("/Users/laurynhigginson/Desktop/sn_raw_data/G146C1/")
G146C1 <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
mydata_final <- CreateSeuratObject(counts = G146C1, project = "G146C1")

setwd("/Users/laurynhigginson/Desktop/sn_raw_data/G146C2/")
G146C2 <- Read10X(data.dir = 'filtered_feature_bc_matrix/')
mydata_final <- CreateSeuratObject(counts = G146C2, project = "G146C2")

counts <- GetAssayData(object = mydata_final, slot = "counts")

# Output a logical vector for every gene on whether there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
mydata_final <- CreateSeuratObject(filtered_counts, meta.data = mydata_final@meta.data)

rm(counts, nonzero, keep_genes, filtered_counts)
gc()



# where gene names are located to find out the ID used for mitochondrial genes
mydata_final@assays$RNA@counts@Dimnames[[1]][grep(pattern = "mt", x = mydata_final@assays$RNA@counts@Dimnames[[1]], fixed = T)]


mydata_final[["percent.mt"]]<-PercentageFeatureSet(mydata_final, pattern = "^mt:")







# orig.ident: this often contains the sample identity if known, but will default to project as we had assigned it
# nCount_RNA: number of UMIs per cell
#UMIs are based on the cell. Each cell gets a UMI and that UMI will get attached to all of the transcripts that are inside that cell.
#This allows us to identify which cell each transcript comes from.
# nFeature_RNA: number of genes detected per cell

# Add number of genes per UMI for each cell to metadata
mydata_final$log10GenesPerUMI <- log10(mydata_final$nFeature_RNA) / log10(mydata_final$nCount_RNA)
metadata<-mydata_final@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
head(metadata$log10GenesPerUMI)
# mean(mydata_final$nCount_RNA)
# sum(mydata_final$nCount_RNA)
#############################################
table(mydata_final$orig.ident)


# orig.ident: this column will contain the sample identity if known. It will default to the value we provided for the project argument when loading in the data
# nCount_RNA: this column represents the number of UMIs per cell
# nFeature_RNA: this column represents the number of genes detected per cell

mydata_final <- subset(mydata_final, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10 & log10GenesPerUMI > .80)

########### Trying the merged method ###############
G11A1<-SCTransform(mydata_final, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)
G11A2<-SCTransform(mydata_final, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)

G146C1<-SCTransform(mydata_final, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)
G146C2<-SCTransform(mydata_final, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)

WT1<-SCTransform(mydata_final, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)
WT2<-SCTransform(mydata_final, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)

G11A1 <- RunPCA(G11A1, features = VariableFeatures(object = G11A1))
G11A2 <- RunPCA(G11A2, features = VariableFeatures(object = G11A2))
G146C1 <- RunPCA(G146C1, features = VariableFeatures(object = G146C1))
G146C2 <- RunPCA(G146C2, features = VariableFeatures(object = G146C2))
WT1 <- RunPCA(WT1, features = VariableFeatures(object = WT1))
WT2 <- RunPCA(WT2, features = VariableFeatures(object = WT2))

G11A1 <- FindNeighbors(G11A1, graph.name = "G11A1", dims = 1:50)
G11A2 <- FindNeighbors(G11A2, graph.name = "G11A2", dims = 1:50)
G146C1 <- FindNeighbors(G146C1, graph.name = "G146C1", dims = 1:50)
G146C2 <- FindNeighbors(G146C2, graph.name = "G146C2", dims = 1:50)
WT1 <- FindNeighbors(WT1, graph.name = "WT1", dims = 1:50)
WT2 <- FindNeighbors(WT2, graph.name = "WT2", dims = 1:50)

G11A1 <- FindClusters(G11A1, resolution = 1, graph.name = "G11A1")
G11A2 <- FindClusters(G11A2, resolution = 1, graph.name = "G11A2")
G146C1 <- FindClusters(G146C1, resolution = 1, graph.name = "G146C1")
G146C2 <- FindClusters(G146C2, resolution = 1, graph.name = "G146C2")
WT1 <- FindClusters(WT1, resolution = 1, graph.name = "WT1")
WT2 <- FindClusters(WT2, resolution = 1, graph.name = "WT2")


G11A1 <- RunUMAP(G11A1, dims = 1:50)
G11A2<- RunUMAP(G11A2, dims = 1:50)
G146C1 <- RunUMAP(G146C1, dims = 1:50)
G146C2<- RunUMAP(G146C2, dims = 1:50)
WT1<- RunUMAP(WT1, dims = 1:50)
WT2<- RunUMAP(WT2, dims = 1:50)

WT1@meta.data$batch<- '1'
G146C1@meta.data$batch<- '1'
G11A1@meta.data$batch <- '1'

WT2@meta.data$batch<- '2'
G146C2@meta.data$batch<- '2'
G11A2@meta.data$batch <- '2'


####Getting FeaturePlots for Arc1####
#G11A combined
G11A<- merge(G11A1, G11A2)
G11A<-SCTransform(G11A, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)

G11A <- ScaleData(G11A, verbose = FALSE)
G11A<- RunPCA(G11A, npcs = 50, verbose = FALSE)
G11A <- RunUMAP(G11A, reduction = "pca", dims = 1:50)

dev.new()
DimPlot(G11A, reduction = "umap", group.by = "orig.ident", shuffle = T, raster = FALSE)
FeaturePlot(G11A, features = 'Arc1', cols = c('azure3', '#06e6e6'))

#G146C combined
G146C<- merge(G146C1, G146C2)
G146C<-SCTransform(G146C, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)

G146C <- ScaleData(G146C, verbose = FALSE)
G146C<- RunPCA(G146C, npcs = 50, verbose = FALSE)
G146C <- RunUMAP(G146C, reduction = "pca", dims = 1:50)

dev.new()
DimPlot(G146C, reduction = "umap", group.by = "orig.ident", shuffle = T, raster = FALSE)
FeaturePlot(G146C, features = 'Arc1', cols = c('azure3', '#06e6e6'))

#WT combined
WT<- merge(WT1, WT2)
WT<-SCTransform(WT, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)

WT <- ScaleData(WT, verbose = FALSE)
WT<- RunPCA(WT, npcs = 50, verbose = FALSE)
WT <- RunUMAP(WT, reduction = "pca", dims = 1:50)

dev.new()
DimPlot(WT, reduction = "umap", group.by = "orig.ident", shuffle = T, raster = FALSE)
FeaturePlot(WT, features = 'Arc1', cols = c('azure3', '#06e6e6'))



#create new group column in the Metadata of each sample
fly_batch_1<- merge(WT1, y = c(G11A1, G11A2, WT2, G146C1, G146C2),
                    project = 'all_samples')
DefaultAssay(fly_batch_1) <- "RNA"
fly_batch_1[['SCT']] <- NULL

setwd('/Users/laurynhigginson/Desktop/')
save(fly_batch_1, file = 'all_samples.RData')
load('all_samples.RData')
fly_batch_1   <- SCTransform(object = fly_batch_1, vars.to.regress = c("nFeature_RNA", "percent.mt"), variable.features.n = 5000)
fly_batch_1 <- RunPCA(fly_batch_1, npcs = 50)
fly_batch_1 <- RunUMAP(fly_batch_1, reduction = 'pca', dims = 1:50)

dev.new()
DimPlot(fly_batch_1, reduction = "umap", group.by = "group", shuffle = T, raster = FALSE)
dev.new()
DimPlot(fly_batch_1, reduction = "umap", group.by = "batch", shuffle = T, raster = FALSE)

fly.list <- SplitObject(fly_batch_1, split.by = "batch")

# normalize and identify variable features for each dataset independently
fly.list <- lapply(X = fly.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = fly.list)

fly.list <- lapply(X = fly.list,
                   FUN = function(x) {
                     x <- ScaleData(x, features = features, verbose = FALSE)
                     x <- RunPCA(x, features = features, verbose = FALSE)
                   })

fly.anchors  <- FindIntegrationAnchors(object.list = fly.list, anchor.features = features, reduction = "rpca", k.anchor = 20)
fly.combined <- IntegrateData(anchorset = fly.anchors)

fly.combined <- ScaleData(fly.combined, verbose = FALSE)
fly.combined <- RunPCA(fly.combined, npcs = 50, verbose = FALSE)
fly.combined <- RunUMAP(fly.combined, reduction = "pca", dims = 1:50)

dev.new()
DimPlot(fly.combined, reduction = "umap", group.by = "batch", shuffle = T, raster = FALSE)
dev.new()
DimPlot(fly_batch_1, reduction = "umap", group.by = "batch", shuffle = T, raster = FALSE)

DefaultAssay(fly.combined) <- "integrated"
fly.combined <- FindNeighbors(fly.combined, dims = 1:50)

fly.combined <- FindClusters(object = fly.combined,
                             resolution = 5,
                             n.start = 100)

#setwd('/Users/laurynhigginson/Desktop/')
#respoint5 only
#save(fly.combined, file = 'all_samples_clean.RData')
#load(file = 'all_samples_clean.RData')

fly.combined <- SetIdent(object = fly.combined, value = 'integrated_snn_res.5')
fly.combined.markers.res1 <- FindAllMarkers(fly.combined, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)


mydata_final_Markers_p_val_adj_0.01<-fly.combined.markers.res1[fly.combined.markers.res1$p_val_adj<=0.01,]

clusters<-c(as.character(unique(mydata_final_Markers_p_val_adj_0.01$cluster)))
mydata_final_top10_markers<-c()
namesofrows<-c()
for (cluster in clusters) {
  namesofrows<-c(namesofrows, paste("Cluster", cluster, sep = "_"))
  mydata_final_top10_markers<-rbind(mydata_final_top10_markers, 
                                    head(mydata_final_Markers_p_val_adj_0.01[mydata_final_Markers_p_val_adj_0.01$cluster==cluster,]$gene, n=10))
}

row.names(mydata_final_top10_markers)<-namesofrows
colnames(mydata_final_top10_markers)<-c("Gene_1", "Gene_2", "Gene_3", "Gene_4", "Gene_5", 
                                        "Gene_6", "Gene_7", "Gene_8", "Gene_9", "Gene_10")

setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/')
write.table(mydata_final_top10_markers, file = "2024May1_G11A_WT_G146C_merge_int_res5_pca50_top10marker_genes.txt", quote = F, sep = "\t", 
            row.names = T, col.names = T, append = F)

dev.new()
sample_colors<- c('WT' = 'azure4', 'G11A' = 'brown', 'G146C' = 'azure4')
sample_colors<- c('WT' = 'blue4', 'G11A' = 'azure4', 'G146C' = 'azure4')
sample_colors<- c('WT' = 'blue4', 'G11A' = 'brown', 'G146C' = 'forestgreen')
DimPlot(object = fly.combined, reduction = "umap", group.by = 'group', cols = sample_colors,raster = FALSE)

#Has res .5 and res5
setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/')
load(file = 'all_samples_clean.RData')

fly.combined[["SCT"]] <- NULL
fly.combined     <- SCTransform(object = fly.combined, vars.to.regress = c("nFeature_RNA", "percent.mt", "batch"), variable.features.n = 5000)
Idents(fly.combined)

dev.new()
VlnPlot(fly.combined, group.by = 'orig.ident', features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster = FALSE)

#res 0.5

fly.combined <- SetIdent(object = fly.combined, value = 'integrated_snn_res.0.5')
fly.combined <- RenameIdents(fly.combined, `0` = "unannotated", `1` = "unannotated", `10` = 'unannotated', `11` = 'unannotated', `12` = 'unannotated', `13` = 'unannotated',
                             `14` = 'Cone cell', `15` = 'Skeletal muscle of head', `16` = 'Pericerebral adult fat mass', `17` = 'Pigment cell', `18` = 'unannotated', `19` = 'Hemocytes', `2` = 'unannotated',
                             `20` = 'Lamina monopolar neuron', `21` = 'unannotated', `22` = 'Ensheathing glial cell', `23` = "unannotated", `24` = "Adult brain perineural glia", `25` = 'unannotated',
                             `26` = 'Columnar neuron T1', `27` = 'Kenyon cell', `28` = 'unannotated', `29` = 'Reticular neuropil associated glia', `3` = 'unannotated', `30` = 'Dopaminergic neuron',
                             `31` = 'Lamina/epithelial marginal glia', `32` = 'T neuron', `33` = 'Centrifugal neuron', `34` = 'unannotated', `35` = 'T neuron',
                             `36` = 'unannotated', `37` = 'unannotated', `38` = 'unannotated', `39` = 'Lamina monopolar neuron', `4` = 'unannotated', `40` = 'Optic lobe associated cortex glia', `41` = 'unannotated', `42` = 'Pericerebral adult fat mass',
                             `5` = 'Photoreceptors', `6` = 'unannotated', `7` = 'unannotated',  `8` = 'Photoreceptors', `9` = 'Photoreceptors')

my_cols <- c( 'unannotated' = 'azure3', 'Photoreceptors' = "red4", 'Skeletal muscle of head' = 'magenta', 'T neuron' = 'brown2', 'Pericerebral adult fat mass' = 'aquamarine', 'Lamina monopolar neuron' = 'cyan4', 'Adult brain perineural glia' = 'plum1',
              'Kenyon cell'= 'royalblue', 'Cone cell' = 'yellow2', 'Ensheathing glial cell'= 'palegreen',  
              'Centrifugal neuron' = 'orchid4', 'Columnar neuron T1' = 'gold2', 'Pigment cell' = 'forestgreen', 'Hemocytes' = 'tomato1',
               'Optic lobe associated cortex glia' = 'indianred', 'Dopaminergic neuron' = 'lightslateblue',
               'Lamina/epithelial marginal glia' = 'lavender', 'Reticular neuropil associated glia' = "salmon2") 


#res 5
fly.combined <- RenameIdents(fly.combined, `0` = "unannotated", `1` = "unannotated", `10` = 'unannotated', `100` = 'Medullary intrinsic neuron Mi1', `101` = 'unannotated', `102` = 'Photoreceptors', `103` = 'unannotated', `104` = 'Skeletal muscle', `105` = 'unannotated',
                             `106` = 'Photoreceptors', `107` = 'unannotated', `108` = 'Johnston organ neuron', `109` = 'Photoreceptors', `11` = 'Photoreceptors', `110` = 'Johnston organ neuron', `12` = 'Pericerebral adult fat mass', `13` = 'unannotated', `14` = 'unannotated', `15` = 'unannotated',
                             `16` = 'unannotated', `17` = 'unannotated', `18` = 'Pericerebral adult fat mass', `19` = 'unannotated', `2` = 'unannotated', `20` = 'unannotated', `21` = 'unannotated', `22` = 'Skeletal muscle', `23` = 'Photoreceptors', `24` = 'Neuron 1', `25` = 'Skeletal muscle', 
                             `26` = 'T neuron T4/5', `27` = 'unannotated', `28` = 'Hemocyte', `29` = 'unannotated', `3` = 'unannotated', `30` = 'Photoreceptors', `31` = 'unannotated', `32` = 'Photoreceptors', `33` = 'Subperineural glia', `34` = 'unannotated', `35` = 'unannotated', `36` = 'unannotated',
                             `37` = 'Adult fat mass', `38` = 'unannotated', `39` = 'unannotated', `4` = 'Pericerebral adult fat mass', `40` = 'unannotated', `41` = 'Columnar neuron T1', `42` = 'T neuron T4/5', `43` = 'Cone cell', `44` = 'Photoreceptors', `45` = 'Lamina monopolar neuron', `46` = 'Photoreceptors', `47` = 'unannotated',
                             `48` = 'Skeletal muscle', `49` = 'Pericerebral adult fat mass', `5`= 'unannotated', `50` = 'Ensheathing glia', `51` = 'unannotated', `52` = 'unannotated', `53` = 'Pigment cell', `54` = 'unannotated', `55` = 'Adult brain cell body glia', `56` = 'unannotated', `57` = 'unannotated',
                             `58` = 'unannotated', `59` = 'Lamina monopolar neuron', `6` = 'unannotated', `60` = 'Kenyon cells', `61` = 'Lamina monopolar neuron', `62` = 'Skeletal muscle', `63` = 'unannotated', `64` = 'Centrifugal neuron', `65` = 'Kenyon cells', `66` = 'Adult optic chiasma glia', `67` = 'Lamina monopolar neuron', 
                             `68` = 'unannotated', `69` = 'Photoreceptors', `7` = 'unannotated', `70` = 'Kenyon cells', `71` = 'Optic-lobe-assoc cortex glia', `72` = 'Pigment cell', `73` = 'Lamina epithelial/marginal glia', `74` = 'unannotated', `75` = 'unannotated', `76` = 'Epithelial cell', 
                             `77` = 'Centrifugal neuron', `78` = 'Ensheathing glia', `79` = 'unannotated', `8` = 'unannotated', `80` = 'Hemocyte', `81` = 'unannotated', `82` = 'unannotated', `83` = 'Columnar neuron T1', `84` = 'unannotated', `85` = 'Adult brain cell body glia', `86` = 'Photoreceptors', `87` = 'Reticular neuropil associated glia',
                             `88` = 'Photoreceptors', `89` = 'Reticular neuropil associated glia', `9` = 'unannotated', `90` = 'Cone cell', `91` = 'Subperineural glia', `92` = 'Ensheathing glia', `93` = 'Reticular neuropil associated glia', `94` = 'unannotated', `95` = 'Pericerebral adult fat mass', `96` = 'Cone cell', 
                             `97` = 'Neuron 2', `98` = 'Transmedullary neuron Tm1', `99` = 'unannotated')
                              
my_cols<- c('unannotated' = 'azure4', 'Medullary intrinsic neuron Mi1' = 'aquamarine', 'Photoreceptors' = 'blue4', 'Kenyon cells' = 'brown2', 'Johnston organ neuron' = 'blueviolet', 'Pericerebral adult fat mass' = 'burlywood2', 'Skeletal muscle' = 'yellow2',
            'Neuron 1' = 'chartreuse1', 'T neuron T4/5' = 'chartreuse4', 'Hemocyte' = 'coral1', 'Subperineural glia' = 'cornflowerblue', 'Adult fat mass' = 'cadetblue','Columnar neuron T1' ='darkmagenta', 'Cone cell' = 'darkolivegreen1', 'Lamina monopolar neuron' = 'darkorange1', 
            'Ensheathing glia' = 'darkslateblue', 'Pigment cell' = 'darkturquoise', 'Centrifugal neuron' = 'deeppink', 'Adult brain cell body glia' = 'deeppink4', 'Adult optic chiasma glia' = 'dodgerblue', 'Optic-lobe-assoc cortex glia' = 'gold2', 'Lamina epithelial/marginal glia' = 'indianred',
            'Epithelial cell' = 'lavender', 'Reticular neuropil associated glia' = 'lightsalmon', 'Neuron 2' = 'palegreen2', 'Transmedullary neuron Tm1' =  'purple1')                             
                             

#res 5 neuronal, glia, head-specific 

fly.combined <- RenameIdents(fly.combined, `0` = "unannotated", `1` = "unannotated", `10` = 'unannotated', `100` = 'Neuron', `101` = 'unannotated', `102` = 'Neuron', `103` = 'unannotated', `104` = 'Head-specific cluster', `105` = 'unannotated',
                             `106` = 'Neuron', `107` = 'unannotated', `108` = 'Neuron', `109` = 'Neuron', `11` = 'Neuron', `110` = 'Neuron', `12` = 'Head-specific cluster', `13` = 'unannotated', `14` = 'unannotated', `15` = 'unannotated',
                             `16` = 'unannotated', `17` = 'unannotated', `18` = 'Head-specific cluster', `19` = 'unannotated', `2` = 'unannotated', `20` = 'unannotated', `21` = 'unannotated', `22` = 'Head-specific cluster', `23` = 'Neuron', `24` = 'Neuron', `25` = 'Head-specific cluster', 
                             `26` = 'Neuron', `27` = 'unannotated', `28` = 'Head-specific cluster', `29` = 'unannotated', `3` = 'unannotated', `30` = 'Neuron', `31` = 'unannotated', `32` = 'Neuron', `33` = 'Glia', `34` = 'unannotated', `35` = 'unannotated', `36` = 'unannotated',
                             `37` = 'Head-specific cluster', `38` = 'unannotated', `39` = 'unannotated', `4` = 'Head-specific cluster', `40` = 'unannotated', `41` = 'Neuron', `42` = 'Neuron', `43` = 'Head-specific cluster', `44` = 'Neuron', `45` = 'Neuron', `46` = 'Neuron', `47` = 'unannotated',
                             `48` = 'Head-specific cluster', `49` = 'Head-specific cluster', `5`= 'unannotated', `50` = 'Glia', `51` = 'unannotated', `52` = 'unannotated', `53` = 'Head-specific cluster', `54` = 'unannotated', `55` = 'Glia', `56` = 'unannotated', `57` = 'unannotated',
                             `58` = 'unannotated', `59` = 'Neuron', `6` = 'unannotated', `60` = 'Neuron', `61` = 'Neuron', `62` = 'Head-specific cluster', `63` = 'unannotated', `64` = 'Neuron', `65` = 'Neuron', `66` = 'Glia', `67` = 'Neuron', 
                             `68` = 'unannotated', `69` = 'Neuron', `7` = 'unannotated', `70` = 'Neuron', `71` = 'Glia', `72` = 'Head-specific cluster', `73` = 'Glia', `74` = 'unannotated', `75` = 'unannotated', `76` = 'Head-specific cluster', 
                             `77` = 'Neuron', `78` = 'Glia', `79` = 'unannotated', `8` = 'unannotated', `80` = 'Head-specific cluster', `81` = 'unannotated', `82` = 'unannotated', `83` = 'Neuron', `84` = 'unannotated', `85` = 'Glia', `86` = 'Neuron', `87` = 'Glia',
                             `88` = 'Neuron', `89` = 'Glia', `9` = 'unannotated', `90` = 'Head-specific cluster', `91` = 'Glia', `92` = 'Glia', `93` = 'Glia', `94` = 'unannotated', `95` = 'Head-specific cluster', `96` = 'Head-specific cluster', 
                             `97` = 'Neuron', `98` = 'Neuron', `99` = 'unannotated')

my_cols<- c('unannotated' = 'azure4', 'Neuron' = 'deepskyblue', 'Glia' = 'lightpink3', 'Head-specific cluster' = 'burlywood3')

sample_cols<-c('G11A' = '#5F75B8', 'G146C'='#4BBDAE', 'WT' = '#EFAB57')
table(fly.combined$group)
DefaultAssay(fly.combined)<- 'RNA'
                           
fly.combined$seurat_clusters<- Idents(fly.combined)
dev.new()
DimPlot(object = fly.combined, reduction = "umap", group.by = 'seurat_clusters', raster = FALSE)
dev.new()
DimPlot(object = fly.combined, reduction = "umap", group.by = 'seurat_clusters', raster = FALSE, cols = my_cols)

my_cols <- c('1' = 'red', '2' = 'blue')
dev.new()
DimPlot(object = fly.combined, reduction = "umap", group.by = 'batch', raster = FALSE, cols = my_cols)

dev.new()
DimPlot(object = fly.combined, reduction = "umap", group.by = 'group', raster = FALSE, cols = sample_cols)

#Featureplot based on neuronal activity
DefaultAssay(fly.combined)<- 'SCT'
#Glutamatergic
dev.new()
FeaturePlot(fly.combined, features = 'VGlut', raster = FALSE)

#Cholinergic 
dev.new()
FeaturePlot(fly.combined, features = 'VAChT', raster = FALSE)

#Gabaergic
dev.new()
FeaturePlot(fly.combined, features = 'Gad1', raster = FALSE)


DefaultAssay(fly.combined)<- 'SCT'
dev.new()
DotPlot(fly.combined, features = c('VGlut', 'VAChT', 'Gad1'), group.by = 'group')
perc_exp <- DotPlot(fly.combined, features = c('VGlut', 'VAChT', 'Gad1'), group.by= 'group')$data[, c("features.plot", "id", "pct.exp")]
avg_exp <- DotPlot(fly.combined, features = c('VGlut', 'VAChT', 'Gad1'), group.by= 'group')$data[, c("features.plot", "id", "avg.exp")]
dp<-DotPlot(fly.combined, features = c('VGlut', 'VAChT', 'Gad1'), group.by = 'group')

avg_exp_scaled <- DotPlot(fly.combined, features = c('VGlut', 'VAChT', 'Gad1'), group.by= 'group')$data[, c("features.plot", "id", "avg.exp.scaled")]

#Cell proportions and basic cell type stats#
#install.packages('devtools')
#devtools::install_github("rpolicastro/scProportionTest")


###10X Genomics LoupeR for DGEA
#remotes::install_github("10xGenomics/loupeR")
loupeR::setup()

setwd('/Users/laurynhigginson/Desktop/')
library(loupeR)
create_loupe_from_seurat(fly.combined)

#Cell Proportion Analysis
library(scProportionTest)
library(ggplot2)
sc_flyhead<- sc_utils(fly.combined)
prop_test<- permutation_test(
  sc_flyhead, cluster_identity = 'seurat_clusters',
  sample_1 = 'WT', sample_2 = 'G11A',
  sample_identity = 'group'
)

dev.new()  
permutation_plot(prop_test,
                 log2FD_threshold = log2(1.25))

install.packages("svglite")
svglite::svglite() 
table(fly.combined$orig.ident)
library('muscat')
#BiocManager::install("muscat")
head(fly.combined)

DefaultAssay(fly.combined)<- 'RNA'
fly.combined[['SCT']] <- NULL
fly.combined<- as.SingleCellExperiment(fly.combined)



fly.combined <- prepSCE(fly.combined,
                        kid = 'seurat_clusters',
                        gid = 'group',
                        sid = 'orig.ident',
                        drop = TRUE)

nk  <- length(kids <- levels(fly.combined$cluster_id))
ns  <- length(sids <- levels(fly.combined$sample_id))
names(kids) <- kids; names(sids) <- sids

t(table(fly.combined$cluster_id, fly.combined$sample_id))
table(fly.combined$group_id)


# Aggregation of single-cell to pseudobulk data
pb <- aggregateData(fly.combined, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

# one sheet per subpopulation
library(SummarizedExperiment)
assayNames(pb)

# Pseudobulk-level MDS plot
pb_mds <- pbMDS(pb)+
  scale_color_manual(values = my_cols)
?pbMDS

?pbMDS
dev.new()
pb_mds

my_cols
dev.off()

# nb. of cells per cluster-sample
cell.per.samp.tab <- t(table(fly.combined$cluster_id, fly.combined$sample_id))
head(cell.per.samp.tab)

# cell types with at least 25 cells from every each sex/cohort sample
celltype.qc <- colnames(cell.per.samp.tab)[apply(cell.per.samp.tab > 25, 2, sum) == 6]
celltype.qc

### extract pseudobulk information for samples that pass the cell number cutoff
counts.pb <- pb@assays@data[celltype.qc]

# some genes have "gene-" in their name: remove for processing
for (i in 1:length(counts.pb)) {
  rownames(counts.pb[[i]]) <- gsub("gene-", "",rownames(counts.pb[[i]]))
}



counts_h <- counts.pb$`Reticular neuropil associated glia`
counts_h<-round(counts_h)
setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/')
META <- read.table(file ="Higginson_metadata.txt", sep="\t", header=TRUE, row.names=1) 
META
#making sure the META data aligns with counts file
head(META)
rownames(META)<-colnames(counts_h)
rownames(META)
all(rownames(META) %in% colnames(counts_h))

colnames(counts_h)

#making DESeq object
library(DESeq2)
dds<- DESeqDataSetFromMatrix(countData = counts_h,
                             colData = META,
                             design = ~Comparison)

dds<-DESeq(dds)
#number of genes in the matrix/dds object
nrow(dds)
#10986
#11070 Neurons

#Filter to remove rows with no or only 1 count in a row
dds1<- dds[rowSums(counts(dds))>1,]
nrow(dds1)
#9684
#10843 Neurons

#Filter to remove rows with no or at least 25 count in a row
dds2<- dds[rowSums(counts(dds))>25,]
nrow(dds2)
#6898

rld<- rlogTransformation(dds1, blind = TRUE)

#Perform and plot PCA based on data from top X most variable genes (default 500)
dev.new()
plotPCA( rld, intgroup = c('Comparison')) + theme_minimal()

plotPCA (rld, intgroup = c('Comparison'), returnData = TRUE)

?plotPCA

PCAcoordinates<-plotPCA(rld, intgroup = c('Comparison'), ntop = 500, returnData = TRUE)
write.table(PCAcoordinates, file = 'Dmel_PCA_coordinates_full_RN5', sep = '\t')

rlogMat<- assay(rld)
rv = apply(rlogMat, 1, var)
select = order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca = prcomp(t(rlogMat[select,]))
write.table( file = 'Dmel_variances_per_component_full_RN5.txt')


#### Heatmaps ###
sampleDists<- dist(t(assay(rld)))
sampleDists

sampleDistMatrix<- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
library(RColorBrewer)
library(pheatmap)
colors<- colorRampPalette(rev(brewer.pal(9, 'Blues')) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col= colors)


#Heatmap of relative rlog-transformed values across samples. Make Top 10 gene heatmaps
#install.packages('pheatmap')
library(pheatmap)
setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/Article Writing/Pub 1/Figures/Figure 5/')
topVarGenes<- head(order(rowVars(assay(rld)), decreasing = TRUE), 15)
mat<- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df<-as.data.frame(colData(rld)[,c('Comparison')])
df
dev.new()
pheatmap(mat)


#Estimate Dispersion
dds1<-estimateSizeFactors(dds1)


#write the normalized count table to a file. 
counts_norm<-as.data.frame(counts(dds1, normalized = TRUE))
write.table(counts_norm, file - 'Dmel_normed_counts_KC', row.names = TRUE, col.names = TRUE, quote = FALSE, sep = '\t')

#write the raw count table to a file
counts_filtered_raw <- counts(dds1, normalized = FALSE)
write.table(counts_filtered_raw, file = 'Dmel_raw_filtered_kc_counts.txt', row.names = TRUE, col.names = TRUE, quote = FALSE, sep= '\t')

boxplot(counts_norm)
boxplot(counts_norm, outline = F)

##core code##
str(dds1)
#relevel dds1
dds1$Comparison<-relevel(dds1$Comparison, "WT")
dds1<-DESeq(dds1)
#write out DEGs
setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/res5/Cluster_Analysis/Reticular neuropil/') 
res1<- results(dds1, contrast = c("Comparison", "G146C", "WT"))
file_out<-'G146C_vs_WT_DEG_results_rnag.txt'
write.table(res1, file = file_out, sep= '\t')
res2<- results(dds1, contrast = c("Comparison", "G11A", "WT"))
file_out<-'G11A_vs_WT_DEG_results_rnag.txt'
write.table(res2, file = file_out, sep= '\t')

#Making volcano plot
#BiocManager::install('clusterProfiler')
library(AnnotationDbi)
library(clusterProfiler)
library(EnhancedVolcano)

setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/res5/Cluster_Analysis/Hemocyte/Pseudobulk/') 

data<-read.table("G146C_vs_WT_DEG_results_hemocyte.txt")

head(data)
colnames(data)
dev.new()
EnhancedVolcano(data, x = 'log2FoldChange', y = 'pvalue', lab = rownames(data))
dev.new()
select<- c('CG16968', 'l(3)mbn','BomBc1', 'BomS1', 'IM4','fit', 'Arc1', 'IM14', 'CG16978')
nina<- c('Hsp23', 'Hsp22')
dev.new()
EnhancedVolcano(data, x = 'log2FoldChange', y = 'pvalue', lab = rownames(data),
                selectLab = select, drawConnectors = TRUE, arrowheads = FALSE)
?EnhancedVolcano
#GO Analysis
library(org.Dm.eg.db)
library(AnnotationDbi)
FC <- 1-3

pvalue <- 0.05
pvalue
genes_to_test<- rownames(data[data$log2FoldChange<FC,])

#determine how many genes are log2FoldChange>0.5

length(genes_to_test) 

GO_results_G11A_up<-enrichGO(gene= genes_to_test, OrgDb = 'org.Dm.eg.db', keyType = "SYMBOL", ont = 'BP') 
head(GO_results)
head(GO_results)
dev.new()


plot(barplot(GO_results_G11A_up, showCategory = 5))


file_out<-'G11A_vs_WT_GO_results_head.txt'
write.table(GO_results@result, file = file_out, sep= '\t')


library(pheatmap)
setwd('/Users/laurynhigginson/Desktop/')
heatmap<-read.table(file ="heatmap.txt", sep="\t", header=TRUE, row.names=1) 

heatmap(heatmap)
dev.new()
heatmap(as.matrix(heatmap), col = c('brown2', 'orange', 'black'))

#Making an UpSetR plot for DEGs
#install.packages('UpSetR')
library(UpSetR)

#head-specific samples
setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/res5/Neuron_Glia_Head_Other/Cluster_Analysis/Head-specific cluster/Pseudobulk/')
G11A_head<-read.table("G11A_vs_WT_DEG_results_head.xls")
G11A_head_genes<- rownames(G11A_head[G11A_head$log2FoldChange>2,])

G11A_head<-read.table("G11A_vs_WT_DEG_results_head.xls")
G11A_head_genes<- rownames(G11A_head[G11A_head$log2FoldChange>2,])

G146C_head<-read.table("G146C_vs_WT_DEG_results_head.xls")
G146C_head_genes<- rownames(G146C_head[G146C_head$log2FoldChange>2,])

#Neuro-specific samples
setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/res5/Neuron_Glia_Head_Other/Cluster_Analysis/Neurons/Psuedobulk/')

G11A_neuron<-read.table("G11A_vs_WT_DEG_results_neuron.xls")
G11A_neuron_genes<- rownames(G11A_neuron[G11A_neuron$log2FoldChange>2,])

G146C_neuron<-read.table("G146C_vs_WT_DEG_results_neuron.xls")
G146C_neuron_genes<- rownames(G146C_neuron[G146C_neuron$log2FoldChange>2,])

#Glial-specific samples
setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/res5/Neuron_Glia_Head_Other/Cluster_Analysis/Glia/Pseudobulk/')

G11A_glia<-read.table("G11A_vs_WT_DEG_results_glia.xls")
G11A_glia_genes<- rownames(G11A_glia[G11A_glia$log2FoldChange>2,])

G146C_glia<-read.table("G146C_vs_WT_DEG_results_glia.xls")
G146C_glia_genes<- rownames(G146C_glia[G146C_glia$log2FoldChange>2,])


x<-list(
  G11A_neuron = G11A_neuron_genes, 
  G146C_neuron = G146C_neuron_genes,
  G11A_glia = G11A_glia_genes, 
  G146C_glia = G146C_glia_genes
)

head(x)
dev.new()
upset(fromList(x), order.by = 'freq', sets.bar.color = 'black', main.bar.color = 'black', matrix.color = 'black')
?upset

#######
#Sub-clustering cannot be done with a pseudobulk object, if this is the format,reload seurat object
table(fly.combined$seurat_clusters, fly.combined$group)
Idents(fly.combined)<-fly.combined$seurat_clusters
Idents(fly.combined)
kc.seu<-subset(fly.combined, idents = 'Kenyon cells')

#Featureplot based on neuronal activity
DefaultAssay(kc)<- 'SCT'
#Glutamatergic
dev.new()
FeaturePlot(kc, features = 'VGlut', raster = FALSE)

#Cholinergic 
dev.new()
FeaturePlot(kc, features = 'VAChT', raster = FALSE)

#Gabaergic
dev.new()
FeaturePlot(kc, features = 'Gad1', raster = FALSE)
DefaultAssay(kc)<-'SCT'
dev.new()
DotPlot(kc, features = c('VGlut', 'VAChT', 'Gad1'), group.by = 'group')
perc_exp <- DotPlot(kc, features = c('VGlut', 'VAChT', 'Gad1'), group.by= 'group')$data[, c("features.plot", "id", "pct.exp")]
avg_exp <- DotPlot(kc, features = c('VGlut', 'VAChT', 'Gad1'), group.by= 'group')$data[, c("features.plot", "id", "avg.exp")]

#Arc1
dev.new()
FeaturePlot(kc, features = 'VAChT', raster = FALSE)


#redo normalization, PCA, find neighbors/clusters

kc <-SCTransform(kc.seu, method = "glmGamPoi",vars.to.regress = "percent.mt", verbose = FALSE)
kc <- RunPCA(kc, features = VariableFeatures(object = kc))
kc <- FindNeighbors(kc, graph.name = "kc", dims = 1:50)
kc <- FindClusters(kc, resolution = 1, graph.name = "kc")
kc <- RunUMAP(kc, dims = 1:50)

sample_cols<-c('G11A' = '#5F75B8', 'G146C'='#4BBDAE', 'WT' = '#EFAB57')
dev.new()
DimPlot(object = kc, reduction = "umap", group.by = 'group', cols = sample_cols,raster = FALSE)
dev.new()
DimPlot(object = kc, reduction = "umap", group.by = 'batch', raster = FALSE)

#Find cluster markers
kc.markers.res1 <- FindAllMarkers(kc, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.3)

kc_Markers_p_val_adj_0.01<-kc.markers.res1[kc.markers.res1$p_val_adj<=0.01,]

clusters<-c(as.character(unique(kc_Markers_p_val_adj_0.01$cluster)))
mydata_final_top10_markers<-c()
namesofrows<-c()
for (cluster in clusters) {
  namesofrows<-c(namesofrows, paste("Cluster", cluster, sep = "_"))
  mydata_final_top10_markers<-rbind(mydata_final_top10_markers, 
                                    head(kc_Markers_p_val_adj_0.01[kc_Markers_p_val_adj_0.01$cluster==cluster,]$gene, n=10))
}

row.names(mydata_final_top10_markers)<-namesofrows
colnames(mydata_final_top10_markers)<-c("Gene_1", "Gene_2", "Gene_3", "Gene_4", "Gene_5", 
                                        "Gene_6", "Gene_7", "Gene_8", "Gene_9", "Gene_10")

setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/res5/Cluster_Analysis/Kenyon_cells/subclustering/')
write.table(mydata_final_top10_markers, file = "2024May29_kc_sub_res1_pca50_top10marker_genes.txt", quote = F, sep = "\t", 
            row.names = T, col.names = T, append = F)


## RES5 initial clustering on all cell types: RES1 on subclustered Kenyon Cell  ##
kc <- RenameIdents(kc, `0` = "unannotated", `1` = "alpha/beta KC", `2` = "gamma KC",
                   `3` = "unannotated",  `4` = "unannotated", `5` = "alpha'/beta' KC", `6` = "unannotated")

## RES1 initial clustering on all cell types (neuronal from neuronal, glial, hsc, other): ##
kc <- RenameIdents(kc, `0` = "Outer photoreceptors", `1` = "Outer photoreceptors", `2` = "T neuron T4/T5",
                   `3` = "Outer photoreceptors",  `4` = "Columnar neuron T1", `5` = "unannotated", `6` = "Kenyon cells", `7` = 'Photoreceptors R7/R8',
                   `8` = 'unannotated', `9` = 'Lamina monopolar neuron L1', `10` = 'Lamina monopolar neuron L1', `11` = 'Lamina monopolar neuron L2', `12` = 'Centrifugal neuron C3',
                   `12` = 'unannotated', `13` = 'Lamina monopolar neuron L2', `14` = 'Kenyon cells', `15` = 'Lamina monopolar neuron L4', `16` = 'Photoreceptors R7/R8', `17` = 'Medullary intrinsic neuron Mi1',
                   `18` = 'Centrifugal neuron C2', `19` = 'Transmedullary neuron Tm1', `20` = 'T neuron T3', `21` = 'Lamina monopolar neuron L5', `22` = 'Outer photoreceptors', `23` = 'unannotated',
                   `24` = 'Ocellus retinula cell', `25` = 'T neuron T2', `26` = 'Johnston organ neuron', `27` = 'Lamina intrinsic amacrine neuron Lai', `28` = 'unannotated', `29` = 'Distal medullary amacrine neuron Dm9') 

my_cols<- c('Outer photoreceptors' = 'magenta4', 'T neuron T4/T5' = 'bisque1', 'Columnar neuron T1' = 'blue4', 'unannotated' = 'azure4', 'Kenyon cells' = 'red2', 'Photoreceptors R7/R8' = 'aquamarine2', 'Lamina monopolar neuron L1' = 'cadetblue',
            'Lamina monopolar neuron L1' = 'cadetblue1', 'Centrifugal neuron C3' = 'chartreuse1', 'Lamina monopolar neuron L4' = 'coral1', 'Medullary intrinsic neuron Mi1' = 'cyan3', 'Centrifugal neuron C2' = 'darkorange1', 'Transmedullary neuron Tm1' = 'darkorange3',
            'T neuron T3' = 'darkslateblue', 'Lamina monopolar neuron L5' = 'deeppink', 'Ocellus retinula cell' = 'gold2', 'T neuron T2' = 'indianred', 'Johnston organ neuron' = 'lightblue2', 
            'Lamina intrinsic amacrine neuron Lai' = 'lightpink', 'Distal medullary amacrine neuron Dm9' = 'lightsalmon')

#RES1 initial clustering on all cell types (neuronal from neuronal, glial, hsc, other): No subtypes 
kc <- RenameIdents(kc, `0` = "Photoreceptors", `1` = "Photoreceptors", `2` = "T neuron",
                   `3` = "Photoreceptors",  `4` = "Columnar neuron", `5` = "unannotated", `6` = "Kenyon cells", `7` = 'Photoreceptors',
                   `8` = 'unannotated', `9` = 'Lamina monopolar neuron', `10` = 'Lamina monopolar neuron', `11` = 'Lamina monopolar neuron', `12` = 'Centrifugal neuron',
                   `12` = 'unannotated', `13` = 'Lamina monopolar neuron', `14` = 'Kenyon cells', `15` = 'Lamina monopolar neuron', `16` = 'Photoreceptors', `17` = 'Medullary intrinsic neuron Mi1',
                   `18` = 'Centrifugal neuron', `19` = 'Transmedullary neuron Tm1', `20` = 'T neuron', `21` = 'Lamina monopolar neuron', `22` = 'Photoreceptors', `23` = 'unannotated',
                   `24` = 'Ocellus retinula cell', `25` = 'T neuron', `26` = 'Johnston organ neuron', `27` = 'Lamina intrinsic amacrine neuron Lai', `28` = 'unannotated', `29` = 'Distal medullary amacrine neuron Dm9') 

my_cols<- c('Photoreceptors' = 'mediumorchid', 'T neuron' = 'bisque1', 'Columnar neuron' = 'blue4', 'unannotated' = 'azure4', 'Kenyon cells' = 'red2', 'Lamina monopolar neuron' = 'cadetblue',
            'Centrifugal neuron' = 'chartreuse1', 'Medullary intrinsic neuron Mi1' = 'cyan3', 'Transmedullary neuron Tm1' = 'darkorange3',
            'Ocellus retinula cell' = 'gold2', 'Johnston organ neuron' = 'lightblue2', 
            'Lamina intrinsic amacrine neuron Lai' = 'lightpink', 'Distal medullary amacrine neuron Dm9' = 'lightsalmon')


kc$seurat_clusters<- Idents(kc)
dev.new()
DimPlot(object = kc,reduction = "umap", group.by = "seurat_clusters")
DimPlot(object = kc,reduction = "umap", group.by = "seurat_clusters", cols = my_cols, label.size = 0.5)

#####
#Cell proportions and basic cell type stats#
#install.packages('devtools')
#devtools::install_github("rpolicastro/scProportionTest")
library(scProportionTest)
library(ggplot2)
#function to make comparisons between groups using the scProportionTest
sc_kc<- sc_utils(kc)

#online code
prop_test <- permutation_test(
  sc_kc, cluster_identity = "seurat_clusters",
  sample_1 = "WT", sample_2 = "G146C",
  sample_identity = "group"
)

dev.new()
?permutation_plot
permutation_plot(prop_test, log2FD_threshold = .32)
DefaultAssay(kc) <- 'RNA'
dev.new()
DotPlot(fly.combined, features = neuro_dev, group.by = 'group')
DotPlot(fly.combined, features = , group.by = 'group')


## Getting percentage of cells are positive for a transcript using DotPlot
DefaultAssay(kc) <- 'SCT'
perc_ten_g11a<- c('Arc1', 'dyw', 'IM14', 'Cpr65Au', 'Achl', 'Proc')
perc_ten_g146c<- c('milt', 'BomBc1', 'fit', 'pro', 'Listericin', 'Obp99a')
perc_exp <- DotPlot(kc, features=perc_ten_g146c, group.by= 'group')$data[, c("features.plot", "id", "pct.exp")]
head(perc_exp)
library("ggplot2")

dev.new()
ggplot(perc_exp, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() +
  facet_wrap(~features.plot)

##FeaturePlot
dev.new()
FeaturePlot(fly.combined, features = '', raster = FALSE)

setwd('/Users/laurynhigginson/Documents/USC Ph.D/Morton Lab-PhD/WT_G11A_G146C/res5/Cluster_Analysis/Kenyon_cells/PseudoBulk/')
data<-read.table("G146C_vs_WT_DEG_results_kc.xls")

##PPI with STRING##

#BiocManager::install("STRINGdb")
#install.packages("rbioapi")
library(STRINGdb)
library(rbioapi)
#fold change threshold, particularly use for negative numbers

FC<-2-5
genes_to_test<- c(rownames(data[data$log2FoldChange>3,]), rownames(data[data$log2FoldChange<FC,]))
genes_to_test
#Drosophila melanogaster taxonomy ID: 7227

proteins_mapped <- rba_string_map_ids(ids = genes_to_test,
                                      species = 7227)

int_net <- rba_string_interactions_network(ids = genes_to_test,
                                           species = 7227,
                                           required_score = 500)
setwd('/Users/laurynhigginson/Desktop/')
graph_1<-rba_string_network_image(ids = genes_to_test,
                                  image_format = "image",
                                  species = 7227,
                                  save_image = TRUE,
                                  required_score = 500,
                                  network_flavor = "confidence")

#Transcript Percentages
FeaturePlot(fly.combined, features = 'Arc1', raster = FALSE, cols = c('grey', 'cyan1'))
DefaultAssay(fly.combined) <- 'SCT'
perc_exp <- DotPlot(fly.combined, features=c('Arc1'), group.by= 'group')$data[, c("features.plot", "id", "pct.exp")]
head(perc_exp)
library("ggplot2")

dev.new()
ggplot(perc_exp, aes(x=id, y=pct.exp, fill=id)) +
  geom_col() +
  facet_wrap(~features.plot)
dev.new()
FeaturePlot(kc.seu, features = 'Arc1')

