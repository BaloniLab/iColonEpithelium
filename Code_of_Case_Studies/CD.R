install.packages("BiocManager", repos = "https://cloud.r-project.org")
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-data")
# library(SeuratData)
library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratWrappers)


###########
#GSE164985#
###########
# setwd("C:/Users/Alvis/Box/Boyu-Gut_microbiome/CyberGut/refinement/GSE164985")
data_CD189 = Read10X('/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/Individual_samples_ALRA/CD_sample_189/')
data_CD299 = Read10X('/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/Individual_samples_ALRA/CD_sample_299/')
data_CD364 = Read10X('/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/Individual_samples_ALRA/CD_sample_364/')
data_HT206 = Read10X('/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/Individual_samples_ALRA/Health_sample_206/')
data_HT214 = Read10X('/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/Individual_samples_ALRA/Health_sample_214/')
data_HT216 = Read10X('/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/Individual_samples_ALRA/Health_sample_216/')
data_HT217 = Read10X('/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/Individual_samples_ALRA/Health_sample_217/')


Object_CD189 = CreateSeuratObject(counts = data_CD189)
Object_CD299 = CreateSeuratObject(counts = data_CD299)
Object_CD364 = CreateSeuratObject(counts = data_CD364)
Object_HT206 = CreateSeuratObject(counts = data_HT206)
Object_HT214 = CreateSeuratObject(counts = data_HT214)
Object_HT216 = CreateSeuratObject(counts = data_HT216)
Object_HT217 = CreateSeuratObject(counts = data_HT217)


# Sample CD_189
# check data and filter data
# calculate the percentage of mitochondria genes, saved in $percent.mt
Object_CD189[["percent.mt"]] <- PercentageFeatureSet(Object_CD189, pattern = "^MT-") 
# visualization 
# VlnPlot(Object_CD189, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(Object_CD189, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(Object_CD189, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
Object_CD189 <- subset(Object_CD189, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
Object_CD189 <- NormalizeData(Object_CD189, normalization.method = "LogNormalize", scale.factor = 10000)
# Run ALRA
Object_CD189 <- RunALRA(Object_CD189)
# Identification of highly variable features
Object_CD189 <- FindVariableFeatures(Object_CD189, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(Object_CD189)
Object_CD189 <- ScaleData(Object_CD189, features = all.genes)
# Perform linear dimensional reduction
Object_CD189 <- RunPCA(Object_CD189, features = VariableFeatures(object = Object_CD189))
# DimPlot(Object_CD189, reduction = "pca")
# DimHeatmap(Object_CD189, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(Object_CD189, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
Object_CD189 <- JackStraw(Object_CD189, num.replicate = 100)
Object_CD189 <- ScoreJackStraw(Object_CD189, dims = 1:20)
JackStrawPlot(Object_CD189, dims = 1:15)
ElbowPlot(Object_CD189)
# Cluster the cells
Object_CD189 <- FindNeighbors(Object_CD189, dims = 1:13)
Object_CD189 <- FindClusters(Object_CD189, resolution = 0.7)
# Run non-linear dimensional reduction
Object_CD189 <- RunUMAP(Object_CD189, dims = 1:13)
# Visualization
DimPlot(Object_CD189, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
Object_CD189.markers <- FindAllMarkers(Object_CD189, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_CD189.markers, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD189_markers_ALRA.csv')
# Annotate cells
new.cluster.ids <- c('Non-colonocyte','Colonocyte','Colonocyte','Colonocyte','Colonocyte',
                     'Colonocyte','Colonocyte','Non-colonocyte','Colonocyte','Colonocyte',
                     'Non-colonocyte','Colonocyte','Colonocyte','Non-colonocyte')
new.cluster.ids <- c('Colonocyte','Non-colonocyte','Colonocyte','Non-colonocyte','Non-colonocyte',
                     'Non-colonocyte','Colonocyte','Non-colonocyte','Colonocyte','Non-colonocyte',
                     'Non-colonocyte','Colonocyte','Non-colonocyte','Non-colonocyte','Non-colonocyte')

names(new.cluster.ids) <- levels(Object_CD189)
Object_CD189 <- RenameIdents(Object_CD189, new.cluster.ids)
Colonocytes <- WhichCells(Object_CD189, idents = "Colonocyte")
DimPlot(Object_CD189, label=T, label.color = "darkred", label.size = 5,  cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey", reduction = "umap")

# Subset of colonocyte
Colonocyte_CD189 <-  subset(Object_CD189, ident = 'Colonocyte')


# calculate colonocytes markers
Object_CD189_colobocyte_marker = FindMarkers(Object_CD189, ident.1 = 'Colonocyte', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_CD189_colobocyte_marker, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD189_CE_log2fc_ALRA.csv')

# Export Colonocyte expression data
exprs <- GetAssay(Colonocyte_CD189)
exprs_counts <- as.matrix(exprs@data)
write.csv(exprs_counts, '~/cybergut_project/Refinement/GSE164985/Individual_samples/CD189_expression.csv')
Row_median = as.data.frame(rowMedians(as.matrix(exprs_counts)))
colnames(Row_median) = 'Medium_Value'
Row_median['hgnc_symbol'] <- rownames(exprs_counts)
Q30_cutoff <-  quantile(Row_median$Medium_Value, probs = 0.3)
Row_median_over <- filter(Row_median, Medium_Value > 0)


symbolID = Row_median_over$hgnc_symbol
mart <- useMart("ensembl","hsapiens_gene_ensembl")
list <- listAttributes(mart)
Symbol2entrez <- getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters = "hgnc_symbol",values = symbolID, mart = mart)
Row_median_over <- merge(Row_median_over, Symbol2entrez, by = 'hgnc_symbol', all.x = TRUE)
expressiondata <- Row_median_over[,c(3,2,1)]
suffix  = '_AT1'
BIGGID <- paste0(expressiondata$entrezgene_id, suffix)
expressiondata['BIGG_ID'] = BIGGID
write.csv(expressiondata, '~/cybergut_project/Refinement/GSE164985/Individual_samples/CD189_expressionSet.csv')
#_______________________________________________________________________________________________________________________________


# Sample_CD_299
# check data and filter data
# calculate the percentage of mitochondria genes, saved in $percent.mt
Object_CD299[["percent.mt"]] <- PercentageFeatureSet(Object_CD299, pattern = "^MT-") 
# visualization 
# VlnPlot(Object_CD299, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(Object_CD299, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(Object_CD299, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
Object_CD299 <- subset(Object_CD299, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
Object_CD299 <- NormalizeData(Object_CD299, normalization.method = "LogNormalize", scale.factor = 10000)
# Run ALRA
Object_CD299 <- RunALRA(Object_CD299)
# Identification of highly variable features
Object_CD299 <- FindVariableFeatures(Object_CD299, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(Object_CD299)
Object_CD299 <- ScaleData(Object_CD299, features = all.genes)
# Perform linear dimensional reduction
Object_CD299 <- RunPCA(Object_CD299, features = VariableFeatures(object = Object_CD299))
# DimPlot(Object_CD299, reduction = "pca")
# DimHeatmap(Object_CD299, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(Object_CD299, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
Object_CD299 <- JackStraw(Object_CD299, num.replicate = 100)
Object_CD299 <- ScoreJackStraw(Object_CD299, dims = 1:20)
JackStrawPlot(Object_CD299, dims = 1:15)
ElbowPlot(Object_CD299)
# Cluster the cells
Object_CD299 <- FindNeighbors(Object_CD299, dims = 1:11)
Object_CD299 <- FindClusters(Object_CD299, resolution = 0.7)
# Run non-linear dimensional reduction
Object_CD299 <- RunUMAP(Object_CD299, dims = 1:11)
# Visualization
DimPlot(Object_CD299, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
Object_CD299.markers <- FindAllMarkers(Object_CD299, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_CD299.markers, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD299_markers_ALRA.csv')
# Annotate cells
new.cluster.ids <-  c('Colonocyte','Non-colonocyte','Colonocyte','Colonocyte','Colonocyte',
                      'Colonocyte','Colonocyte','Colonocyte','Non-colonocyte')
names(new.cluster.ids) <- levels(Object_CD299)
Object_CD299 <- RenameIdents(Object_CD299, new.cluster.ids)
Colonocytes <- WhichCells(Object_CD299, idents = "Colonocyte")
DimPlot(Object_CD299, label=T, label.color = "darkred", label.size = 5, cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey")

# calculate colonocytes markers
Object_CD299_colobocyte_marker = FindMarkers(Object_CD299, ident.1 = 'Colonocyte', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_CD299_colobocyte_marker, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD299_CE_log2fc_ALRA.csv')

# Subset of colonocyte
Colonocyte_CD299 <-  subset(Object_CD299, ident = 'Colonocyte')

# Export Colonocyte expression data
exprs <- GetAssay(Colonocyte_CD299)
exprs_counts <- as.matrix(exprs@data)
write.csv(exprs_counts, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD299_expression_original.csv')
Row_median = as.data.frame(rowMedians(as.matrix(exprs_counts)))
colnames(Row_median) = 'Medium_Value'
Row_median['hgnc_symbol'] <- rownames(exprs_counts)
Q30_cutoff <-  quantile(Row_median$Medium_Value, probs = 0.3)
Row_median_over <- filter(Row_median, Medium_Value > 0)

symbolID = Row_median_over$hgnc_symbol
mart <- useMart("ensembl","hsapiens_gene_ensembl")
list <- listAttributes(mart)
Symbol2entrez <- getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters = "hgnc_symbol",values = symbolID, mart = mart)
Row_median_over <- merge(Row_median_over, Symbol2entrez, by = 'hgnc_symbol', all.x = TRUE)
expressiondata <- Row_median_over[,c(3,2,1)]
write.csv(expressiondata, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD299_expressionSet.csv')
#_______________________________________________________________________________________________________________________________


# Sample CD_364
# check data and filter data
# calculate the percentage of mitochondria genes, saved in $percent.mt
Object_CD364[["percent.mt"]] <- PercentageFeatureSet(Object_CD364, pattern = "^MT-") 
# visualization 
# VlnPlot(Object_CD364, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(Object_CD364, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(Object_CD364, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
Object_CD364 <- subset(Object_CD364, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
Object_CD364 <- NormalizeData(Object_CD364, normalization.method = "LogNormalize", scale.factor = 10000)
# run alra
Object_CD364 <- RunALRA(Object_CD364)
# Identification of highly variable features
Object_CD364 <- FindVariableFeatures(Object_CD364, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(Object_CD364)
Object_CD364 <- ScaleData(Object_CD364, features = all.genes)
# Perform linear dimensional reduction
Object_CD364 <- RunPCA(Object_CD364, features = VariableFeatures(object = Object_CD364))
# DimPlot(Object_CD364, reduction = "pca")
# DimHeatmap(Object_CD364, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(Object_CD364, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
Object_CD364 <- JackStraw(Object_CD364, num.replicate = 100)
Object_CD364 <- ScoreJackStraw(Object_CD364, dims = 1:20)
JackStrawPlot(Object_CD364, dims = 1:15)
ElbowPlot(Object_CD364)
# Cluster the cells
Object_CD364 <- FindNeighbors(Object_CD364, dims = 1:13)
Object_CD364 <- FindClusters(Object_CD364, resolution = 0.7)
# Run non-linear dimensional reduction
Object_CD364 <- RunUMAP(Object_CD364, dims = 1:13)
# Visualization
DimPlot(Object_CD364, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
Object_CD364.markers <- FindAllMarkers(Object_CD364, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_CD364.markers, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD364_markers_ALRA.csv')
# Annotate cells
new.cluster.ids <- c('Non-colonocyte','Colonocyte','Colonocyte','Colonocyte','Colonocyte',
                     'Colonocyte','Colonocyte','Non-colonocyte','Colonocyte','Non-colonocyte',
                     'Colonocyte','Colonocyte','Non-colonocyte','Non-colonocyte')
names(new.cluster.ids) <- levels(Object_CD364)
Object_CD364 <- RenameIdents(Object_CD364, new.cluster.ids)
Colonocytes <- WhichCells(Object_CD364, idents = "Colonocyte")
DimPlot(Object_CD364, label=T, label.color = "darkred", label.size = 5, cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey")


# calculate colonocytes markers
Object_CD364_colobocyte_marker = FindMarkers(Object_CD364, ident.1 = 'Colonocyte', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_CD364_colobocyte_marker, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD364_CE_log2fc_ALRA.csv')

# Subset of colonocyte
Colonocyte_CD364 <-  subset(Object_CD364, ident = 'Colonocyte')

# Export Colonocyte expression data
exprs <- GetAssay(Colonocyte_CD364)
exprs_counts <- as.matrix(exprs@data)
write.csv(exprs_counts, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD364_expression_original.csv')
Row_median = as.data.frame(rowMedians(as.matrix(exprs_counts)))
colnames(Row_median) = 'Medium_Value'
Row_median['hgnc_symbol'] <- rownames(exprs_counts)
Q30_cutoff <-  quantile(Row_median$Medium_Value, probs = 0.3)
Row_median_over <- filter(Row_median, Medium_Value > 0)

symbolID = Row_median_over$hgnc_symbol
mart <- useMart("ensembl","hsapiens_gene_ensembl")
list <- listAttributes(mart)
Symbol2entrez <- getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters = "hgnc_symbol",values = symbolID, mart = mart)
Row_median_over <- merge(Row_median_over, Symbol2entrez, by = 'hgnc_symbol', all.x = TRUE)
expressiondata <- Row_median_over[,c(3,2,1)]
suffix  = '_AT1'
BIGGID <- paste0(expressiondata$entrezgene_id, suffix)
expressiondata['BIGG_ID'] = BIGGID
write.csv(expressiondata, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/CD364_expressionSet.csv')

#_______________________________________________________________________________________________________________________________

# Sample HT_206
# check data and filter data
# calculate the percentage of mitochondria genes, saved in $percent.mt
Object_HT206[["percent.mt"]] <- PercentageFeatureSet(Object_HT206, pattern = "^MT-") 
# visualization 
# VlnPlot(Object_HT206, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(Object_HT206, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(Object_HT206, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
Object_HT206 <- subset(Object_HT206, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
Object_HT206 <- NormalizeData(Object_HT206, normalization.method = "LogNormalize", scale.factor = 10000)
#run alra
Object_HT206 <- RunALRA(Object_HT206)
# Identification of highly variable features
Object_HT206 <- FindVariableFeatures(Object_HT206, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(Object_HT206)
Object_HT206 <- ScaleData(Object_HT206, features = all.genes)
# Perform linear dimensional reduction
Object_HT206 <- RunPCA(Object_HT206, features = VariableFeatures(object = Object_HT206))
# DimPlot(Object_HT206, reduction = "pca")
# DimHeatmap(Object_HT206, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(Object_HT206, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
Object_HT206 <- JackStraw(Object_HT206, num.replicate = 100)
Object_HT206 <- ScoreJackStraw(Object_HT206, dims = 1:20)
JackStrawPlot(Object_HT206, dims = 1:15)
ElbowPlot(Object_HT206)
# Cluster the cells
Object_HT206 <- FindNeighbors(Object_HT206, dims = 1:13)
Object_HT206 <- FindClusters(Object_HT206, resolution = 0.7)
# Run non-linear dimensional reduction
Object_HT206 <- RunUMAP(Object_HT206, dims = 1:13)
# Visualization
DimPlot(Object_HT206, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
Object_HT206.markers <- FindAllMarkers(Object_HT206, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_HT206.markers, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT206_markers_ALRA.csv')
# Annotate cells
new.cluster.ids <- c('Colonocyte','Non-colonocyte','Colonocyte','Colonocyte','Colonocyte',
                     'Colonocyte','Colonocyte','Non-colonocyte','Non-colonocyte','Non-colonocyte',
                     'Colonocyte','Colonocyte','Colonocyte','Non-colonocyte','Non-colonocyte',
                     'Non-colonocyte','Non-colonocyte','Non-colonocyte')
new.cluster.ids <- c('Colonocyte','Non-colonocyte','Colonocyte','Non-colonocyte','Colonocyte',
                     'Colonocyte','Non-colonocyte','Non-colonocyte','Non-colonocyte','Non-colonocyte',
                     'Colonocyte','Non-colonocyte','Colonocyte','Non-colonocyte',
                     'Non-colonocyte','Non-colonocyte','Non-colonocyte','Colonocyte')
names(new.cluster.ids) <- levels(Object_HT206)
Object_HT206 <- RenameIdents(Object_HT206, new.cluster.ids)
Colonocytes <- WhichCells(Object_HT206, idents = "Colonocyte")
DimPlot(Object_HT206, label=T, label.color = "darkred", label.size = 5, cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey")

# calculate colonocytes markers
Object_HT206_colobocyte_marker = FindMarkers(Object_HT206, ident.1 = 'Colonocyte', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_HT206_colobocyte_marker, 'HT206_CE_log2fc_ALRA.csv')

# Subset of colonocyte
Colonocyte_HT206 <-  subset(Object_HT206, ident = 'Colonocyte')

# Export Colonocyte expression data
exprs <- GetAssay(Colonocyte_HT206)
exprs_counts <- as.matrix(exprs@data)
write.csv(exprs_counts, '~/cybergut_project/Refinement/GSE164985/Individual_samples/HT206_expression_original.csv')
Row_median['hgnc_symbol'] <- rownames(exprs_counts)
Q30_cutoff <-  quantile(Row_median$Medium_Value, probs = 0.3)
Row_median_over <- filter(Row_median, Medium_Value > 0)

symbolID = Row_median_over$hgnc_symbol
mart <- useMart("ensembl","hsapiens_gene_ensembl")
list <- listAttributes(mart)
Symbol2entrez <- getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters = "hgnc_symbol",values = symbolID, mart = mart)
Row_median_over <- merge(Row_median_over, Symbol2entrez, by = 'hgnc_symbol', all.x = TRUE)
expressiondata <- Row_median_over[,c(3,2,1)]
suffix  = '_AT1'
BIGGID <- paste0(expressiondata$entrezgene_id, suffix)
expressiondata['BIGG_ID'] = BIGGID
write.csv(expressiondata, '~/cybergut_project/Refinement/GSE164985/Individual_samples/HT206_expressionSet.csv')


#_______________________________________________________________________________________________________________________________


# Sample HT_214
# check data and filter data
# calculate the percentage of mitochondria genes, saved in $percent.mt
Object_HT214[["percent.mt"]] <- PercentageFeatureSet(Object_HT214, pattern = "^MT-") 
# visualization 
# VlnPlot(Object_HT214, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(Object_HT214, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(Object_HT214, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
Object_HT214 <- subset(Object_HT214, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
Object_HT214 <- NormalizeData(Object_HT214, normalization.method = "LogNormalize", scale.factor = 10000)
# run alra
Object_HT214 <- RunALRA(Object_HT214)
# Identification of highly variable features
Object_HT214 <- FindVariableFeatures(Object_HT214, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(Object_HT214)
Object_HT214 <- ScaleData(Object_HT214, features = all.genes)
# Perform linear dimensional reduction
Object_HT214 <- RunPCA(Object_HT214, features = VariableFeatures(object = Object_HT214))
# DimPlot(Object_HT214, reduction = "pca")
# DimHeatmap(Object_HT214, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(Object_HT214, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
Object_HT214 <- JackStraw(Object_HT214, num.replicate = 100)
Object_HT214 <- ScoreJackStraw(Object_HT214, dims = 1:20)
JackStrawPlot(Object_HT214, dims = 1:15)
ElbowPlot(Object_HT214)
# Cluster the cells
Object_HT214 <- FindNeighbors(Object_HT214, dims = 1:13)
Object_HT214 <- FindClusters(Object_HT214, resolution = 0.7)
# Run non-linear dimensional reduction
Object_HT214 <- RunUMAP(Object_HT214, dims = 1:13)
# Visualization
DimPlot(Object_HT214, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
Object_HT214.markers <- FindAllMarkers(Object_HT214, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_HT214.markers, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT214_markers_ALRA.csv')
# Annotate cells
new.cluster.ids <- c('Colonocyte','Colonocyte','Colonocyte','Colonocyte','Colonocyte',
                     'Colonocyte','Colonocyte','Non-colonocyte','Colonocyte','Non-colonocyte',
                     'Colonocyte','Colonocyte','Non-colonocyte','Colonocyte','Colonocyte')
names(new.cluster.ids) <- levels(Object_HT214)
Object_HT214 <- RenameIdents(Object_HT214, new.cluster.ids)
Colonocytes <- WhichCells(Object_HT214, idents = "Colonocyte")
DimPlot(Object_HT214, label=T, label.color = "darkred", label.size = 5, cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey")


# calculate colonocytes markers
Object_HT214_colobocyte_marker = FindMarkers(Object_HT214, ident.1 = 'Colonocyte', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_HT214_colobocyte_marker, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT214_CE_log2fc_ALRA.csv')

# Subset of colonocyte
Colonocyte_HT214 <-  subset(Object_HT214, ident = 'Colonocyte')

# Export Colonocyte expression data
exprs <- GetAssay(Colonocyte_HT214)
exprs_counts <- as.matrix(exprs@data)
write.csv(exprs_counts, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT214_expression_original.csv')
Row_median = as.data.frame(rowMedians(as.matrix(exprs_counts)))
colnames(Row_median) = 'Medium_Value'
Row_median['hgnc_symbol'] <- rownames(exprs_counts)
Q30_cutoff <-  quantile(Row_median$Medium_Value, probs = 0.3)
Row_median_over <- filter(Row_median, Medium_Value > 0)

symbolID = Row_median_over$hgnc_symbol
mart <- useMart("ensembl","hsapiens_gene_ensembl")
list <- listAttributes(mart)
Symbol2entrez <- getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters = "hgnc_symbol",values = symbolID, mart = mart)
Row_median_over <- merge(Row_median_over, Symbol2entrez, by = 'hgnc_symbol', all.x = TRUE)
expressiondata <- Row_median_over[,c(3,2,1)]
suffix  = '_AT1'
BIGGID <- paste0(expressiondata$entrezgene_id, suffix)
expressiondata['BIGG_ID'] = BIGGID
write.csv(expressiondata, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT214_expressionSet.csv')

#_______________________________________________________________________________________________________________________________


# Sample HT_216
# check data and filter data
# calculate the percentage of mitochondria genes, saved in $percent.mt
Object_HT216[["percent.mt"]] <- PercentageFeatureSet(Object_HT216, pattern = "^MT-") 
# visualization 
# VlnPlot(Object_HT216, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(Object_HT216, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(Object_HT216, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
Object_HT216 <- subset(Object_HT216, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
Object_HT216 <- NormalizeData(Object_HT216, normalization.method = "LogNormalize", scale.factor = 10000)
# RUN alra
Object_HT216 <- RunALRA(Object_HT216)
# Identification of highly variable features
Object_HT216 <- FindVariableFeatures(Object_HT216, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(Object_HT216)
Object_HT216 <- ScaleData(Object_HT216, features = all.genes)
# Perform linear dimensional reduction
Object_HT216 <- RunPCA(Object_HT216, features = VariableFeatures(object = Object_HT216))
# DimPlot(Object_HT216, reduction = "pca")
# DimHeatmap(Object_HT216, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(Object_HT216, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
Object_HT216 <- JackStraw(Object_HT216, num.replicate = 100)
Object_HT216 <- ScoreJackStraw(Object_HT216, dims = 1:20)
JackStrawPlot(Object_HT216, dims = 1:15)
ElbowPlot(Object_HT216)
# Cluster the cells
Object_HT216 <- FindNeighbors(Object_HT216, dims = 1:13)
Object_HT216 <- FindClusters(Object_HT216, resolution = 0.7)
# Run non-linear dimensional reduction
Object_HT216 <- RunUMAP(Object_HT216, dims = 1:13)
# Visualization
DimPlot(Object_HT216, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
Object_HT216.markers <- FindAllMarkers(Object_HT216, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_HT216.markers, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT216_markers_ALRA.csv')
# Annotate cells
new.cluster.ids <- c('Non-colonocyte','Colonocyte','Colonocyte','Colonocyte','Colonocyte',
                     'Colonocyte','Colonocyte','Colonocyte','Non-colonocyte','Colonocyte')
names(new.cluster.ids) <- levels(Object_HT216)
Object_HT216 <- RenameIdents(Object_HT216, new.cluster.ids)
Colonocytes <- WhichCells(Object_HT216, idents = "Colonocyte")
DimPlot(Object_HT216, label=T, label.color = "darkred", label.size = 5, cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey")


# calculate colonocytes markers
Object_HT216_colobocyte_marker = FindMarkers(Object_HT216, ident.1 = 'Colonocyte', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_HT216_colobocyte_marker, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT216_CE_log2fc_ALRA.csv')

# Subset of colonocyte
Colonocyte_HT216 <-  subset(Object_HT216, ident = 'Colonocyte')

# Export Colonocyte expression data
exprs <- GetAssay(Colonocyte_HT216)
exprs_counts <- as.matrix(exprs@data)
write.csv(exprs_counts, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT216_expression_original.csv')
Row_median = as.data.frame(rowMedians(as.matrix(exprs_counts)))
colnames(Row_median) = 'Medium_Value'
Row_median['hgnc_symbol'] <- rownames(exprs_counts)
Q30_cutoff <-  quantile(Row_median$Medium_Value, probs = 0.3)
Row_median_over <- filter(Row_median, Medium_Value > 0)

symbolID = Row_median_over$hgnc_symbol
mart <- useMart("ensembl","hsapiens_gene_ensembl")
list <- listAttributes(mart)
Symbol2entrez <- getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters = "hgnc_symbol",values = symbolID, mart = mart)
Row_median_over <- merge(Row_median_over, Symbol2entrez, by = 'hgnc_symbol', all.x = TRUE)
expressiondata <- Row_median_over[,c(3,2,1)]
suffix  = '_AT1'
BIGGID <- paste0(expressiondata$entrezgene_id, suffix)
expressiondata['BIGG_ID'] = BIGGID
write.csv(expressiondata, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT216_expressionSet.csv')

#_______________________________________________________________________________________________________________________________


# Sample HT_217
# check data and filter data
# calculate the percentage of mitochondria genes, saved in $percent.mt
Object_HT217[["percent.mt"]] <- PercentageFeatureSet(Object_HT217, pattern = "^MT-") 
# visualization 
# VlnPlot(Object_HT217, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(Object_HT217, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(Object_HT217, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
Object_HT217 <- subset(Object_HT217, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
Object_HT217 <- NormalizeData(Object_HT217, normalization.method = "LogNormalize", scale.factor = 10000)
# Run alra
Object_HT217 <- RunALRA(Object_HT217)
# Identification of highly variable features
Object_HT217 <- FindVariableFeatures(Object_HT217, selection.method = "vst", nfeatures = 2000)
# Scaling the data
all.genes <- rownames(Object_HT217)
Object_HT217 <- ScaleData(Object_HT217, features = all.genes)
# Perform linear dimensional reduction
Object_HT217 <- RunPCA(Object_HT217, features = VariableFeatures(object = Object_HT217))
# DimPlot(Object_HT217, reduction = "pca")
# DimHeatmap(Object_HT217, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(Object_HT217, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
Object_HT217 <- JackStraw(Object_HT217, num.replicate = 100)
Object_HT217 <- ScoreJackStraw(Object_HT217, dims = 1:20)
JackStrawPlot(Object_HT217, dims = 1:15)
ElbowPlot(Object_HT217)
# Cluster the cells
Object_HT217 <- FindNeighbors(Object_HT217, dims = 1:9)
Object_HT217 <- FindClusters(Object_HT217, resolution = 0.7)
# Run non-linear dimensional reduction
Object_HT217 <- RunUMAP(Object_HT217, dims = 1:9)
# Visualization
DimPlot(Object_HT217, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
Object_HT217.markers <- FindAllMarkers(Object_HT217, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_HT217.markers, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT217_markers_ALRA.csv')
# Annotate cells
new.cluster.ids <- c('Non-colonocyte','Colonocyte','Non-colonocyte','Colonocyte','Colonocyte',
                     'Colonocyte','Colonocyte','Non-colonocyte','Non-colonocyte','Colonocyte',
                     'Colonocyte','Colonocyte','Colonocyte','Colonocyte','Colonocyte')
names(new.cluster.ids) <- levels(Object_HT217)
Object_HT217 <- RenameIdents(Object_HT217, new.cluster.ids)
Colonocytes <- WhichCells(Object_HT217, idents = "Colonocyte")
DimPlot(Object_HT217, label=T, label.color = "darkred", cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey")


# calculate colonocytes markers
Object_HT217_colobocyte_marker = FindMarkers(Object_HT217, ident.1 = 'Colonocyte', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Object_HT217_colobocyte_marker, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT217_CE_log2fc_ALRA.csv')

# Subset of colonocyte
Colonocyte_HT217 <-  subset(Object_HT217, ident = 'Colonocyte')

# Export Colonocyte expression data
exprs <- GetAssay(Colonocyte_HT217)
exprs_counts <- as.matrix(exprs@data)
write.csv(exprs_counts, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT217_expression_original.csv')
Row_median = as.data.frame(rowMedians(as.matrix(exprs_counts)))
colnames(Row_median) = 'Medium_Value'
Row_median['hgnc_symbol'] <- rownames(exprs_counts)
Q30_cutoff <-  quantile(Row_median$Medium_Value, probs = 0.3)
Row_median_over <- filter(Row_median, Medium_Value > 0)

symbolID = Row_median_over$hgnc_symbol
mart <- useMart("ensembl","hsapiens_gene_ensembl")
list <- listAttributes(mart)
Symbol2entrez <- getBM(attributes=c("hgnc_symbol","entrezgene_id"),filters = "hgnc_symbol",values = symbolID, mart = mart)
Row_median_over <- merge(Row_median_over, Symbol2entrez, by = 'hgnc_symbol', all.x = TRUE)
expressiondata <- Row_median_over[,c(3,2,1)]
suffix  = '_AT1'
BIGGID <- paste0(expressiondata$entrezgene_id, suffix)
expressiondata['BIGG_ID'] = BIGGID
write.csv(expressiondata, '~/cybergut_project/Refinement/GSE164985/Individual_samples_ALRA/HT217_expressionSet.csv')

#_______________________________________________________________________________________________________________________________

# CD_merge
# Merging Seurat Objects
CD_merge <- merge(Object_CD189, y = c(Object_CD299, Object_CD364), add.cell.ids = c('CD189', 'CD299', 'CD364'), project = 'GSE164985')
head(colnames(CD_merge))
tail(colnames(CD_merge))
unique(sapply(X = strsplit(colnames(CD_merge), split = "_"), FUN = "[", 1))
table(CD_merge$orig.ident)

CD_merge[["percent.mt"]] <- PercentageFeatureSet(CD_merge, pattern = "^MT-")
# plot1 <- FeatureScatter(CD_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(CD_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
CD_merge <- subset(CD_merge, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
CD_merge <- NormalizeData(CD_merge, normalization.method = "LogNormalize", scale.factor = 10000)
# Run ALRA
CD_merge <- RunALRA(CD_merge)
# Identification of highly variable features
CD_merge <- FindVariableFeatures(CD_merge, selection.method = "vst", nfeatures = 3000)
# Scaling the data
all.genes <- rownames(CD_merge)
CD_merge <- ScaleData(CD_merge, features = all.genes)
# Perform linear dimensional reduction
CD_merge <- RunPCA(CD_merge, features = VariableFeatures(object = CD_merge))
# DimPlot(CD_merge, reduction = "pca")
# DimHeatmap(CD_merge, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(CD_merge, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
CD_merge <- JackStraw(CD_merge, num.replicate = 100)
CD_merge <- ScoreJackStraw(CD_merge, dims = 1:20)
JackStrawPlot(CD_merge, dims = 1:15)
ElbowPlot(CD_merge)
# Cluster cells
CD_merge <- FindNeighbors(CD_merge, dims = 1:15)
CD_merge <- FindClusters(CD_merge, resolution = 0.7)
# Run non-linear dimensional reduction
CD_merge <- RunUMAP(CD_merge, dims = 1:15)
# Visualization
DimPlot(CD_merge, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
CD_merge.markers <- FindAllMarkers(CD_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CD_merge.markers, '~/cybergut_project/Case_Studies/GSE164985/Merged_sample/feature_3000/CD_merge_markers_feature_3k.csv')
# Annotate cells

new.cluster.ids <- c("CD_Colonocyte","CD_Colonocyte","CD_Colonocyte","CD_Colonocyte","CD_Colonocyte",
                     "CD_Colonocyte","Non-Colonocyte","CD_Colonocyte","Non-Colonocyte","Non-Colonocyte",
                     "Non-Colonocyte","CD_Colonocyte","CD_Colonocyte","Non-Colonocyte","Non-Colonocyte",
                     "CD_Colonocyte","Non-Colonocyte","CD_Colonocyte","Non-Colonocyte")

names(new.cluster.ids) <- levels(CD_merge)
CD_merge <- RenameIdents(CD_merge, new.cluster.ids)
CD_Colonocytes <- WhichCells(CD_merge, idents = "CD_Colonocyte")
DimPlot(CD_merge, label=T, label.color = "darkred", cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey")

CD_colonocyte = subset(CD_merge,  idents = 'CD_Colonocyte')

# calculate colonocytes markers
CD_merge_colonocyte_marker = FindMarkers(CD_merge, ident.1 = 'Colonocyte', only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CD_merge_colonocyte_marker, '~/cybergut_project/Case_Studies/GSE164985/Merged_sample/feature_3000/CD_merge_colonocyte_markter_feature_3k.csv')
#_______________________________________________________________________________________________________________________________

# HT_merge
#Merging Seurat Objects
HT_merge<- merge(Object_HT206, y = c(Object_HT214, Object_HT216, Object_HT217), add.cell.ids = c('HT206', 'HT214', 'HT216', 'HT217'), project = 'GSE164985')
head(colnames(HT_merge))
tail(colnames(HT_merge))
unique(sapply(X = strsplit(colnames(HT_merge), split = "_"), FUN = "[", 1))
table(HT_merge$orig.ident)

HT_merge[["percent.mt"]] <- PercentageFeatureSet(HT_merge, pattern = "^MT-")
# plot1 <- FeatureScatter(HT_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(HT_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# filter cell with too few/many genes and high precent age of mitochondria genes
HT_merge<- subset(HT_merge, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
# normalization 
HT_merge<- NormalizeData(HT_merge, normalization.method = "LogNormalize", scale.factor = 10000)
# Run ALRA
HT_merge<- RunALRA(HT_merge)
# Identification of highly variable features
HT_merge<- FindVariableFeatures(HT_merge, selection.method = "vst", nfeatures = 3000)
# Scaling the data
all.genes <- rownames(HT_merge)
HT_merge<- ScaleData(HT_merge, features = all.genes)
# Perform linear dimensional reduction
HT_merge<- RunPCA(HT_merge, features = VariableFeatures(object = HT_merge))
# DimPlot(HT_merge, reduction = "pca")
# DimHeatmap(HT_merge, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(HT_merge, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the dimensionality of the dataset
HT_merge<- JackStraw(HT_merge, num.replicate = 100)
HT_merge<- ScoreJackStraw(HT_merge, dims = 1:20)
JackStrawPlot(HT_merge, dims = 1:15)
ElbowPlot(HT_merge)
# Cluster cells
HT_merge<- FindNeighbors(HT_merge, dims = 1:12)
HT_merge<- FindClusters(HT_merge, resolution = 0.7)
# Run non-linear dimensional reduction
HT_merge<- RunUMAP(HT_merge, dims = 1:12)
# Visualization
DimPlot(HT_merge, reduction = "umap")
# Finding differentially expressed features 
# find markers for every cluster compared to all remaining cells, report only the positive one
HT_merge.markers <- FindAllMarkers(HT_merge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(HT_merge.markers, '~/cybergut_project/Case_Studies/GSE164985/Merged_sample/feature_3000/HT_merge_markers_feature_3k.csv')
# Annotate cells
new.cluster.ids <- c("Non-Colonocyte","HT_Colonocyte","HT_Colonocyte","HT_Colonocyte","Non-Colonocyte","HT_Colonocyte",
                     "Non-Colonocyte","HT_Colonocyte","Non-Colonocyte","HT_Colonocyte","HT_Colonocyte",
                     "HT_Colonocyte","HT_Colonocyte","HT_Colonocyte","Non-Colonocyte","HT_Colonocyte",
                     "Non-Colonocyte","Non-Colonocyte","Non-Colonocyte","Non-Colonocyte")

names(new.cluster.ids) <- levels(HT_merge)
HT_merge<- RenameIdents(HT_merge, new.cluster.ids)
Colonocytes <- WhichCells(HT_merge, idents = "Colonocyte")
Colonocytes <- subset(HT_merge, idents = "Colonocyte")
DimPlot(HT_merge, label=T, label.color = "darkred", cells.highlight= Colonocytes, cols.highlight = c("darkblue"), cols= "grey")
HT_colonocyte = subset(HT_merge,  idents = 'HT_Colonocyte')

# calculate colonocytes markers
HT_merge_colobocyte_marker = FindMarkers(HT_merge, ident.1 = 'Colonocyte', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(HT_merge_colobocyte_marker, '~/cybergut_project/Case_Studies/GSE164985/Merged_sample/feature_3000/HT_merge_colonocyte_marker_feature_3k_POSITIVE.csv')
#_______________________________________________________________________________________________________________________________
# merge CD and HT colonocyte
colonocyte_merge  <- merge(CD_colonocyte, y = HT_colonocyte, add.cell.ids = c("CD", "HT"), project = "GSE164985")
# without normalization
DEG1 = FindMarkers(colonocyte_merge, ident.1 = 'CD_Colonocyte', only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DEG1, '~/cybergut_project/Case_Studies/GSE164985/HT_VS_CD_colonocyte/DEG_HTCD_NO_norm.csv')
#_______________________________________________________________________________________________________________________________
# with normalizatio
HT_merge<- merge(Object_HT206, y = c(Object_HT214, Object_HT216, Object_HT217), add.cell.ids = c('HT206', 'HT214', 'HT216', 'HT217'), project = 'GSE164985')
CD_merge <- merge(Object_CD189, y = c(Object_CD299, Object_CD364), add.cell.ids = c('CD189', 'CD299', 'CD364'), project = 'GSE164985')

seur = merge(HT_merge, y = CD_merge)
seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-") 
seur <- subset(seur, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
seur <- RunALRA(seur)
saveRDS(seur, '/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/HT_VS_CD_colonocyte/CD_seur_alra.rds')
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(seur)
seur <- ScaleData(seur, features = all.genes)
seur <- RunPCA(seur, features = VariableFeatures(object = seur))
seur <- JackStraw(seur, num.replicate = 100)
seur <- ScoreJackStraw(seur, dims = 1:20)
saveRDS(seur, '/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/HT_VS_CD_colonocyte/CD_seur_alra_jackstraw.rds')
JackStrawPlot(seur, dims = 1:15)
ElbowPlot(seur)
seur <- FindNeighbors(seur, dims = 1:15)
seur <- FindClusters(seur, resolution = 0.5)
seur <- RunUMAP(seur, dims = 1:15)
DimPlot(seur, reduction = "umap")

HT_colonocyte <- readRDS('~/cybergut_project/Case_Studies/GSE164985/HT_VS_CD_colonocyte/HT_colonocyte.rds')
CD_colonocyte <- readRDS('~/cybergut_project/Case_Studies/GSE164985/HT_VS_CD_colonocyte/CD_colonocyte.rds')

seur[["CellName"]] <- colnames(seur)
Name_CD_colonocyte <- colnames(CD_colonocyte)
Name_HT_colonocyte <- colnames(HT_colonocyte)

DimPlot(seur, label=T, label.color = 'grey', label.size = 5,  cells.highlight= list(Name_CD_colonocyte, Name_HT_colonocyte), cols.highlight = c("darkblue", "darkred"), cols= "grey")

HT_colonocytes <- subset(seur, subset = CellName %in% Name_HT_colonocyte)
CD_cxolonocytes <- subset(seur, subset = CellName %in% Name_CD_colonocyte)
saveRDS(HT_colonocytes, '/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/HT_VS_CD_colonocyte/HT_colonocytes_norm.rds')
saveRDS(CD_cxolonocytes, '/scratch/negishi/jiang817/CyberGut_project/Case_studies/GSE164985/HT_VS_CD_colonocyte/CD_colonocytes_norm.rds')


Idents(CD_cxolonocytes) <- "CD_colonocyte"
Idents(HT_colonocytes) <- "HT_colonocyte"
colonocyte_merge  <- merge(CD_cxolonocytes, y = HT_colonocytes)
DEG1 = FindMarkers(colonocyte_merge, ident.1 = 'CD_colonocyte', only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(DEG1, '~/cybergut_project/Case_Studies/GSE164985/HT_VS_CD_colonocyte/DEG_HTCD_with_norm.csv')

# t-test 
library(dplyr)
data <- read.csv('CD_FVA_TEST.csv', row.names = 1)
data <- data %>% 
  mutate(P_value=NA)

for(i in 1: 4628){
  print(rownames(data)[i])
  fluxes <- data[i,1:7]
  num_na <- sum(is.na(fluxes))
  if(num_na < 4){
    fluxes[is.na(fluxes)] <- 0
    fluxes_value <- as.vector(t(fluxes))
    duplicata_num <- sum(duplicated(fluxes_value))
    if(duplicata_num >=6){data[i,8] <- 'Na'}
    else{
      fluxes <- as.data.frame(t(fluxes))
      colnames(fluxes) <- c('value')
      fluxes$group <- c('CD', 'CD', 'CD', 'HT', 'HT', 'HT', 'HT')
      fluxes$group <- factor(fluxes$group)
      t_test <- t.test(value~group, fluxes, paired = FALSE, alternative = 'two.sided')
      p_value <- t_test$p.value
      data[i,8] <- p_value
    }
    
  }else{
    data[i,8] <- 'Na'}
}

write.csv(data, 't_test_outputs.csv')

#_______________________________________________________________________________________________________________________________
# DEG
# add meta data

Object_CD189@meta.data$Sample = 'CD189'
Object_CD299@meta.data$Sample = 'CD299'
Object_CD364@meta.data$Sample = 'CD364'
Object_HT206@meta.data$Sample = 'HT206'
Object_HT214@meta.data$Sample = 'HT214'
Object_HT216@meta.data$Sample = 'HT216'
Object_HT217@meta.data$Sample = 'HT217'

CD_merge <- merge(Object_CD189, y = c(Object_CD299, Object_CD364), add.cell.ids = c('CD189', 'CD299', 'CD364'), project = 'GSE164985')
HT_merge<- merge(Object_HT206, y = c(Object_HT214, Object_HT216, Object_HT217), add.cell.ids = c('HT206', 'HT214', 'HT216', 'HT217'), project = 'GSE164985')
CD_merge[["CellName"]] <- colnames(CD_merge)
HT_merge[["CellName"]] <- colnames(HT_merge)

HT_colonocyte <- readRDS('~/cybergut_project/Case_Studies/GSE164985/HT_VS_CD_colonocyte/HT_colonocyte.rds')
CD_colonocyte <- readRDS('~/cybergut_project/Case_Studies/GSE164985/HT_VS_CD_colonocyte/CD_colonocyte.rds')

Name_CD_colonocyte <- colnames(CD_colonocyte)
Name_HT_colonocyte <- colnames(HT_colonocyte)

HT_colonocyte <- subset(HT_merge, subset = CellName %in% Name_HT_colonocyte )
CD_colonocyte <- subset(CD_merge, subset = CellName %in% Name_CD_colonocyte )

CD_colonocyte@meta.data$Health = 'Inflammed'
HT_colonocyte@meta.data$Health = 'Healthy'

colonocyte_merge  <- merge(CD_colonocyte, y = HT_colonocyte)
pseudo_ifnb <- AggregateExpression(colonocyte_merge, assays = "RNA", return.seurat = T, group.by = c("Health", "Sample"))
DEG = FindMarkers(pseudo_ifnb, ident.1 = 'Inflammed', ident.2 = "Healthy",test.use = "DESeq2")

