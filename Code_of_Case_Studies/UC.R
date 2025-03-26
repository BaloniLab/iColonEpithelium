library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratWrappers)
library(Matrix)
library(Biobase)

## Manually load ata
library(Matrix)
library(Seurat)
epi.counts = readMM('/scratch/negishi/jiang817/CyberGut_project/Case_studies/Broad_study/data/gene_sorted-Epi.matrix.mtx')
rownames(epi.counts) = readLines('/scratch/negishi/jiang817/CyberGut_project/Case_studies/Broad_study/data/Epi.genes.tsv')
colnames(epi.counts) = readLines('/scratch/negishi/jiang817/CyberGut_project/Case_studies/Broad_study/data/Epi.barcodes2.tsv')
meta = read.table('/scratch/negishi/jiang817/CyberGut_project/Case_studies/Broad_study/data/all.meta2.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)

# load seurat object 
batch_epi <- readLines('/scratch/negishi/jiang817/CyberGut_project/Case_studies/Broad_study/data/Epi.barcodes2.tsv')
meta_epi <- meta[batch_epi, ]
seur = CreateSeuratObject(counts=epi.counts, meta.data = meta_epi)

# EXTEACT individual samples
health_sample_id <- meta_epi[meta_epi$Health == 'Healthy',]
health_sample_id <- unique(health_sample_id$Sample)

for (i in health_sample_id){
  filename = paste0('Health_', i, '.rsd')
  individual_sample <- subset( x = seur, subset = Sample == i)
  save_path = paste0('~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Healthy/', filename)
  saveRDS(individual_sample, save_path)
}

inflame_sample_id <- meta_epi[meta_epi$Health == 'Inflamed',]
inflame_sample_id <- unique(inflame_sample_id$Sample)

for (i in inflame_sample_id){
  filename = paste0('Inflame_', i, '.rsd')
  individual_sample <- subset( x = seur, subset = Sample == i)
  save_path = paste0('~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Inflamed/', filename)
  saveRDS(individual_sample, save_path)
}


## impute data 
# health_sample
setwd("~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Healthy")
All_health_samples <- list.files('~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Healthy/')
All_health_samples <- All_health_samples[1:24]
for (i in All_health_samples){
  object <- readRDS(i)
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object <- subset(object, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- RunALRA(object)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 5000)
  filepath <- paste0('~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Healthy/Imputed_sample_data/', 'Imputed_',i)
  saveRDS(object, filepath)
}

# inflame_sample
setwd("~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Inflamed")
All_inflame_samples <- list.files('~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Inflamed/')
All_inflame_samples <- All_inflame_samples[2:17]
for (i in All_inflame_samples){
  object <- readRDS(i)
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object <- subset(object, subset = nFeature_RNA > 1000  & nFeature_RNA < 50000 & percent.mt < 25)
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- RunALRA(object)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 5000)
  filepath <- paste0('~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Inflamed/Imputed_sample_data/', 'Imputed_',i)
  saveRDS(object, filepath)
}

# extract rawdata
#Healtjh
setwd("~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Healthy/Imputed_sample_data")
health_imputed_files <- list.files('~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Healthy/Imputed_sample_data/')
for (i in health_imputed_files){
  filename = strsplit(i, '.rsd')[[1]]
  object <- readRDS(i)
  object <- try(subset(object, subset = Cluster == 'Enterocytes'))
  counts <- GetAssayData(object)
  savepath <- paste0('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/Counts_sets/Imputed_Health/', filename, '_counts.csv')
  write.csv(counts, savepath)
}

#Inflame
setwd("~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Inflamed/Imputed_sample_data")
inflame_imputed_files <- list.files('~/cybergut_project/Case_Studies/Broad_study/data/individual_sample_data/Inflamed/Imputed_sample_data/')
for (i in inflame_imputed_files){
  filename = strsplit(i, '.rsd')[[1]]
  object <- readRDS(i)
  object <- try(subset(object, subset = Cluster == 'Enterocytes'))
  counts <- GetAssayData(object)
  savepath <- paste0('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/Counts_sets/Imputed_Inflame/', filename, '_counts.csv')
  write.csv(counts, savepath)
}

# convert symbol to BIGG
# HEALTH
setwd("~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/Counts_sets/Imputed_Health/Update_id")
symbol2entrez <- read.csv('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/symtol2entrez.csv')
health_count_files = list.files('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/Counts_sets/Imputed_Health/Update_id/')
for (i in health_count_files){
  exprs <- read.csv(i, header = 1)
  exprs_count <- exprs[, 2:ncol(exprs)]
  Row_median = as.data.frame(rowMedians(as.matrix(exprs_count)))
  colnames(Row_median) = 'Medium_Value'
  Row_median['hgnc_symbol'] <- exprs[, 1]
  expressionset <- merge(x = Row_median, y = symbol2entrez, by.x = 'hgnc_symbol', all.x =TRUE, no.dups = TRUE )
  suffix  = '_AT1'
  BIGGID <- paste0(expressionset$entrezgene_id, suffix)
  expressionset['BIGG_ID'] = BIGGID
  filename <- substr(i, 9, nchar(i) -10)
  savepath <- paste0('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/expressionSets/Health/', filename, 'expression.csv')
  write.csv(expressionset, savepath, row.names = FALSE)
}


#INFLAME
setwd("~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/Counts_sets/Imputed_Inflame/Update_id")
symbol2entrez <- read.csv('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/symtol2entrez.csv')
inflame_count_files = list.files('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/Counts_sets/Imputed_Inflame/Update_id/')
for (i in inflame_count_files){
  exprs <- read.csv(i, header = 1)
  exprs_count <- exprs[, 2:ncol(exprs)]
  Row_median = as.data.frame(rowMedians(as.matrix(exprs_count)))
  colnames(Row_median) = 'Medium_Value'
  Row_median['hgnc_symbol'] <- exprs[, 1]
  expressionset <- merge(x = Row_median, y = symbol2entrez, by.x = 'hgnc_symbol', all.x =TRUE, no.dups = TRUE )
  suffix  = '_AT1'
  BIGGID <- paste0(expressionset$entrezgene_id, suffix)
  expressionset['BIGG_ID'] = BIGGID
  filename <- substr(i, 9, nchar(i) -10)
  savepath <- paste0('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/expressionSets/Inflame/', filename, 'expression.csv')
  write.csv(expressionset, savepath, row.names = FALSE)
}


# t-test 
library(dplyr)
data <- read.csv('~/cybergut_project/Case_Studies/Broad_study/FAV_analysis/Fva_outputs/UC_FVA_TEST.csv', row.names = 1)
data <- data %>% 
  mutate(P_value=NA)

for(i in 1: 5087){
  print(rownames(data)[i])
  fluxes <- data[i,1:36]

  fluxes <- as.data.frame(t(fluxes))
  colnames(fluxes) <- c('value')
  fluxes$group <- c('HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT', 'HT',
                    'UC', 'UC', 'UC', 'UC', 'UC', 'UC', 'UC', 'UC', 'UC', 'UC', 'UC', 'UC', 'UC')
  fluxes$group <- factor(fluxes$group)


  result <- try(t.test(value~group, fluxes, paired = FALSE, alternative = 'two.sided'))
  if('try-error' %in% class(result)){
    p_value <- "NAN"
    data[i,37] <- p_value
    
  }else{
    p_value <- result$p.value
    data[i,37] <- p_value
  }
  
}

    
write.csv(data, 't_test_outputs.csv')



