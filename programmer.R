library(Seurat)
library(patchwork)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(fishpond)
library(tximport)
library(Matrix)
library(tidyverse)
library(EnsDb.Hsapiens.v79)


# read in alevin data
files <- file.path("salmon_output/salmon_output_combined/alevin/quants_mat.gz")

# import using txi
txi <- tximport(files, type="alevin")
t <- txi$counts


# get ensembl ids, and remove extensions
ensembl <-  rownames(t) #get ensmbl ids
ensembl <- as.character(genes)
ensembl <- sub("[.][0-9]*","",ensembl) #remove ends
rownames(t) <- ensembl

# convert gene names of ensembl ids 
symbols <- select(EnsDb.Hsapiens.v79, keys= ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID")) 

# subset to rows that have the matches, and rename matrix with genenames
t <- t[rownames(t) %in% symbols$GENEID,]
rownames(t) <- symbols$SYMBOL


#-------------Seurat -----------------------
# create object
pbmc <- CreateSeuratObject(counts = t, min.cells = 3, min.features = 200, project = "10X_PBMC")

# Quality control - mt genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

png("vlnplot.png")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


# Quality control - view features 
png("scatter.png", width = 800)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# get rid of low quality either high (explain why you chose)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4200 & percent.mt < 20 )

# normalize
pbmc <- NormalizeData(pbmc)


# feature selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
 
# plot variable features with and without labels
png("feature_selection.png", width = 1000)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

# scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Linear Reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
png("pca.png")
DimPlot(pbmc, reduction = "pca")
dev.off()

# determine residuals
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

png("jackstraw.png")
JackStrawPlot(pbmc, dims = 1:15)
dev.off()

png("elbow.png")
ElbowPlot(pbmc) 
dev.off()

# Cluster
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


# Run non-linear reduction
pbmc <- RunUMAP(pbmc, dims = 1:10)
png("umap.png")
DimPlot(pbmc, reduction = "umap")
dev.off()

# save
saveRDS(pbmc, file = "pbmc_prelim.rds")

# plot counts, proportions
counts <- as.vector(table(Idents(pbmc)))
names <- names(table(Idents(pbmc)))

png('barplot.png')
bar <- barplot(height = table(Idents(pbmc)), names.arg = (names), xlab = 'Cluster', ylab = 'total cells')
text(x = bar, y = table(Idents(pbmc))-60, label = table(Idents(pbmc)), pos = 3, cex = 0.8, col = "blue")
dev.off()

png('pie.png')
slices <- as.vector(prop.table(table(Idents(pbmc))))
lbls <- names(prop.table(table(Idents(pbmc))))
pie(slices, labels = lbls, main="relative proportions")
dev.off()
