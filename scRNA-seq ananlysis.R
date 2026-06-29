#------------single sample-----------

data_dir <- '/path_data/data'
counts<-readMM("matrix.mtx.gz")
barcodes = fread("barcodes.tsv.gz")
features = fread("features.tsv.gz")
expression_matrix <- Read10X(data.dir = data_dir, gene.column = 2, cell.column = 1, unique.features = TRUE)
seurat_object1 = CreateSeuratObject(counts = expression_matrix, project = "sample1", min.cells = 3, min.features = 200)
seurat_object1[["percent.mt"]] <- PercentageFeatureSet(seurat_object1, pattern = "^MT-")
seurat_object1[["percent.rb"]] <- PercentageFeatureSet(seurat_object1, pattern = "^RP[SL]")
seurat_object1[["percent.hb"]] <- PercentageFeatureSet(seurat_object1, pattern = "^HB[^(P)]")
VlnPlot(seurat_object1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb", "percent.hb"), ncol = 5)

seurat_object1 <- subset(seurat_object1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15 & percent.rb > 0 & percent.hb < 0.2)
seurat_object1 <- NormalizeData(seurat_object1)
seurat_object1 <- FindVariableFeatures(seurat_object1)
seurat_object1 <- ScaleData(seurat_object1)
seurat_object1 <- RunPCA(seurat_object1)
stdv <- seurat_object1[["pca"]]@stdev
sum.stdv <- sum(seurat_object1[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min.pc <- min(co1, co2)
seurat_object1 <- RunUMAP(seurat_object1, dims = 1:min.pc)
seurat_object1 <- FindNeighbors(object = seurat_object1, dims = 1:min.pc)
seurat_object1 <- FindClusters(object = seurat_object1, resolution = 0.1)
sweep.list <- paramSweep(seurat_object1, PCs = 1:min.pc)
sweep.stats <- summarizeSweep(sweep.list)
bcmvn <- find.pK(sweep.stats)
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
annotations <- seurat_object1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp.poi <- round(ncol(seurat_object1) * 0.04)
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
seurat_object1 <- doubletFinder(seu = seurat_object1, 
                                PCs = 1:min.pc,
                                pN = 0.25,
                                pK = optimal.pk,
                                nExp = nExp.poi.adj)
DF.name = colnames(seurat_object1@meta.data)[grepl("DF.classification", colnames(seurat_object1@meta.data))]
cowplot::plot_grid(ncol = 2, DimPlot(seurat_object1, group.by = "orig.ident") + NoAxes(),
                   DimPlot(seurat_object1, group.by = DF.name) + NoAxes())
VlnPlot(seurat_object1, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

seurat_object1 = seurat_object1[,seurat_object1@meta.data[, DF.name] == "Singlet"]
save(seurat_object1,file='seurat_object1.Rdata')

#------------merge-----------
load("/path_to_data/seurat_object1.Rdata")
load("/path_to_data/seurat_object2.Rdata")
load("/path_to_data/seurat_object3.Rdata")
load("/path_to_data/seurat_object4.Rdata")
load("/path_to_data/seurat_object5.Rdata")
load("/path_to_data/seurat_object6.Rdata")

ifnb.list <- list(seurat_object1, seurat_object2, seurat_object3, seurat_object4,
                  seurat_object5, seurat_object6)

ifnb <- Merge_Seurat_List(ifnb.list, add.cell.ids = c("sample1", "sample2", "sample3", "sample4",
                                                      "sample5", "sample6"), merge.data = TRUE, project = "SeuratProject")
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)

ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:40)
ifnb <- FindClusters(ifnb, resolution = 0.4)
ifnb <- RunUMAP(ifnb, dims = 1:40, reduction = "integrated.cca")

#------------subpopulations-----------


load("/path_to_data/ifnb.RData")
sce <- ifnb
Idents(sce) = "celltype"
Idents(sce)

myogenic_sce <- subset(sce, idents="Myogenic cells")

myogenic_sce[["RNA"]] <- split(
  myogenic_sce[["RNA"]],
  f = myogenic_sce$orig.ident
)

myogenic_sce <- NormalizeData(myogenic_sce)
myogenic_sce <- FindVariableFeatures(myogenic_sce)
myogenic_sce <- ScaleData(myogenic_sce)
myogenic_sce <- RunPCA(myogenic_sce)

myogenic_sce <- IntegrateLayers(
  myogenic_sce,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
)

myogenic_sce[["RNA"]] <- JoinLayers(myogenic_sce[["RNA"]])

ElbowPlot(myogenic_sce, ndims = 50)

myogenic_sce <- FindNeighbors(
  myogenic_sce,
  reduction="integrated.cca",
  dims=1:40
)

myogenic_sce <- FindClusters(myogenic_sce, resolution=0.4)

myogenic_sce <- RunUMAP(
  myogenic_sce,
  reduction="integrated.cca",
  dims=1:40
)



