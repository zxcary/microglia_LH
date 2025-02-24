setwd("./")
res_dir <- "."
library(Seurat)
set.seed(100)
Idents(seurat_obj) <- "RNA_snn_res.1"
seurat_anno <- RenameIdents(object = seurat_obj,
                            "0"="Microgila",
                            "1"="Microgila",
                            "2"="Microgila",
                            "3"="Microgila",
                            "4"="Microgila",
                            "5"="Microgila",
                            "6"="Microgila",
                            "7"="Microgila",
                            "8"="Oligodendrocyte precursor cell",
                            "9"="Microgila",
                            "10"="BAM",
                            "11"="Astrocyte",
                            "12"="Microgila",
                            "13"="BAM",
                            "14"="BAM",
                            "15"="Microgila",
                            "16"="T",
                            "17"="Neutrophil",
                            "18"="Endothelial cell",
                            "19"="Late activated neural stem cell",
                            "20"="DC",
                            "21"="B",
                            "22"="Monocyte",
                            "23"="Smooth muscle cell"
                            
)

seurat_anno$anno_rough <- seurat_anno@active.ident
DimPlot(object = seurat_anno,reduction = "tsne",label = T)
saveRDS(object = seurat_anno,file = "seurat_anno.rds")
seurat_anno <- readRDS(file = "seurat_anno.rds")


############ 细分t cell

t_cell <- subset(seurat_anno,idents = "T")

DefaultAssay(t_cell) <- "RNA"
t_cell = NormalizeData(t_cell)
t_cell = FindVariableFeatures(t_cell,nfeatures = 2000)%>% ScaleData()


t_cell = RunPCA(t_cell,verbose = F)
t_cell <- RunHarmony(t_cell, group.by.vars = "sample",plot_convergence = T)

t_cell <-  t_cell%>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = c(0.4,0.5,0.6,0.8,1,1.2,1.4,1.5))
t_cell <- RunTSNE(t_cell,reduction = "harmony",dims = 1:30)
t_cell <- RunUMAP(t_cell,reduction = "harmony",dims = 1:30)

Idents(object = t_cell) <- "RNA_snn_res.0.2"
DimPlot(object = t_cell,reduction = "tsne",label = T)
t_cell <- RenameIdents(object = t_cell,"0"="T","1"="T","2"="ILC","3"="T")

markers_t <- FindAllMarkers(object = t_cell,assay = "RNA",logfc.threshold = 0.5,only.pos = T)


VlnPlot(object = t_cell,features = c("Cd3e","Cd3d","Mki67"),split.by = "sample")
FeaturePlot(object = t_cell,features = c("Cd3d"),reduction = "tsne")
FeaturePlot(object = seurat_anno,features = c("Cd3e"),reduction = "tsne")

data_info <- data.frame(t_cell@active.ident)
data_info[["cell_id"]] <- rownames(data_info)
data_info <- data.frame(t(data_info[data_info$t_cell.active.ident!="T cell",]))
data_meta <- data.frame(seurat_anno@active.ident)
data_meta[["cell_id"]] <- rownames(data_meta)
data_meta <- data.frame(t(data_meta))
data_meta[,colnames(data_info)] <- data_info
new_anno <- factor(unlist(data_meta[1,]))
names(new_anno) <- unlist(data_meta[2,])
seurat_anno@active.ident <- new_anno
seurat_anno$anno_detail <- new_anno

################# 细分 neutrophil
neut <- subset(seurat_anno,idents = c("Neutrophil"))

DefaultAssay(neut) <- "RNA"
neut = NormalizeData(neut)
neut = FindVariableFeatures(neut,nfeatures = 2000)%>% ScaleData()


neut = RunPCA(neut,verbose = F)
neut <- RunHarmony(neut, group.by.vars = "sample",plot_convergence = T)

neut <-  neut%>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    Findclusters(resolution = c(0.4,0.5,0.6,0.8,1,1.2,1.4,1.5))
neut <- RunTSNE(neut,reduction = "harmony",dims = 1:30)
neut <- RunUMAP(neut,reduction = "harmony",dims = 1:30)

Idents(object = neut) <- "RNA_snn_res.0.5"

DimPlot(object = neut,reduction = "tsne",label = T)
DimPlot(object = neut,reduction = "tsne",label = F,split.by = "re_names",ncol = 2)

DefaultAssay(object = neut) <- "RNA"
t_markers <- FindAllMarkers(object = neut,assay = "RNA",logfc.threshold = 0.25,min.pct = 0.1)
t_markers <- t_markers[t_markers$p_val_adj<0.05,]
t_markers <- t_markers[order(t_markers$cluster,-t_markers$avg_log2FC),]

top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5,(length(y)-4):length(y))]
})))
top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:10)]
})))
library(MySeuratWrappers)
DefaultAssay(object = neut) <- "RNA"
VlnPlot(neut,features = top_5,assay = "RNA",stacked = T,pt.size = 0)
FeaturePlot(object = neut,features = c("Siglech","Ccr7","Ly6c2","Axl","Ccr9","Xcr1"),
            label = T,reduction = "tsne")
VlnPlot(neut,features = c("Siglech","Ccr7","Ly6c2","Ly6d","Axl","Ccr9","Xcr1"),
        assay = "RNA",stacked = T,pt.size = 0)

neut_anno <- RenameIdents(neut,"0"="Neu_1","1"="Neu_1","2"="Neu_2","3"="Neu_2")
neut_anno[["anno_rough"]] <- neut_anno@active.ident
# Idents(neut_anno) <- "RNA_snn_res.0.5"
markers_neu <- FindAllMarkers(object = neut_anno,assay = "RNA",logfc.threshold = 0.5,only.pos = T)


data_info <- data.frame(neut_anno$anno_rough)
data_info[["cell_id"]] <- rownames(data_info)
data_info <- data.frame(t(data_info))
# data_info <- data.frame(t(data_info[data_info$neut_anno.active.ident!="Neutrophil",]))
data_meta <- data.frame(seurat_anno$anno_rough)
data_meta[["cell_id"]] <- rownames(data_meta)
data_meta <- data.frame(t(data_meta))
data_meta[,colnames(data_info)] <- data_info
new_anno <- factor(unlist(data_meta[1,]))
names(new_anno) <- unlist(data_meta[2,])
seurat_anno@active.ident <- new_anno
seurat_anno$anno_detail <- new_anno
DimPlot(object = seurat_anno,reduction = "tsne",label = T)

################# 细分 DC
dc <- subset(seurat_anno,idents = c("DC"))

DefaultAssay(dc) <- "RNA"
dc = NormalizeData(dc)
dc = FindVariableFeatures(dc,nfeatures = 2000)%>% ScaleData()


dc = RunPCA(dc,verbose = F)
dc <- RunHarmony(dc, group.by.vars = "sample",plot_convergence = T)

dc <-  dc%>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = c(0.4,0.5,0.6,0.8,1,1.2,1.4,1.5))
dc <- RunTSNE(dc,reduction = "harmony",dims = 1:30)
dc <- RunUMAP(dc,reduction = "harmony",dims = 1:30)

Idents(object = dc) <- "RNA_snn_res.0.5"

DimPlot(object = dc,reduction = "tsne",label = T)
DimPlot(object = dc,reduction = "tsne",label = F,split.by = "re_names",ncol = 2)

DefaultAssay(object = dc) <- "RNA"
t_markers <- FindAllMarkers(object = dc,assay = "RNA",logfc.threshold = 0.25,min.pct = 0.1)
t_markers <- t_markers[t_markers$p_val_adj<0.05,]
t_markers <- t_markers[order(t_markers$cluster,-t_markers$avg_log2FC),]

top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5,(length(y)-4):length(y))]
})))
top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5)]
})))
library(MySeuratWrappers)
DefaultAssay(object = dc) <- "RNA"
VlnPlot(dc,features = top_5,assay = "RNA",stacked = T,pt.size = 0)
FeaturePlot(object = dc,features = c("Siglech","Ccr7","Ly6c2","Axl","Itgam","Xcr1"),
            label = T,reduction = "tsne")
VlnPlot(dc,features = c("Siglech","Ccr7","Ly6c2","Ly6d","Axl","Ccr9","Xcr1","Flt3"),
        assay = "RNA",stacked = T,pt.size = 0)

dc_anno <- RenameIdents(dc,"0"="cDC","1"="cDC","2"="migDC")
dc_anno[["anno_rough"]] <- dc_anno@active.ident
# Idents(dc_anno) <- "RNA_snn_res.0.5"
markers_dc <- FindAllMarkers(object = dc_anno,assay = "RNA",logfc.threshold = 0.5,only.pos = T)


data_info <- data.frame(dc_anno$anno_rough)
data_info[["cell_id"]] <- rownames(data_info)
data_info <- data.frame(t(data_info))
# data_info <- data.frame(t(data_info[data_info$dc_anno.active.ident!="dcrophil",]))
data_meta <- data.frame(seurat_anno$anno_rough)
data_meta[["cell_id"]] <- rownames(data_meta)
data_meta <- data.frame(t(data_meta))
data_meta[,colnames(data_info)] <- data_info
new_anno <- factor(unlist(data_meta[1,]))
names(new_anno) <- unlist(data_meta[2,])
seurat_anno@active.ident <- new_anno
seurat_anno$anno_detail <- new_anno
DimPlot(object = seurat_anno,reduction = "tsne",label = T)


################# 细分 mono
mono <- subset(seurat_anno,idents = c("Monocyte"))

DefaultAssay(mono) <- "RNA"
mono = NormalizeData(mono)
mono = FindVariableFeatures(mono,nfeatures = 2000)%>% ScaleData()


mono = RunPCA(mono,verbose = F)
mono <- RunHarmony(mono, group.by.vars = "sample",plot_convergence = T)

mono <-  mono%>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    Findclusters(resolution = c(0.4,0.5,0.6,0.8,1,1.2,1.4,1.5))
mono <- RunTSNE(mono,reduction = "harmony",dims = 1:30)
mono <- RunUMAP(mono,reduction = "harmony",dims = 1:30)

Idents(object = mono) <- "RNA_snn_res.0.5"

DimPlot(object = mono,reduction = "tsne",label = T)
DimPlot(object = mono,reduction = "tsne",label = F,split.by = "re_names",ncol = 2)

DefaultAssay(object = mono) <- "RNA"
t_markers <- FindAllMarkers(object = mono,assay = "RNA",logfc.threshold = 0.25,min.pct = 0.1)
t_markers <- t_markers[t_markers$p_val_adj<0.05,]
t_markers <- t_markers[order(t_markers$cluster,-t_markers$avg_log2FC),]

top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5,(length(y)-4):length(y))]
})))
top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5)]
})))
library(MySeuratWrappers)
DefaultAssay(object = mono) <- "RNA"
VlnPlot(mono,features = top_5,assay = "RNA",stacked = T,pt.size = 0)
FeaturePlot(object = mono,features = c("Siglech","Cd74","Cd38","Malat1","Neat1","Lyve1"),
            label = T,reduction = "tsne")
FeaturePlot(object = mono,features = c("Ace","Ly6c2","Fcgr1","Hp","Ly6c2","F10"),
            label = T,reduction = "tsne")

VlnPlot(mono,features = c("Siglech","Cd74","Ccr7","Ly6c2","Ly6d","Axl","Xcr1","Flt3"),
        assay = "RNA",stacked = T,pt.size = 0)

mono_anno <- RenameIdents(mono,"0"="c_Mono","1"="nc_Mono")
mono_anno[["anno_rough"]] <- mono_anno@active.ident
# Idents(mono_anno) <- "RNA_snn_res.0.5"
markers_mono <- FindAllMarkers(object = mono_anno,assay = "RNA",logfc.threshold = 0.5,only.pos = T)


data_info <- data.frame(mono_anno$anno_rough)
data_info[["cell_id"]] <- rownames(data_info)
data_info <- data.frame(t(data_info))
# data_info <- data.frame(t(data_info[data_info$mono_anno.active.ident!="monorophil",]))
data_meta <- data.frame(seurat_anno$anno_rough)
data_meta[["cell_id"]] <- rownames(data_meta)
data_meta <- data.frame(t(data_meta))
data_meta[,colnames(data_info)] <- data_info
new_anno <- factor(unlist(data_meta[1,]))
names(new_anno) <- unlist(data_meta[2,])
seurat_anno@active.ident <- new_anno
seurat_anno$anno_detail <- new_anno
DimPlot(object = seurat_anno,reduction = "tsne",label = T)


################# 细分 B
b <- subset(seurat_anno,idents = c("B cell"))

DefaultAssay(b) <- "RNA"
b = NormalizeData(b)
b = FindVariableFeatures(b,nfeatures = 2000)%>% ScaleData()


b = RunPCA(b,verbose = F)
b <- RunHarmony(b, group.by.vars = "sample",plot_convergence = T)

b <-  b%>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    Findclusters(resolution = c(0.4,0.5,0.6,0.8,1,1.2,1.4,1.5))
b <- RunTSNE(b,reduction = "harmony",dims = 1:30)
b <- RunUMAP(b,reduction = "harmony",dims = 1:30)

Idents(object = b) <- "RNA_snn_res.0.5"

DimPlot(object = b,reduction = "tsne",label = T)
DimPlot(object = b,reduction = "tsne",label = F,split.by = "re_names",ncol = 2)

DefaultAssay(object = b) <- "RNA"
t_markers <- FindAllMarkers(object = b,assay = "RNA",logfc.threshold = 0.25,min.pct = 0.1)
t_markers <- t_markers[t_markers$p_val_adj<0.05,]
t_markers <- t_markers[order(t_markers$cluster,-t_markers$avg_log2FC),]

top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5,(length(y)-4):length(y))]
})))
top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5)]
})))
library(MySeuratWrappers)
DefaultAssay(object = b) <- "RNA"
VlnPlot(b,features = top_5,assay = "RNA",stacked = T,pt.size = 0)
FeaturePlot(object = seurat_anno,features = c("Ms4a1","Mzb1","Igkc","Ly6d","Cd79b"),
            label = F,reduction = "tsne")
FeaturePlot(object = seurat_anno,features = c("Ly6c2"),
            label = T,reduction = "tsne")


VlnPlot(b,features = c("Siglech","Cd74","Ccr7","Ly6c2","Ly6d","Axl","Xcr1","Flt3"),
        assay = "RNA",stacked = T,pt.size = 0)

VlnPlot(seurat_anno,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        assay = "RNA",stacked = T,pt.size = 0)



b_anno <- RenameIdents(b,"0"="B_1","1"="B_1","2"="B_2")
b_anno[["anno_rough"]] <- b_anno@active.ident
# Idents(b_anno) <- "RNA_snn_res.0.5"


data_info <- data.frame(b_anno$anno_rough)
data_info[["cell_id"]] <- rownames(data_info)
data_info <- data.frame(t(data_info))
# data_info <- data.frame(t(data_info[data_info$b_anno.active.ident!="brophil",]))
data_meta <- data.frame(seurat_anno$anno_rough)
data_meta[["cell_id"]] <- rownames(data_meta)
data_meta <- data.frame(t(data_meta))
data_meta[,colnames(data_info)] <- data_info
new_anno <- factor(unlist(data_meta[1,]))
names(new_anno) <- unlist(data_meta[2,])
seurat_anno@active.ident <- new_anno
seurat_anno$anno_detail <- new_anno
DimPlot(object = seurat_anno,reduction = "tsne",label = F)





################################################################## bam
bam <- subset(immu_seurat,idents = c("BAM"))

DefaultAssay(bam) <- "RNA"
bam = NormalizeData(bam)
bam = FindVariableFeatures(bam,nfeatures = 2000)%>% ScaleData()


bam = RunPCA(bam,verbose = F)
bam <- RunHarmony(bam, group.by.vars = "sample",plot_convergence = T)

bam <-  bam%>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    Findclusters(resolution = c(0.4,0.5,0.6,0.8,1,1.2,1.4,1.5))
bam <- RunTSNE(bam,reduction = "harmony",dims = 1:30)
bam <- RunUMAP(bam,reduction = "harmony",dims = 1:30)

Idents(object = bam) <- "RNA_snn_res.0.5"

DimPlot(object = bam,reduction = "tsne",label = T)
DimPlot(object = bam,reduction = "umap",label = T)

DimPlot(object = bam,reduction = "tsne",label = F,split.by = "re_names",ncol = 2)
DimPlot(object = bam,reduction = "umap",label = F,split.by = "re_names",ncol = 2)



DefaultAssay(object = bam) <- "RNA"
t_markers <- FindAllMarkers(object = bam,assay = "RNA",logfc.threshold = 0.25,min.pct = 0.1)
t_markers <- t_markers[t_markers$p_val_adj<0.05,]
t_markers <- t_markers[order(t_markers$cluster,-t_markers$avg_log2FC),]
t_markers25 <- FindMarkers(object = bam,assay = "RNA",ident.1 = "2",ident.2 = "5",logfc.threshold = 0.25,min.pct = 0.1)
t_markers02 <- FindMarkers(object = bam,assay = "RNA",ident.1 = "0",ident.2 = "2",logfc.threshold = 0.25,min.pct = 0.1)

write.csv(x = t_markers,"bam_markers.csv")


top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5,(length(y)-4):length(y))]
})))
top_5 <- unique(unlist(lapply(X = unique(t_markers$cluster),FUN = function(x){
    y <- t_markers[t_markers$cluster==x,"gene"]
    y <- y[c(1:5)]
})))
library(MySeuratWrappers)

VlnPlot(bam,features = top_5,assay = "RNA",stacked = T,pt.size = 0)
VlnPlot(bam,features =  c("Apoe","Ccl5","Rpl13","Pik3ip1","Rsrp1","Cd52","P2ry12","Malat1","Gm26917"),stacked = T,pt.size = 0)
VlnPlot(bam,features =  c("Tpt1","Rps14","Rps7","Cd52","Ccl5","Ifi27l2a","Apoe"),stacked = T,pt.size = 0)

FeaturePlot(object = bam,features = c("P2ry12","Rnase4","Slc2a5","Aif1","Cd52","Rsrp1"),
            label = T,reduction = "umap")
FeaturePlot(object = bam,features = c("Pf4"),
            label = T,reduction = "tsne")

aa <- AverageExpression(object =bam_anno,features = c("Cbr2","Wfdc17","Pf4","Malat1","Rsrp1") )


# 查看某一簇
res_rate <- "RNA_snn_res.0.5"
cluster_num <- "8"
cellscd4 <- unique(bam@assays[["RNA"]]@data@Dimnames[[2]][bam@meta.data[[res_rate]] == cluster_num])
DimPlot(bam, reduction = "umap", cells.highlight = cellscd4 ,raster=FALSE,label = TRUE)


VlnPlot(object = bam,assay = "RNA",features = gene0,stacked = T,pt.size = 0)
VlnPlot(object = bam,assay = "RNA",features = c("Cd3d","Cd3e","Cd3g","Cd8a","Il1rl1","Rnf128","Il2ra","Icos","Gata3","Emb","Ltb4r1","Pard3"),pt.size = 0)

VlnPlot(object = seurat_RNA,assay = "RNA",features = c("Apoe","Ms4a7"),pt.size = 0)


#######################注释

bam_anno <- RenameIdents(bam,"0"="BAM1_Pf4hi_Malat1lo","1"="BAM2_Pf4hi_Malat1hi_Lyve1lo","2"="BAM3_Cd74hi_Sparclo",
                         "3"="BAM4_Pf4hi_Malat1hi_Lyve1hi","4"="BAM5_Pf4lo_P2ry12hi","5"="BAM6_Cd74hi_Sparchi")
bam_anno[["anno_full"]] <- bam_anno@active.ident
Idents(bam_anno) <- "RNA_snn_res.0.5"
bam_anno <- RenameIdents(bam,"0"="BAM1","1"="BAM2","2"="BAM3",
                         "3"="BAM4","4"="BAM5","5"="BAM6")
bam_anno[["anno_abbr"]] <- bam_anno@active.ident
bam_anno$re_names <- factor(bam_anno$re_names,levels = c( "Con","LPS","HYP","L/H"))
saveRDS(object = bam_anno,file = "bam_anno.rds")
DimPlot(object = bam_anno,reduction = "tsne",label = T)

Idents(bam_anno) <- "anno_abbr"

FeaturePlot(object = bam_anno,features = c("P2ry12","Cd52","Ifitm3","Gm26917"),ncol = 2,cols = c("yellow","red"),reduction = "umap",label = T)

###################3
bam_genes <- c("Pf4","Malat1","Cd74","Lyve1","P2ry12","Sparc")
Idents(bam_anno) <- "anno_abbr"
VlnPlot(object = bam_anno,assay = "RNA",features = bam_genes,pt.size = 0,stacked = T)
all.genes <- rownames(bam_anno)
bam_anno <- ScaleData(bam_anno, features = all.genes)

DoHeatmap(bam_anno,features = bam_genes,assay = "RNA") 

ave_exp_markers <- AverageExpression(bam_anno,features = bam_genes,assay = "RNA")

library(pheatmap)
bk <- c(seq(-1,-0.1,by=0.02),seq(0,1,by=0.02))
pheatmap(mat = ave_exp_markers$RNA,scale = "row",
         cluster_cols = F,angle_col = 45,cluster_rows = F,
         breaks = bk,main = "Average Expression of Marker Genes")


##############333高亮图

res_rate <- "anno_abbr"
cluster_num1 <- "BAM3"
cluster_num2 <- "BAM4"
cluster_num3 <- "BAM5"
cluster_num4 <- "BAM6"

cells_highlight <- list("BAM3" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num1]),
                        "BAM4" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num2]),
                        "BAM5" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num3]),
                        "BAM6" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num4])
)
cells_highlight <- list("BAM3" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num1])
                        
)
DimPlot(object = bam_anno,split.by = "re_names",
        cells.highlight = cells_highlight,
        sizes.highlight = 0.5,ncol = 4,reduction = "tsne",
        cols.highlight = c("#20b2aa","#009900"))

DimPlot(object = bam_anno,split.by = "re_names",
        cells.highlight = cells_highlight,
        sizes.highlight = 0.5,ncol = 4,reduction = "tsne",
        cols.highlight = c("#ff00ff","#6495ed","#20b2aa","#009900"))

cells_highlight <- list("BAM3" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num1]))

DimPlot(object = bam_anno,split.by = "re_names",
        cells.highlight = cells_highlight,
        sizes.highlight = 0.5,reduction = "tsne",
        cols.highlight = c("BAM3"="#00BA38"))

cells_highlight <- list("BMA4" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num2]))

DimPlot(object = bam_anno,split.by = "re_names",
        cells.highlight = cells_highlight,
        sizes.highlight = 0.5,reduction = "tsne",
        cols.highlight = c("BAM4"="#20b2aa"))

cells_highlight <- list("BMA5" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num3]))
DimPlot(object = bam_anno,split.by = "re_names",
        cells.highlight = cells_highlight,
        sizes.highlight = 0.5,reduction = "tsne",
        cols.highlight = c("BAM4"="#6495ed"))

cells_highlight <- list("BMA6" = unique(bam_anno@assays[["RNA"]]@data@Dimnames[[2]][bam_anno@meta.data[[res_rate]] == cluster_num4]))
DimPlot(object = bam_anno,split.by = "re_names",
        cells.highlight = cells_highlight,
        sizes.highlight = 0.5,reduction = "tsne",
        cols.highlight = c("BAM4"="#ff00ff"))

###############比例统计
anno_info <- data.frame(sample=bam_anno$re_names,
                        cell_type=bam_anno@active.ident)

cell_number <- data.frame(table(anno_info[c("sample","cell_type")]))
cell_num0 <- data.frame(table(anno_info[c("sample")]))
cell_number[["sample_num"]] <- rep(cell_num0$Freq,length(unique(bam_anno$anno_abbr)))
cell_number[["percentage"]] <- round(cell_number$Freq/cell_number$sample_num*100,2)
cell_number$sample <- factor(cell_number$sample,levels = names(col_model))

cell_names <- unique(cell_number$cell_type)
write.csv(x = cell_number,file = "percentage_in_bam.csv")

library(scater)
list_pics <- list()
for (cell_name in cell_names){
    cell_number1 <- cell_number[cell_number$cell_type==cell_name,]
    # pdf(file = paste0(cell_name,".pdf"),width =5,height = 5)
    list_pics[[cell_name]] <- ggplot(data = cell_number1,mapping = aes(x = sample,y = percentage,color=sample)) +
        geom_bar(stat = "identity",mapping = aes(x=sample,y = percentage,fill=sample))+
        geom_text(mapping = aes(label = percentage,y = percentage),vjust = 0,position = position_dodge(0.9),size=3) +
        # geom_label(mapping = aes(colour = "black",label = percentage)) +
        # theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,margin = margin(t = 0,r = -3,b = 0)),) +
        scale_colour_manual(values = alpha(colour = col_model,alpha = 0.8)) + 
        scale_fill_manual(values = alpha(colour = col_model,alpha = 0.8)) +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_rect(colour = "black",fill = NA),
              plot.title = element_text(hjust = 0.5,size = 8))+
        ylab(label = "Percentage") + 
        # ylim(0,100)+
        ggtitle(paste0(cell_name," in bam cells"))
}


multiplot(plotlist = list_pics,cols=3)

#############差异基因
library(ggplot2)
library(ggrepel)

Idents(object = bam_anno) <- "anno_abbr"
deg_bam3 <- FindMarkers(object = bam_anno,ident.1 = "BAM3",assay = "RNA",logfc.threshold = 0,min.pct = 0.1)
# deg_bam3 <- deg_bam3[deg_bam3$p_val_adj<0.05,]
deg_bam4 <- FindMarkers(object = bam_anno,ident.1 = "BAM4",assay = "RNA",logfc.threshold = 0,min.pct = 0.1)
# deg_bam4 <- deg_bam4[deg_bam4$p_val_adj<0.05,]
deg_bam5 <- FindMarkers(object = bam_anno,ident.1 = "BAM5",assay = "RNA",logfc.threshold = 0,min.pct = 0.1)
# deg_bam5 <- deg_bam5[deg_bam5$p_val_adj<0.05,]
deg_bam6 <- FindMarkers(object = bam_anno,ident.1 = "BAM6",assay = "RNA",logfc.threshold = 0,min.pct = 0.1)
# deg_bam6 <- deg_bam6[deg_bam6$p_val_adj<0.05,]



write.csv(deg_bam3,"deg_bam3.csv")
write.csv(deg_bam4,"deg_bam4.csv")
write.csv(deg_bam3,"deg_bam5.csv")
write.csv(deg_bam4,"deg_bam6.csv")


Dat <- deg_bam6
#确定是上调还是下调，用于给图中点上色）
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 0.25, ifelse(Dat$avg_log2FC>= 0.25 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
Dat$Gene = rownames(Dat)
for (i in 1:nrow(Dat)){
    if(Dat[i,"p_val_adj"]<10**(-300)){
        Dat[i,"p_val_adj"] = 10**(-(300+sample(x = 1:20,size = 1)))
    }
}
dd <- Dat$p_val_adj<10**(-300)

ggplot(Dat,aes(x=avg_log2FC,y=-log10(p_val_adj),color=threshold))+
    geom_point()+
    scale_color_manual(values=c("Up"="#CC0033", "Down"="#3366CC","NoSignifi"="#808080"))+#确定点的颜色
    geom_text_repel(
        data = Dat[Dat$p_val_adj<0.05&abs(Dat$avg_log2FC)>0.25,],
        aes(label = Gene),
        box.padding = 0.01,
        point.padding = 0,
        seed = 243,
        size = 3,
        segment.color = "black", show.legend = FALSE )+#添加关注的点的基因名
    theme_bw()+#修改图片背景
    theme(
        legend.title = element_text("BAM3")#不显示图例标题
    )+
    ylab('-log10 (p_val_adj)')+#修改y轴名称
    xlab('avg_log2FC')+#修改x轴名称
    geom_vline(xintercept=c(-0.25,0.25),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
    geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05


######################富集分析


library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)

deg_bam3 <- read.csv()
deg_list <- list(deg_bam3,deg_bam4,deg_bam5,deg_bam6)
names(deg_list) <- c("BAM3","BAM4","BAM5","BAM6")
for (i in names(deg_list)){
    # deg_list[[i]] <- deg_list[[i]][deg_list[[i]]$avg_log2FC>0.25,]
    deg_list[[i]]$gene <- rownames(deg_list[[i]])
    deg_list[[i]]$cluster <- rep(x = i,nrow(deg_list[[i]]))
}

markers <- do.call(what = rbind,args = deg_list)

list_pic_count <- list()
list_pic_qvalue <- list()

for (i in unique(markers$cluster)){
    genes0 <- markers[markers$cluster==i & markers$avg_log2FC>0.25 & markers$p_val_adj < 0.05,"gene"]
    genes1 <- markers[markers$cluster==i & markers$avg_log2FC< -0.25 & markers$p_val_adj < 0.05,"gene"]
    enrich_go_bp_up <- enrichGO(gene = genes0,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP")
    enrich_go_bp_down <- enrichGO(gene = genes1,OrgDb = org.Mm.eg.db,keyType = "SYMBOL",ont = "BP")
    go_int <- AnnotationDbi::select(GO.db,keys=enrich_go_bp_up$ID,columns = c("DEFINITION"),keytype = "GOID")
    enrich_go_bp_up@result <- merge(x = enrich_go_bp_up@result,y=go_int,by.x="ID",by.y="GOID",sort = F)
    go_int <- AnnotationDbi::select(GO.db,keys=enrich_go_bp_down$ID,columns = c("DEFINITION"),keytype = "GOID")
    enrich_go_bp_down@result <- merge(x = enrich_go_bp_down@result,y=go_int,by.x="ID",by.y="GOID",sort = F)
    
    write.csv(x = enrich_go_bp_up,file = paste0("go_up_bp_enrich_t_",i,".csv"))
    write.csv(x = enrich_go_bp_down,file = paste0("go_down_bp_enrich_t_",i,".csv"))
    
    
    ############################# count
    tmp <- enrich_go_bp_down
    tmp@result <- tmp@result[1:15,]
    tmp@result <- tmp@result[order(tmp@result$Count,decreasing = T,method = "radix"),]
    tmp@result$Count <- -tmp@result$Count
    
    enrich_trans <- enrich_go_bp_up
    
    enrich_trans@result <- enrich_trans@result[1:15,]
    enrich_trans@result <- enrich_trans@result[order(enrich_trans@result$Count,decreasing = F,method = "radix"),]
    # enrich_trans@result$Count <- -log10(enrich_trans@result$p.adjust)
    
    ####重复展示
    int_err <- intersect(x = tmp@result$Description,y = enrich_trans@result$Description)
    
    raw_des <- enrich_trans@result$Description
    for (err0 in int_err){
        raw_des <- gsub(pattern = err0,replacement = paste0(err0,"_1"),x = raw_des)
    }
    enrich_trans@result$Description <- raw_des
    
    raw_des <- tmp@result$Description
    for (err0 in int_err){
        raw_des <- gsub(pattern = err0,replacement = paste0(err0,"_2"),x = raw_des)
    }
    tmp@result$Description <- raw_des
    
    enrich_trans@result <- rbind(tmp@result,enrich_trans@result)
    
    enrich_trans@result[["regulation"]] <- rep(c("down","up"),c(15,15))
    aa <- enrich_trans@result
    aa$Description <- factor(x = aa$Description,levels = unique(aa$Description))
    
    p <- ggplot(aa, aes(x = Count, y = Description, fill = regulation)) +
        scale_fill_manual(values = c("up"="#CC0033", "down"="#3366CC"))+
        geom_col() + # geom_bar(stat = "identity") + coord_flip() +
        ggtitle(label = i)+
        labs(x="Count",y=NULL) + 
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              # axis.text.x = element_text(angle = 90, hjust = 1),
              # axis.text.y = element_text(angle = 90, hjust = 1),
              panel.border = element_rect(colour = "black",fill = NA),
              plot.title = element_text(hjust = 0.5))  
    # xlim(c(-15,35))
    
    # xlim(c(-15,35))
    list_pic_count[[i]] <- p
    
    ############################# pvalue
    tmp <- enrich_go_bp_down
    tmp@result <- tmp@result[1:15,]
    tmp@result <- tmp@result[order(tmp@result$p.adjust,decreasing = F,method = "radix"),]
    tmp@result$Count <- log10(tmp@result$p.adjust)
    
    enrich_trans <- enrich_go_bp_up
    
    enrich_trans@result <- enrich_trans@result[1:15,]
    enrich_trans@result <- enrich_trans@result[order(enrich_trans@result$p.adjust,decreasing = T,method = "radix"),]
    enrich_trans@result$Count <- -log10(enrich_trans@result$p.adjust)
    
    ####重复展示
    int_err <- intersect(x = tmp@result$Description,y = enrich_trans@result$Description)
    
    raw_des <- enrich_trans@result$Description
    for (err0 in int_err){
        raw_des <- gsub(pattern = err0,replacement = paste0(err0,"_1"),x = raw_des)
    }
    enrich_trans@result$Description <- raw_des
    
    raw_des <- tmp@result$Description
    for (err0 in int_err){
        raw_des <- gsub(pattern = err0,replacement = paste0(err0,"_2"),x = raw_des)
    }
    tmp@result$Description <- raw_des
    
    enrich_trans@result <- rbind(tmp@result,enrich_trans@result)
    
    enrich_trans@result[["regulation"]] <- rep(c("down","up"),c(15,15))
    aa <- enrich_trans@result
    aa$Description <- factor(x = aa$Description,levels = unique(aa$Description))
    
    p <- ggplot(aa, aes(x = Count, y = Description, fill = regulation)) +
        scale_fill_manual(values = c("up"="#CC0033", "down"="#3366CC"))+
        geom_col() + # geom_bar(stat = "identity") + coord_flip() +
        ggtitle(label = i) +
        labs(x="-log10(p.adjust)",y=NULL) + 
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              # axis.text.x = element_text(angle = 90, hjust = 1),
              # axis.text.y = element_text(angle = 90, hjust = 1),
              panel.border = element_rect(colour = "black",fill = NA),
              plot.title = element_text(hjust = 0.5))  
    # xlim(c(-15,35))
    list_pic_qvalue[[i]] <- p
    
}

list_pic_count$BAM3
list_pic_count$BAM4
list_pic_count$BAM5
list_pic_count$BAM6


list_pic_qvalue$BAM3
list_pic_qvalue$BAM4
list_pic_qvalue$BAM5
list_pic_qvalue$BAM6

#########################
bam3 <- subset(x = bam_anno,idents = c("BAM3"))
DefaultAssay(bam3) <- "RNA"
genes_all <- rownames(bam3@assays$RNA)
bam3 <- ScaleData(object = bam3,features = genes_all)
Idents(bam3) <- "re_names"
markers_4 <- FindAllMarkers(object = bam3,assay = "RNA",logfc.threshold = 0.25,min.pct = 0.1)
markers_4 <- markers_4[markers_4$p_val_adj<0.05,]
markers_4_up <- markers_4[markers_4$avg_log2FC>0,]
markers_4_up <- markers_4_up[order(markers_4_up$avg_log2FC,decreasing = T),]
markers_4_down <- markers_4[markers_4$avg_log2FC<0,]
markers_4_down <- markers_4_down[order(markers_4_down$avg_log2FC,decreasing = F),]

top_5_up <- unique(unlist(lapply(X = c("Con","LPS","HYP","L/H"),FUN = function(x){
    y <- markers_4_up[markers_4_up$cluster==x,"gene"]
    y <- y[c(1:10)]
})))
top_5_down <- unique(unlist(lapply(X = c("Con","LPS","HYP","L/H"),FUN = function(x){
    y <- markers_4_down[markers_4_down$cluster==x,"gene"]
    y <- y[c(1:5)]
})))
DoHeatmap(object = bam3,assay = "RNA",features = top_5_up,group.by = "re_names")
exp_ave <- AverageExpression(object = bam3,assays = "RNA",features = top_5_up)
DoHeatmap(object = bam3,assay = "RNA",features = top_5_down,group.by = "re_names")
pheatmap(mat = exp_ave$RNA,scale = "row",cluster_cols = F,cluster_rows = F)

markers_4_up_list <- list()
for (i in unique(markers_4_up$cluster)){
    markers_4_up_list[[i]] <- rownames(markers_4_up[markers_4_up$cluster==i,])
}
markers_4_down_list <- list()
for (i in unique(markers_4_down$cluster)){
    markers_4_down_list[[i]] <- rownames(markers_4_down[markers_4_down$cluster==i,])
}
res_up <- compareCluster(geneClusters = markers_4_up_list,fun = enrichGO,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont="BP")
res_down <- compareCluster(geneClusters = markers_4_down_list,fun = enrichGO,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont="BP")

dotplot(res_up)
dotplot(res_down)
write.csv(res_up,"./res_up.csv")
write.csv(res_down,"./res_down.csv")



##############################3
#################3
seurat_anno <- readRDS(file = "./data/seurat_anno_final.rds")
Idents(seurat_anno) <- "anno_detail"
DimPlot(object = seurat_anno,reduction = "tsne",label = T)
Idents(seurat_anno) <- "anno_plot"
DimPlot(object = seurat_anno,reduction = "tsne",label = T)
col_seq <- c("#9489fa","#f47a75","#96d7f9","#5690dd","#bd88f5","#765005","#f7af59","#f0da49","#71c16f","#024b51","#009db2","#f06464","#0780cf","blue","yellow","green")
DimPlot(object = seurat_anno,reduction = "tsne",label = T,cols = col_seq)

##################################### clean
# 去除B_2粘连细胞
seurat_anno_clean <- subset(x = seurat_anno,idents = c("B_2"),invert= TRUE)
unique(seurat_anno_clean$anno_detail)

# 重新命名
tmp <- seurat_anno_clean$anno_detail
# tmp <- gsub(pattern = "cDC",replacement = "DC",x = tmp)
tmp <- gsub(pattern = "B_1",replacement = "B",x = tmp)
# tmp <- gsub(pattern = "Mono_1",replacement = "c_Mono",x = tmp)
# tmp <- gsub(pattern = "Mono_2",replacement = "nc_Mono",x = tmp)
# tmp <- gsub(pattern = "T cell",replacement = "T",x = tmp)
# tmp <- gsub(pattern = "Neutrophil_1",replacement = "Neu_1",x = tmp)
# tmp <- gsub(pattern = "Neutrophil_2",replacement = "Neu_2",x = tmp)



seurat_anno_clean$anno_detail <- tmp
Idents(seurat_anno_clean) <- "anno_detail"

DimPlot(seurat_anno_clean,reduction = "tsne",label = T)
FeaturePlot(object = mono,features = c("Ly6c2","C1qa"),
            label = T,reduction = "tsne")




# Idents(object = seurat_anno) <- "anno_rough"
# seurat_anno$anno_rough2 <- gsub(pattern = "NKT cell",replacement = "Cycling T",seurat_anno$anno_rough)

list_dic <- list("Y-1"="Con","Y-2"="LPS","Y-3"="HYP","Y-4"="L/H")
re_names <- seurat_anno_clean$sample
re_names2 <- c()
for (i in re_names){
    re_names2 <- c(re_names2,list_dic[[i]])
}
names(re_names2) <- rownames(seurat_anno_clean@meta.data)
seurat_anno_clean[["re_names"]] <- re_names2

cell_number <- data.frame(table(seurat_anno_clean$anno_detail))
cell_number <- cell_number[order(cell_number$Freq,decreasing = T),]
write.csv(cell_number,paste0(res_dir,"cell_number.csv"))
seurat_anno_clean$anno_detail <- factor(x = seurat_anno_clean$anno_detail,levels = cell_number$Var1)
saveRDS(object = seurat_anno_clean,file = paste0(res_dir, "seurat_anno_clean.rds"))
Idents(seurat_anno_clean) <- "anno_detail"
col_seq <- c("#9489fa","#f06464","#f47a75","#5690dd","#f7af59","#f0da49","#71c16f","#2aaaef","#96d7f9","#bd88f5","#009db2","#024b51","#0780cf")

DimPlot(object = seurat_anno_clean,reduction = "umap",label = T,cols = col_seq)
DimPlot(object = seurat_anno_clean,reduction = "tsne",label = T)
DimPlot(object = seurat_anno_clean,reduction = "tsne",label = F)
DimPlot(object = seurat_anno_clean,reduction = "tsne",label = F,split.by = "re_names",ncol = 2)
DimPlot(object = seurat_anno_clean,reduction = "umap",label = F,split.by = "re_names",ncol = 2)
DefaultAssay(seurat_anno_clean) <- "RNA"

FeaturePlot(object = seurat_anno_clean,features = c("Mbp","Mobp","Slc1a2","Sparcl1","Cldn5","Hist1h2ap","Top2a","Cenpf","Rgs5"),ncol = 3,cols = c("yellow","red"),reduction = "tsne")
FeaturePlot(object = seurat_anno_clean,features = c("Mbp","Slc1a2","Cldn5","Top2a","Rgs5"),ncol = 3,cols = c("yellow","red"),reduction = "tsne")


################免疫细胞

seurat_anno <- readRDS("./seurat_anno_clean.rds")
seurat_anno_clean <- readRDS("./seurat_anno_clean.rds")
# col_seq2 <- c("#9489fa","#f06464","#f7af59","#f0da49","#96d7f9","#bd88f5","#009db2","#0780cf")
genes_other <- c("Mbp","Mobp","Slc1a2","Sparcl1","Cldn5","Rgs5","Hist1h2ap","Top2a","Cenpf")
FeaturePlot(object = seurat_obj,features =genes_other,ncol = 3,cols = c("yellow","red"),reduction = "tsne")

genes_other <- c("Mbp","Slc1a2","Cldn5","Top2a","Rgs5")
FeaturePlot(object = seurat_obj,features =genes_other,ncol = 3,cols = c("yellow","red"),reduction = "tsne")



FeaturePlot(object = seurat_anno,features ="Rgs5",cols = c("yellow","red"),reduction = "tsne")

DimPlot(object = seurat_anno,reduction = "tsne",label = T)

bam <- subset(x = seurat_anno,idents = "BAM")
Idents(object = bam) <- "anno_detail"
DimPlot(object = bam,reduction = "tsne",label = T)
DefaultAssay(bam) <- "RNA"
DEGs_bam <- FindAllMarkers(object = bam,assay = "RNA",logfc.threshold = 0.5)
DEGs_bam <- DEGs_bam[DEGs_bam$p_val_adj<0.05,]
DEGs_bam <- DEGs_bam[order(DEGs_bam$cluster,-DEGs_bam$avg_log2FC),]

######
col_seq_immu <- c("#9489fa","#f47a75","#5690dd","#f7af59","#f0da49","#71c16f","#96d7f9","#bd88f5","#009db2","#024b51","#765005")

immu_seurat <- subset(seurat_anno_clean,idents = c("Smooth muscle cell","Endothelial cell","Astrocyte",
                                                   "Oligodendrocyte precursor cell",
                                                   "Late activated neural stem cell"),invert = T)
immu_seurat <- subset(seurat_anno,idents = c("Smooth muscle cell","Endothelial cell","Astrocyte",
                                             "Oligodendrocyte precursor cell",
                                             "Late activated neural stem cell"),invert = T)


# col_seq_immu <- c("#9489fa","#f47a75","#5690dd","#009db2","#f0da49","#71c16f","#96d7f9","#bd88f5","#765005","#024b51","#f7af59")
cell_types <- c("Microgila","BAM","Neu_1","Neu_2","T","ILC","cDC","migDC","B","c_Mono","nc_Mono")
col_seq_immu <- c("#9489fa","#f47a75","#96d7f9","#5690dd","#bd88f5","#765005","#f7af59","#f0da49","#71c16f","#024b51","#009db2")
immu_seurat$detail_anno <- factor(immu_seurat$anno_detail,levels = cell_types)
Idents(immu_seurat) <- "detail_anno"
DimPlot(object = immu_seurat,reduction = "tsne",label = T,cols = col_seq_immu)
DimPlot(object = immu_seurat,reduction = "tsne",label = F,cols = col_seq_immu)

# Idents(immu_seurat) <- "anno_rough"
DimPlot(object = immu_seurat,reduction = "tsne",split.by = "re_names",label = F,ncol = 2,cols = col_seq_immu)

################## 小胶之外

immu_seurat2 <- subset(immu_seurat,idents = c("Microgila"), invert=T)
# list_dic <- list("Y-1"="Con","Y-2"="LPS","Y-3"="HYP","Y-4"="L/H")
# re_names <- immu_seurat2$sample
# re_names2 <- c()
# for (i in re_names){
#     re_names2 <- c(re_names2,list_dic[[i]])
# }
# names(re_names2) <- rownames(immu_seurat2@meta.data)
# immu_seurat2[["re_names"]] <- re_names2


DefaultAssay(object = immu_seurat2) <- "RNA"
ElbowPlot(immu_seurat2,ndims = 50)

immu_seurat2 <- FindNeighbors(immu_seurat2,dims = 1:30) #2mins
immu_seurat2 <- FindClusters(immu_seurat2,
                             resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6,1.8,2))
# immu_seurat2 <- FindClusters(immu_seurat2,resolution = c(0.6))

immu_seurat2 <- RunTSNE(immu_seurat2,dims = 1:30)
immu_seurat2 <- RunUMAP(immu_seurat2,dims = 1:30)



cell_number <- data.frame(table(immu_seurat2$anno_detail))
cell_number <- cell_number[order(cell_number$Freq,decreasing = T),]

cell_types2 <- c("BAM","Neu_1","Neu_2","T","ILC","cDC","migDC","B","c_Mono","nc_Mono")

immu_seurat2$anno_detail <- factor(x = immu_seurat2$anno_detail,levels = cell_types2)

Idents(immu_seurat2) <- "anno_detail"
# col_seq3 <- c("#f06464","#f7af59","#f0da49","#96d7f9","#bd88f5","#009db2","#0780cf")
# col_seq_immu2 <- c("#f47a75","#5690dd","#009db2","#f0da49","#71c16f","#96d7f9","#bd88f5","#765005","#024b51","#f7af59")
col_seq_immu2 <- c("#f47a75","#96d7f9","#5690dd","#bd88f5","#765005","#f7af59","#f0da49","#71c16f","#024b51","#009db2")

DimPlot(object = immu_seurat2,reduction = "tsne",label = T,cols = col_seq_immu2)

DimPlot(object = immu_seurat2,reduction = "tsne",label = F,split.by = "re_names",ncol = 2,cols = col_seq_immu2)
# DimPlot(object = immu_seurat2,reduction = "umap",label = F,split.by = "re_names",ncol = 2)

############################ markers
res_dir <- ""
immu_markers <- FindAllMarkers(object = immu_seurat,logfc.threshold = 0.25,assay = "RNA")
immu_markers <- immu_markers[immu_markers$p_val_adj<0.05,]
immu_markers2 <- immu_markers[order(immu_markers$avg_log2FC,decreasing = T),]
immu_markers2 <- immu_markers2[order(immu_markers2$cluster,decreasing = T),]
write.csv(x = immu_markers2,file = paste(res_dir,"DEGs_immu_markers.csv"))
immu_markers2 <- read.csv(file = "./DEGs_immu_markers.csv",header = T,row.names = 1)
heatmap_genes <- unique(unlist(lapply(X = unique(immu_markers2$cluster),FUN = function(x){
    y <- immu_markers2[immu_markers2$cluster==x,"gene"]
    y <- y[c(1:5)]
})))
heatmap_genes <- c(heatmap_genes,c("F10","C1qa","Siglech"))
genes_all <- rownames(immu_seurat@assays$RNA)
immu_seurat <- ScaleData(object = immu_seurat,features = heatmap_genes,assay = "RNA")
pdf(file = "./Doheatmap.pdf",width = 8,height = 8)
Seurat::DoHeatmap(immu_seurat,features = heatmap_genes,assay = "RNA",group.by = "detail_anno")+scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

p_data <- data.frame(immu_seurat@assays$RNA@scale.data[heatmap_genes,])

cells_id <- data.frame(rownames(immu_seurat@meta.data),immu_seurat$detail_anno)
cells_id$immu_seurat.detail_anno <- factor(x = cells_id$immu_seurat.detail_anno,levels = levels(immu_seurat$detail_anno))
p_data <- p_data[,order(cells_id$immu_seurat.detail_anno)]
bk <- c(seq(-1,-0.1,by=0.02),seq(0,1,by=0.02))
gaps_col <- c(5,10,17,22,27,36,41)
pheatmap(mat = p_data,scale = "none",gaps_row = gaps_col,labels_col = F,
         cluster_cols = F,angle_col = 45,cluster_rows = F,
         breaks = bk,main = "Average Expression of Marker Genes",axis.text.y.left = element_text())

ave_exp <- AverageExpression(object = immu_seurat, assays = "RNA",features = heatmap_genes,group.by = "anno_detail")
library(pheatmap)
p_data <- ave_exp$RNA
p_data <- data.frame(t(p_data[,cell_types]))
bk <- c(seq(-1,-0.1,by=0.02),seq(0,1,by=0.02))
gaps_col <- c(5,10,17,22,27,36,41)
pheatmap(mat = p_data,scale = "column",gaps_col = gaps_col,
         cluster_cols = F,angle_col = 45,cluster_rows = F,
         breaks = bk,main = "Average Expression of Marker Genes",axis.text.y.left = element_text())

pheatmap(mat = p_data,scale = "row",cluster_rows = F,cluster_cols = F)
library(MySeuratWrappers)
DefaultAssay(immu_seurat) <- "RNA"

VlnPlot(immu_seurat,features = top_5,stacked = T,pt.size = 0)
VlnPlot(immu_seurat,features = c("H2-Aa"),pt.size = 0)

cell_number <- data.frame(table(immu_seurat$anno_detail))
cell_number <- cell_number[order(cell_number$Freq,decreasing = T),]
write.csv(cell_number,"cell_number.csv")


list_markers <- list("Microgila" = c("Tmem119","Sparc","Hexb"),
                     # "oligodendrocyte precursor cell" = c("Mbp","Mobp"),
                     # "Bergmann glial cell" = c("Slc1a2","Sparcl1"),
                     "BAM" = c("Pf4","Apoe","Ms4a7"),
                     "Neutrophils" = c("S100a9","S100a8","Retnlg","Ltf"),
                     "T"=c("Cd3d","Cd8a"),# 
                     "ILC" = c("Gata3","Il1rl1","Rnf128"),
                     "DC" = c("Itgax","Flt3","Cd74"),
                     "B"=c("Igkc","Ly6d","Cd79b"),
                     "Monocytes" = c("Hp","Ly6c2","F10")
                     # "Cycling T cell" = c("Mki67","Gzmb"),
                     # "Endothelial cell"=c("Cldn5"),
                     # "Type II spiral ganglion neuron"=c("H2-Aa"),#
                     # "Late activated neural stem cell" = c("Hist1h2ap","Top2a","Cenpf"),
                     # "Smooth muscle cell" = c("Rgs5")
)

list_markers <- list("Microgila" = c("Tmem119","Sparc","Hexb"),
                     "BAM" = c("Apoe","Pf4","Ms4a7"),
                     "Neutrophils" = c("S100a9","S100a8","Retnlg","Il1b","Ccl6","Camp","Ltf"),
                     "T"=c("Cd3d","Cd8a"),#
                     "ILC" = c("Gata3","Il1rl1","Rnf128"),
                     "cDC" = c("Itgax","Flt3","Cd74","Gm2a","Ciita"),
                     "migDC"= c("Ccr7","Ccl22"),
                     "B"=c("Igkc","Ly6d","Cd79b"),
                     "Monocytes" = c("Hp","Ly6c2","F10","Vcan","S100a4","C1qa","Siglech")
)
genes0 <- unlist(list_markers)
names(genes0) <- NULL
cell_types <- c("Microgila","BAM","Neu_1","Neu_2","T","ILC","cDC","migDC","B","c_Mono","nc_Mono")

DefaultAssay(immu_seurat) <- "RNA"
immu_seurat$anno_rough2 <- factor(x = immu_seurat$anno_rough,levels = cell_types[length(cell_types):1])
Idents(immu_seurat) <- "anno_rough2"
DotPlot(immu_seurat,features =genes0, cols = c("blue","#CC0033"),
        dot.scale = 8,idents = cell_types) + RotatedAxis() 

FeaturePlot(object = immu_seurat,features = c("Tmem119","Pf4","S100a9","Cd3d","Cd74","Cd79b","Gata3","Mki67"),ncol = 3,cols = c("yellow","red"),reduction = "tsne")
FeaturePlot(object = immu_seurat,features = c("Sparc"),ncol = 1,cols = c("yellow","red"),reduction = "tsne")

DefaultAssay(immu_seurat2) <- "RNA"
FeaturePlot(object = immu_seurat2,features = c("Apoe","Cd74","Mgl2","Hp","S100a9","Cd3d","Gata3","Cd79b"),ncol = 4,cols = c("yellow","red"),reduction = "tsne")
FeaturePlot(object = immu_seurat2,features = c("Trbc2","Gzmb","Ms4a4b","Nkg7"),cols = c("yellow","red"),reduction = "tsne")
FeaturePlot(object = immu_seurat2,features = c("Apoe","Ccl6","Ltf","Cd3d","Gata3","Gm2a","Ccr7","Cd79b","S100a4","C1qa"),ncol = 5,cols = c("yellow","red"),reduction = "tsne")
FeaturePlot(object = immu_seurat2,features = c("Apoe","S100a8","Cd3d","Gata3","Cd74","Cd79b","F10"),ncol = 4,cols = c("yellow","red"),reduction = "umap")


FeaturePlot(object = immu_seurat2,features = c("Pf4","Apoe"),ncol = 2,cols = c("yellow","red"),reduction = "tsne")



#################################### cell percentage
cell_number <- data.frame(table(immu_seurat$anno_detail,immu_seurat$re_names) )
# cell_number <- cell_number[order(cell_number$Freq,decreasing = T),]
colnames(cell_number) <- c("Cell_Type","Group","Freq")
cell_others <- c("Oligodendrocyte precursor cell", "Astrocyte", "Late activated neural stem cell", "Endothelial cell", "Smooth muscle cell")

cell_number <- cell_number[!cell_number$Cell_Type %in% cell_others,]
cell_types <- c("Microgila","BAM","Neu_1","Neu_2","T","ILC","cDC","migDC","B","c_Mono","nc_Mono")
cell_number$Cell_Type <- factor(x = cell_number$Cell_Type,levels = cell_types)
group_order <- c("Con","LPS","HYP","L/H")
cell_number$Group <- factor(x = cell_number$Group,levels = group_order)

cell_number[["Cell_Number"]] <- rep(table(immu_seurat$re_names),rep(length(cell_types),length(group_order)))
cell_number[["Group_Number"]] <- rep(table(cell_number$Group),rep(length(cell_types),length(group_order)))

cell_number[["Percentage"]] <- round(cell_number[["Freq"]]/cell_number[["Cell_Number"]]*100,2)
cell_number[["Percentage2"]] <- round(cell_number[["Percentage"]]/cell_number[["Group_Number"]],2)
write.csv(x = cell_number,file = "cell_number_detail.csv")

library(gg.gap)
p1 <- ggplot(data = cell_number,mapping = aes(x = Group,y = Percentage,fill=Cell_Type)) +
    geom_col() +  #coord_flip() +
    scale_fill_manual(values=col_seq_immu) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(hjust = 1,angle = 45),
          panel.border = element_rect(colour = "white",fill = NA),
          plot.title = element_text(hjust = 0.5,size = 8))
p2 <- gg.gap(plot = p1,segments = c(15,90),ylim = c(0,100))
p2

