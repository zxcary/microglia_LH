setwd("~/microgila/")
library(Seurat)
library(DoubletFinder)
library(dplyr)
library(scater)
library(ggplot2)

input_dir = "~/microgila/EXP2/"
des_dir <- "./"
sample_names <- c("Y-1","Y-2","Y-3","Y-4" )
sample_renames <- c("Con","LPS","HYP","L/H" )
group_names <-  c("Con","LPS","HYP","L_H" )

hiv.list0 <- lapply(X = 1:length(sample_names),FUN = function(i) {
    file_path <- file.path(input_dir, sample_names[i],"outs/filtered_feature_bc_matrix")
    print(file_path)
    d10x <- Read10X(file_path)
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x), split='-'), '[[', 1L), sample_names[i], sep='-')
    pbmc <- CreateSeuratObject(counts = d10x, project = sample_names[i], min.cells = ncol(d10x)*0.001, min.features = 200)
    pbmc[["sample"]] <- rep(sample_renames[i],nrow(pbmc@meta.data))
    pbmc[["group"]] <- rep(group_names[i],nrow(pbmc@meta.data))
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
    pbmc[["log10GenesPerUMI"]] <- log10(pbmc[["nFeature_RNA"]]) / log10(pbmc[["nCount_RNA"]])
    pbmc
})

names(hiv.list0) <- sample_names
plot_list1 <- lapply(hiv.list0,function(x){
    VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0)
})
plot_list1[[4]]

hiv.list_clean <- list()
dim <- 30 # 维度
res <- 1 # 分辨率
pK_bcmvn <- 0.02
list_d <- list()
for (i in names(hiv.list0)){
    # hiv.list, fun = function(x){
    x <- hiv.list0[[i]]
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
    x <- ScaleData(x, verbose = T) # a few seconds
    x <- RunPCA(x,npcs = dim, verbose = T)  # 1 min
    # ElbowPlot(x,ndims = 50)
    x <- FindNeighbors(x,dims = 1:dim) #2mins
    x <- FindClusters(x,resolution = c(res))
    # x$seurat_clusters
    x <- RunUMAP(x, dims = 1:dim) # 1 min
    #x <- RunTSNE(x, dims = 1:dim) # 1 min
    sweep.res.list <- paramSweep_v3(x, PCs = 1:dim, sct = F)
    #使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
    bcmvn <- find.pK(sweep.stats) #可以看到最佳参数的点
    #pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
    
    DoubletRate = ncol(x)*5*1e-6 #更通用
    #估计同源双细胞比例，根据modelHomotypic()中的参数人为混合双细胞。这里是从seurat_clusters中来混双细胞 
    homotypic.prop <- modelHomotypic(x$seurat_clusters) #最好提供celltype，而不是seurat_clusters。
    # 计算双细胞比例
    nExp_poi <- round(DoubletRate*ncol(x)) 
    # 使用同源双细胞比例对计算的双细胞比例进行校正 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## 使用确定好的参数鉴定doublets
    x <- doubletFinder_v3(x, PCs = 1:dim, pN = 0.25, pK = pK_bcmvn, 
                          nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    
    x$doubFind_res = x@meta.data %>% select(contains('DF.classifications'))
    x$doubFind_score = x@meta.data %>% select(contains('pANN'))
    hiv.list_clean[[i]] <- x
    list_d[[i]] <- bcmvn
}

saveRDS(object = list_d,file = "list_d.rds")
saveRDS(object = hiv.list_clean,file = "brain_list_clean.rds")

list_dd <- list()
for (i in names(hiv.list_clean)){
    x <- hiv.list_clean[[i]]
    ff <- data.frame(table(x@meta.data %>% select(contains('doubFind_res'))))
    ff <- data.frame("Freq"=c(ff[,2],ff[1,2]/c(ff[1,2]+ff[2,2])),row.names = c(ff[,1],"percentage"))
    colnames(ff) <- i
    list_dd[[i]] <- ff
}
res_d <- do.call(what = cbind,args = list_dd)
res_d

hiv.list <- lapply(X = hiv.list_clean,FUN = function(x){
    pbmc <- subset(x, subset = nFeature_RNA >= 500  & nCount_RNA <= 20000 & (nCount_RNA >=1000) & (percent.mt <= 10))
})
names(hiv.list) <- sample_renames

# cell numbers count
cell_numbers0 <- sapply(hiv.list0,FUN = function(x)dim(x)[2])
cell_numbers1 <- sapply(hiv.list,FUN = function(x)dim(x)[2])
cell_numbers2 <- data.frame(Sample_names=rep(sample_names,2),
                            Sample_renames=rep(sample_renames,2),
                            Condition=rep(group_names,2),
                            Group=factor(rep(c("raw","clean"),c(length(cell_numbers0),length(cell_numbers1))),levels = c("raw","clean")),
                            Cell_numbers=c(cell_numbers0,cell_numbers1)
)
cell_numbers2$Sample_renames <- factor(x = cell_numbers2$Sample_renames,levels = sample_renames)
write.csv(cell_numbers2,file = paste0(des_dir,"cell_numbers.csv"))

pdf(file = "quality_cell_number.pdf",width = 8,height = 6)
ggplot(data = cell_numbers2,mapping = aes(x= Sample_renames,y = Cell_numbers,group = Group,fill = Group))+
    geom_col(position = position_dodge(width=0.9))+
    geom_text(mapping = aes(y= Cell_numbers+1000,label = Cell_numbers),position = position_dodge(width=0.9))+
    theme_bw()
dev.off()

pdf("quality_sample.pdf",width = 6,height =5)
for (x in hiv.list){
    p <- VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0)
    print(p)
}
dev.off()

sce.all <- merge(hiv.list[[1]], 
                 y = c(hiv.list[[2]],hiv.list[[3]],hiv.list[[4]]),
                 add.cell.ids = sample_names, #添加样本名
                 project = "scRNA")


seurat_obj = NormalizeData(sce.all)
seurat_obj = FindVariableFeatures(seurat_obj,nfeatures = 2000)%>% ScaleData()


seurat_obj = RunPCA(seurat_obj,verbose = F)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "sample",plot_convergence = T)

seurat_obj <-  seurat_obj%>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = c(0.2,0.4,0.5,0.6,0.8,1,1.2,1.4,1.5,1.6,1.8,2))


seurat_obj <- seurat_obj %>% RunUMAP(reduction = "harmony", dims = 1:30)
seurat_obj <- seurat_obj %>% RunTSNE(reduction = "harmony", dims = 1:30)

#Idents(seurat_obj) <- "RNA_snn_res.1"
Idents(seurat_obj) <- "doubFind_res"
pdf(file = paste0(des_dir,"dimplot_doubFind.pdf"),width =10,height = 8 )
DimPlot(seurat_obj,label = T,reduction = "umap",raster = F)
dev.off()
Idents(seurat_obj) <- "RNA_snn_res.1"
pdf(file = paste0(des_dir,"dimplot_cluster_res.1.pdf"),width =10,height = 8 )
DimPlot(seurat_obj,label = T,reduction = "umap",raster = F)
dev.off()
saveRDS(object = seurat_obj,file = paste0(des_dir,"brain_obj_dim30.rds"))
#Idents(seurat_obj) <- "RNA_snn_res.0.5"
#all_markers_cluster <- FindAllMarkers(object = seurat_obj,assay = "RNA",logfc.threshs = T)

list_doublet <- list()
#i <- "34"
for (i in unique(seurat_obj$RNA_snn_res.1)){
    T_PRKCH <- subset(x = seurat_obj,idents = i)
    list_doublet[[i]] <- c(i,table(T_PRKCH$doubFind_res))
}
#list_doublet[["28"]] <- c(28,0,348)
for (i in names(list_doublet)){
    tt <- length(list_doublet[[i]])
    if(tt ==2 ){
        list_doublet[[i]] <- c(i,0,list_doublet[[i]][2])
    }
}

list_doublet <- data.frame(list_doublet)
list_doublet <- data.frame(t(list_doublet))

list_doublet$percentage <- round(as.numeric(list_doublet$Doublet)/(as.numeric(list_doublet$Doublet)+as.numeric(list_doublet$Singlet))*100,2)
colnames(list_doublet)[1] <- "cluster"
list_doublet$cluster <- factor(x = list_doublet$cluster,levels = c(0:(nrow(list_doublet)-1)))

pdf("doublet_percent.pdf",width = 8,height = 6)
ggplot(data = list_doublet,mapping = aes(x = cluster,y = percentage))+
    geom_col(mapping = aes(fill = cluster))+
    geom_text(mapping = aes(y = percentage+1,label = percentage))+
    ggtitle(label = "Percentage of doublet cells in each cluster")+
    theme_bw()
dev.off()

pdf(file = paste0(des_dir,"dimplot_cluster_res.1_sample_split.pdf"),width =20,height = 18 )
DimPlot(seurat_obj,label = T,reduction = "umap",split.by = "sample",raster = F,ncol = 2)
dev.off()
