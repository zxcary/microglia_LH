setwd(dir = "./")
library(Seurat)
library(harmony)
data_seurat <- readRDS(file = "./data/seurat_anno_final.rds")
DimPlot(object = data_seurat,reduction = "umap",label = T)

microgila <- subset(x = data_seurat,idents = "Microgila")


#########################################

DefaultAssay(microgila) <- "RNA"
microgila = NormalizeData(microgila)
microgila = FindVariableFeatures(microgila,nfeatures = 2000)%>% ScaleData()


microgila = RunPCA(microgila,verbose = F)
microgila <- RunHarmony(microgila, group.by.vars = "sample",plot_convergence = T)

microgila <-  microgila%>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = c(0.4,0.5,0.6,0.8,1,1.2,1.4,1.5))
microgila <- RunTSNE(microgila,reduction = "harmony",dims = 1:30)
microgila <- RunUMAP(microgila,reduction = "harmony",dims = 1:30)
############
DefaultAssay(object = microgila) <- "RNA"
genes_rp <- grep(pattern = "^Rpl|^Rps",x = row.names(x = microgila@assays$RNA),value = T)
#microgila = ScaleData(microgila,vars.to.regress=genes_rp)

# microgila <- subset(x = microgila,idents = "1",invert = T)

Idents(object = microgila) <- "RNA_snn_res.0.5"
DimPlot(object = microgila,reduction = "umap",label = T)
DimPlot(object = microgila,reduction = "tsne",label = T)

microgila <- RenameIdents(object = microgila,"0" = "MG1","1" = "MG2","2" = "MG3","3" = "MG4","4" = "MG5","5" = "MG6","6" = "MG7","7" = "MG8","8" = "MG9")
microgila$anno_detail <- factor(x = microgila@active.ident,levels = paste0("MG",1:9))
Idents(microgila) <- "anno_detail"
color_model=c("#ffdd00","#ae63e4","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
pdf(file = "dimplot.pdf",width = 8,height = 7)
DimPlot(object = microgila,reduction = "umap",label = T,cols = color_model)
dev.off()
pdf(file = "dimplot_tsne.pdf",width = 8,height = 7)
DimPlot(object = microgila,reduction = "tsne",label = T,cols = color_model)
dev.off()
###################
Idents(object = microgila) <- "RNA_snn_res.0.5"
microgila <- RenameIdents(object = microgila,"0" = "Alternative homeostatic","1" = "Others","2" = "Alternative homeostatic","3" = "Alternative homeostatic","4" = "Classical homeostatic","5" = "DAM-like","6" = "DAM-like","7" = "DAM-like","8" = "Classical homeostatic")
microgila$anno_rename <- factor(x = microgila@active.ident,levels = c("Classical homeostatic","Alternative homeostatic","DAM-like","Others"))
#Idents(microgila) <- "anno_detail"
# sss <- FindMarkers(object = microgila,ident.1 = c("MG6","MG7","MG8"),logfc.threshold = 0.5)
# sss <- sss[order(sss$avg_log2FC,decreasing = T),]
# degs_use <- AverageExpression(object = microgila,assays = "RNA",features = rownames(sss))$RNA
# 
# pheatmap(mat = degs_use,scale = "row",cluster_cols = T,cluster_rows = T,main = "degs_use",filename = "pheatmap_degs_use.pdf",width = 8,height = 16)

pdf(file = paste0("dimplot_tsne.pdf"),width = 8,height = 6)
DimPlot(object = microgila,reduction = "tsne",label = T)
dev.off()

Idents(microgila) <- "anno_rename"
pdf(file = paste0("dimplot_rename_tsne.pdf"),width = 8.5,height = 6)
DimPlot(object = microgila,reduction = "tsne",cols = c("#ff0000","#00a4e4","#c1d82f","#6a67ce"))
dev.off()
pdf(file = paste0("dimplot_rename_umap.pdf"),width = 8.5,height = 6)
DimPlot(object = microgila,reduction = "umap",cols = c("#ff0000","#00a4e4","#c1d82f","#6a67ce"))
dev.off()
microgila_clean <- subset(x = microgila,idents = "Others",invert = T)
# DefaultAssay(object = microgila_clean) <- "RNA"
Idents(microgila_clean) <- "anno_rename"
degs_3 <- FindAllMarkers(object = microgila_clean,assay = "RNA",logfc.threshold = 0.5,only.pos = T)
write.csv(x = degs_3,file = "degs_3_0.5.csv",row.names = F)
DefaultAssay(object = microgila) <- "RNA"
FeaturePlot(object = microgila,features = "mt-Co2")
Idents(microgila) <- "anno_detail"
gene_hd <- read.xlsx(xlsxFile = "PAM DAM homeostatic genes.xlsx",sheet = 1)
c_genes <- intersect(x = gene_hd$Homeostatic,degs_3[degs_3$cluster == "Classical homeostatic","gene"])
a_genes <- intersect(x = gene_hd$Homeostatic,degs_3[degs_3$cluster == "Alternative homeostatic","gene"])
d_genes <- intersect(x = gene_hd$DAM,degs_3[degs_3$cluster == "DAM-like","gene"])

classical_degs <- AverageExpression(object = microgila,assays = "RNA",features = c_genes)$RNA[,-2]
alternative_degs <- AverageExpression(object = microgila,assays = "RNA",features = a_genes)$RNA[,-2]
dam_degs <- AverageExpression(object = microgila,assays = "RNA",features = d_genes)$RNA[,-2]

classical_degs <- AverageExpression(object = microgila,assays = "RNA",features = degs_3[degs_3$cluster == "Classical homeostatic","gene"])$RNA[1:200,]
alternative_degs <- AverageExpression(object = microgila,assays = "RNA",features = degs_3[degs_3$cluster == "Alternative homeostatic","gene"])$RNA
dam_degs <- AverageExpression(object = microgila,assays = "RNA",features = degs_3[degs_3$cluster == "DAM-like","gene"])$RNA[1:200,]
library(openxlsx)

library(pheatmap)
pheatmap(mat = classical_degs,scale = "row",cluster_cols = T,cluster_rows = T,main = "Classical homeostatic",filename = "pheatmap_Classical_homeostatic.pdf",width = 8,height = 6)
pheatmap(mat = alternative_degs,scale = "row",cluster_cols = T,cluster_rows = T,main = "Alternative homeostatic",filename = "pheatmap_Alternative_homeostatic.pdf",width = 8,height = 6)
pheatmap(mat = dam_degs,scale = "row",cluster_cols = T,cluster_rows = T,main = "DAM-like",filename = "pheatmap_DAM-like.pdf",width = 8,height = 8)
# microgila <- microgila %>% RunUMAP(reduction = "harmony", dims = 1:20)
# microgila <- microgila %>% RunTSNE(reduction = "harmony", dims = 1:20)
Idents(microgila) <- "RNA_snn_res.0.5"
Idents(microgila) <- "doubFind_res"
pdf(file = paste0(des_dir,"microgila_dimplot_cluster_res.1.pdf"),width =10,height = 8 )
DimPlot(microgila,label = T,reduction = "tsne",raster = F)
dev.off()
FeaturePlot(object = microgila,features = c(c("Ccl4","Ifit2","Malat1","Jun","Rps7","P2ry12")),reduction = "tsne",raster = F,label = F,ncol = 3)
library(MySeuratWrappers)
VlnPlot(object = microgila,features = c(c("Ccl4","Ifit2","Malat1","Jun","Rps7","P2ry12")),stacked = T,pt.size = 0)

###########################3 MG2

list_degs <- list("Classical homeostatic" = c("Ube2d3","Gm47283","Hnrnpdl","Lacc1","Cflar","Tet2","Zc3h7a","Dhx9","Fnbp4","Frmd4b","Ankrd44","Atrx","Ssh2","Maml3","Srsf11","Rbm39","Ash1l","Dleu2","Stab1","Atp8a1","Pkn1","Mef2c","Fgd2","Arhgap17","Setd5","Malat1","Ptbp2","Eif4a2","R3hdm1","Rhoh","Gm43462"),
                  "Alternative homeostatic"=c("Fcrls","Olfml3","Marcks","Cx3cr1","P2ry12","Gpr34","Slc2a5","Selplg","P2ry13"),
                  "DAM−like" = c("Saa3","Marcksl1","Apoe","Cxcl16","Cd74","Bcl2a1b","Spp1","Bst2","Ifitm3","Ms4a6c","Tspo","Ms4a6d","Srgn","Mt1","Mt2","Sod2","Ccl5","AW112010")
)
ggg <- unlist(list_degs)
Idents(microgila) <- "anno_detail"
sep_degs <- AverageExpression(object = microgila,assays = "RNA",features = ggg)$RNA
sep_degs <- sep_degs[,c("MG5","MG9","MG1","MG3","MG4","MG7","MG6","MG8")]
# sep_degs <- sep_degs[,]
pheatmap(mat = sep_degs,scale = "row",cluster_cols = T,cluster_rows = T,clustering_method = "complete",main = "genes",filename = "pheatmap_genes.pdf",width = 8,height = 16)

write.csv(x = sep_degs,file = "pheatmap_genes.csv")


################33 degs three condition 0.5
DefaultAssay(object =microgila ) <- "RNA"
Idents(microgila) <- "re_names"
degs1 <- FindMarkers(object = microgila,ident.1 = "LPS",ident.2 = "Con",logfc.threshold = 0.5)
degs2 <- FindMarkers(object = microgila,ident.1 = "HYP",ident.2 = "Con",logfc.threshold = 0.5)
degs3 <- FindMarkers(object = microgila,ident.1 = "L/H",ident.2 = "Con",logfc.threshold = 0.5)
write.xlsx(x = list("LPS"=degs1,"HYP"=degs2,"L/H"=degs3),file = "degs_groups_0.5.xlsx",rowNames = T)
degs_number <- data.frame("Group"=rep(c("LPS","HYP","L/H"),c(2,2,2)),"Regulated"=rep(c("Up","Down"),3))
degs_number$Number <- c(nrow(degs1[degs1$avg_log2FC>0,]),nrow(degs1[degs1$avg_log2FC<0,]),
                        nrow(degs2[degs2$avg_log2FC>0,]),nrow(degs2[degs2$avg_log2FC<0,]),
                        nrow(degs3[degs3$avg_log2FC>0,]),nrow(degs3[degs3$avg_log2FC<0,]))
degs_number$Group <- factor(x = degs_number$Group,levels = c("LPS","HYP","L/H"))
degs_number$Regulated <- factor(x = degs_number$Regulated,levels = c("Up","Down"))
degs_number$Number2 <- c(degs_number$Number[1]/2+degs_number$Number[2],degs_number$Number[2]/2,degs_number$Number[3]/2+degs_number$Number[4],degs_number$Number[4]/2,
                         degs_number$Number[5]/2+degs_number$Number[6],degs_number$Number[6]/2)
library(ggplot2)
pdf(file = "DEGs_number_0.5.pdf",width = 5,height = 6)
ggplot(data = degs_number,mapping = aes(x = Group,y = Number,group = Regulated))+
    geom_col(mapping = aes(fill = Regulated))+
    xlab(label = "")+
    ylab(label = "Number of DEGs")+
    ggtitle(label = "DEGs of three groups compared with Con")+
    geom_text(mapping = aes(y = Number2,label = Number))+
    theme_bw()
dev.off()

genes_inflam <- read.xlsx(xlsxFile = "小鼠炎症因子.xlsx",sheet = 1)
degs_number <- data.frame("Group"=rep(c("LPS","HYP","L/H"),c(2,2,2)),"Regulated"=rep(c("Up","Down"),3))
degs_number$Number <- c(nrow(degs1[degs1$avg_log2FC>0 & (rownames(degs1) %in% genes_inflam$Symbol),]),nrow(degs1[degs1$avg_log2FC<0 & (rownames(degs1) %in% genes_inflam$Symbol),]),
                        nrow(degs2[degs2$avg_log2FC>0 & (rownames(degs2) %in% genes_inflam$Symbol),]),nrow(degs2[degs2$avg_log2FC<0 & (rownames(degs2) %in% genes_inflam$Symbol),]),
                        nrow(degs3[degs3$avg_log2FC>0 & (rownames(degs3) %in% genes_inflam$Symbol),]),nrow(degs3[degs3$avg_log2FC<0 & (rownames(degs3) %in% genes_inflam$Symbol),]))
degs_number$Group <- factor(x = degs_number$Group,levels = c("LPS","HYP","L/H"))
degs_number$Regulated <- factor(x = degs_number$Regulated,levels = c("Up","Down"))
degs_number$Number2 <- c(degs_number$Number[1]/2+degs_number$Number[2],degs_number$Number[2]/2,degs_number$Number[3]/2+degs_number$Number[4],degs_number$Number[4]/2,
                         degs_number$Number[5]/2+degs_number$Number[6],degs_number$Number[6]/2)
library(ggplot2)
pdf(file = "Inflams_number_0.5.pdf",width = 5,height = 6)
ggplot(data = degs_number,mapping = aes(x = Group,y = Number,group = Regulated))+
    geom_col(mapping = aes(fill = Regulated))+
    xlab(label = "")+
    ylab(label = "Number of inflammatory factors")+
    ggtitle(label = "Inflammatory factors of three DEGs")+
    geom_text(mapping = aes(y = Number2,label = Number))+
    theme_bw()
dev.off()



################33 degs three condition 0.5
DefaultAssay(object =microgila ) <- "RNA"
Idents(microgila) <- "re_names"
degs1 <- FindMarkers(object = microgila,ident.1 = "LPS",ident.2 = "Con",logfc.threshold = 0.25)
degs2 <- FindMarkers(object = microgila,ident.1 = "HYP",ident.2 = "Con",logfc.threshold = 0.25)
degs3 <- FindMarkers(object = microgila,ident.1 = "L/H",ident.2 = "Con",logfc.threshold = 0.25)
write.xlsx(x = list("LPS"=degs1,"HYP"=degs2,"L/H"=degs3),file = "degs_groups_0.25.xlsx",rowNames = T)
degs_number <- data.frame("Group"=rep(c("LPS","HYP","L/H"),c(2,2,2)),"Regulated"=rep(c("Up","Down"),3))
degs_number$Number <- c(nrow(degs1[degs1$avg_log2FC>0,]),nrow(degs1[degs1$avg_log2FC<0,]),
                        nrow(degs2[degs2$avg_log2FC>0,]),nrow(degs2[degs2$avg_log2FC<0,]),
                        nrow(degs3[degs3$avg_log2FC>0,]),nrow(degs3[degs3$avg_log2FC<0,]))
degs_number$Group <- factor(x = degs_number$Group,levels = c("LPS","HYP","L/H"))
degs_number$Regulated <- factor(x = degs_number$Regulated,levels = c("Up","Down"))
degs_number$Number2 <- c(degs_number$Number[1]/2+degs_number$Number[2],degs_number$Number[2]/2,degs_number$Number[3]/2+degs_number$Number[4],degs_number$Number[4]/2,
                         degs_number$Number[5]/2+degs_number$Number[6],degs_number$Number[6]/2)
library(ggplot2)
pdf(file = "DEGs_number_0.25.pdf",width = 5,height = 6)
ggplot(data = degs_number,mapping = aes(x = Group,y = Number,group = Regulated))+
    geom_col(mapping = aes(fill = Regulated))+
    xlab(label = "")+
    ylab(label = "Number of DEGs")+
    ggtitle(label = "DEGs of three groups compared with Con")+
    geom_text(mapping = aes(y = Number2,label = Number))+
    theme_bw()
dev.off()

genes_inflam <- read.xlsx(xlsxFile = "小鼠炎症因子.xlsx",sheet = 1)
degs_number <- data.frame("Group"=rep(c("LPS","HYP","L/H"),c(2,2,2)),"Regulated"=rep(c("Up","Down"),3))
degs_number$Number <- c(nrow(degs1[degs1$avg_log2FC>0 & (rownames(degs1) %in% genes_inflam$Symbol),]),nrow(degs1[degs1$avg_log2FC<0 & (rownames(degs1) %in% genes_inflam$Symbol),]),
                        nrow(degs2[degs2$avg_log2FC>0 & (rownames(degs2) %in% genes_inflam$Symbol),]),nrow(degs2[degs2$avg_log2FC<0 & (rownames(degs2) %in% genes_inflam$Symbol),]),
                        nrow(degs3[degs3$avg_log2FC>0 & (rownames(degs3) %in% genes_inflam$Symbol),]),nrow(degs3[degs3$avg_log2FC<0 & (rownames(degs3) %in% genes_inflam$Symbol),]))
degs_number$Group <- factor(x = degs_number$Group,levels = c("LPS","HYP","L/H"))
degs_number$Regulated <- factor(x = degs_number$Regulated,levels = c("Up","Down"))
degs_number$Number2 <- c(degs_number$Number[1]/2+degs_number$Number[2],degs_number$Number[2]/2,degs_number$Number[3]/2+degs_number$Number[4],degs_number$Number[4]/2,
                         degs_number$Number[5]/2+degs_number$Number[6],degs_number$Number[6]/2)
library(ggplot2)
pdf(file = "Inflams_number_0.25.pdf",width = 5,height = 6)
ggplot(data = degs_number,mapping = aes(x = Group,y = Number,group = Regulated))+
    geom_col(mapping = aes(fill = Regulated))+
    xlab(label = "")+
    ylab(label = "Number of inflammatory factors")+
    ggtitle(label = "Inflammatory factors of three DEGs")+
    geom_text(mapping = aes(y = Number2,label = Number))+
    theme_bw()
dev.off()
##################

inflam_degs <- unique(c(rownames(degs1[rownames(degs1) %in% genes_inflam$Symbol,]),rownames(degs2[rownames(degs2) %in% genes_inflam$Symbol,]),rownames(degs3[rownames(degs3) %in% genes_inflam$Symbol,])))
inflam_degs_exp <- data.frame(AverageExpression(object = microgila,assays = "RNA",features = inflam_degs)$RNA)
# inflam_degs_exp <- data.frame(AverageExpression(object = microgila_clean,assays = "RNA",features = c("Ccl5"))$RNA)
inflam_degs_exp2 <- data.frame("Gene" = rep(rownames(inflam_degs_exp),4),"Group"=rep(c("Con","LPS","HYP","L/H"),rep(nrow(inflam_degs_exp),4)),"exp" = unlist(inflam_degs_exp))
write.csv(x = inflam_degs_exp,file = "inflam_degs_exp.csv")
inflam_degs_exp2$Group <- factor(x = inflam_degs_exp2$Group,levels = c("Con","LPS","HYP","L/H"))
my_comparisons <- list(c("Moderate", "Severe"), c("Conv", "Severe"))
color_vec <- c("Con"="#1B9E77","LPS"="#D95F02","HYP"="#FDBF6F","L/H"="#7570B3")
pdf(file = "Inflams_gene_expression.pdf",width = 20,height = 8)
ggplot(data = inflam_degs_exp2,mapping = aes(x = Group, y = exp, group = Group))+
    geom_col(mapping = aes(x = Group, y = exp, fill = Group))+
    geom_text(mapping = aes( y = exp+0.1, color = Group,label= round(exp,2)))+
    # geom_point(mapping = aes(x = group, y = percentage2, fill = group),size=0.1) +
    # stat_summary(mapping = aes(x = group, y = percentage2, group = group),
    #              fun.data="mean_cl_normal", fun.args = list(mult=1),
    #              geom = "errorbar", color = "black", width=0.01) +
    # stat_compare_means(method = "wilcox.test",comparisons=my_comparisons, 
    #                    hide.ns = TRUE,vjust = 0.5,
    #                    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
    theme_classic() +
    ylab(label = "Average Expression")+
    facet_grid(~Gene,scales = "free")+
    theme(plot.title = element_text(hjust = 0.5),
          strip.text.x = element_text(size = 15),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_fill_manual(values = color_vec)+
    scale_color_manual(values = color_vec)
dev.off()

############################## 4 gene heatmap
microgila_clean <- readRDS(file = "./data/microgila.rds")
data_ccl5 = AverageExpression(object = microgila_clean,assays = "RNA",features = c(paste0("Ccl",1:10)),group.by = c("anno_detail","re_names"))$RNA
degs1 <- read.xlsx(xlsxFile = "result/degs/degs_groups_0.25.xlsx",sheet = "LPS",rowNames = T)
degs2 <- read.xlsx(xlsxFile = "result/degs/degs_groups_0.25.xlsx",sheet = "HYP",rowNames = T)
degs3 <- read.xlsx(xlsxFile = "result/degs/degs_groups_0.25.xlsx",sheet = "L/H",rowNames = T)

# genes_use <- unique(c(rownames(degs1[order(abs(degs1$avg_log2FC),decreasing = T),][1:50,]),
#                       rownames(degs2[order(abs(degs2$avg_log2FC),decreasing = T),][1:50,]),
#                       rownames(degs3[order(abs(degs3$avg_log2FC),decreasing = T),][1:50,])))
genes_use <- unique(c(rownames(degs1[order(degs1$avg_log2FC,decreasing = T),][1:50,]),
                      rownames(degs2[order(degs2$avg_log2FC,decreasing = T),][1:50,]),
                      rownames(degs3[order(degs3$avg_log2FC,decreasing = T),][1:50,])))
microgila_clean <- ScaleData(object = microgila_clean,features = genes_use)
png(filename = "pheatmap_4_conditions.png",width = 3000,height = 4000,res = 300)
Seurat::DoHeatmap(object = microgila_clean,features =  genes_use,label = T,draw.lines = T,raster = F,assay = "RNA")
dev.off()
pdf(file = "pheatmap_4_conditions.pdf",width = 12,height = 16)
Seurat::DoHeatmap(object = microgila_clean,features =  genes_use,label = T,draw.lines = T,raster = F,assay = "RNA")
dev.off()
DimPlot(object = microgila_clean,reduction = "umap")

###############################

microgila_clean <- readRDS(file = "./data/microgila.rds")

genes_use <- list("LPS"=rownames(degs1),
                  "IF" = read.xlsx(xlsxFile = "小鼠炎症因子.xlsx",sheet = 1)[,1],
                  "HYP"=rownames(degs2),
                  "L/H"=rownames(degs3))
#genes_use[["IF"]] <- read.xlsx(xlsxFile = "小鼠 炎症因子 孝昌.xlsx",sheet = 1)[,1]

library(VennDiagram)
color_model <- c("#FF8080","#C5944E","#099963","#87CEFA")
# color_model <- c("#F6AAA7","#ADD2A9","#A5A7FB")
list_sig_genes2 <- genes_use
venn.plot <- venn.diagram(
    x = list_sig_genes2,imagetype = "svg",
    filename = paste0("venn_3_degs_inflame.svg"),
    # col = "transparent",
    fill = color_model,
    alpha = 0.5,
    # main = "Down regulated",main.cex = 3,main.fontfamily = "arial", 
    # label.col = c("darkred", "white", "darkblue", "white",
    #               "white", "white", "darkgreen"),
    cex = 2,#内标签的字体大小
    fontfamily = "arial",
    # fontface = "bold",
    cat.default.pos = "outer",#设置标签在圆外面
    #cat.col = c("darkred", "darkblue", "darkgreen"),
    cat.cex = 2.5,#外标签的字体大小
    cat.fontfamily = "arial",
    # cat.dist = c(0.05, 0.05, 0.05,0.05),#相对圆圈的位置
    # cat.pos = c(-20,20,180,180)  #相对12点方向旋转的角度
)



venn.plot <- venn.diagram(
    x = list_sig_genes2,
    filename = NULL,
    # col = "transparent",
    fill = color_model,
    alpha = 0.5,
    # main = "Down regulated",main.cex = 3,main.fontfamily = "arial", 
    # label.col = c("darkred", "white", "darkblue", "white",
    #               "white", "white", "darkgreen"),
    cex = 2,#内标签的字体大小
    #fontfamily = "arial",
    # fontface = "bold",
    cat.default.pos = "outer",#设置标签在圆外面
    #cat.col = c("darkred", "darkblue", "darkgreen"),
    cat.cex = 2.5,#外标签的字体大小
    #cat.fontfamily = "arial",
    # cat.dist = c(0.05, 0.05, 0.05,0.05),#相对圆圈的位置
    # cat.pos = c(-20,20,180,180)  #相对12点方向旋转的角度
)

pdf("venn_3_degs_inflame.pdf")
grid.draw(venn.plot)
dev.off()

inter <- get.venn.partitions(list_sig_genes2)

for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(6)], 'venn_inter_4.txt', row.names = FALSE, sep = '\t', quote = FALSE)


microgila_clean <- ScaleData(object = microgila_clean,features = genes_use)
png(filename = "pheatmap_4_conditions.png",width = 3000,height = 4000,res = 300)
Seurat::DoHeatmap(object = microgila_clean,features =  genes_use,label = T,draw.lines = T,raster = F,assay = "RNA")
dev.off()
pdf(file = "pheatmap_4_conditions.pdf",width = 12,height = 16)
Seurat::DoHeatmap(object = microgila_clean,features =  genes_use,label = T,draw.lines = T,raster = F,assay = "RNA")
dev.off()
DimPlot(object = microgila_clean,reduction = "umap")
###################################
color_model <- c("#ffdd00","#ae63e4","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
Idents(microgila_clean) <- "anno_detail"
pdf(file = "Vlnplot_markers.pdf",width = 8,height = 8)
MySeuratWrappers::VlnPlot(object = microgila_clean,features = c(c("P2ry12","Gpr34","Cd52","Tpt1","Jun","Ddx17","Eef1b2","Ccl4","Fth1","Ifit2","Malat1","Rps7")),
                          stacked = T,pt.size = 0,cols = color_model)
dev.off()
##################################
genes_use2 <- intersect(Reduce(f  = union,x = genes_use[c("LPS","HYP","L/H")]),genes_use$IF)

microgila_clean$subtype_condition <- paste0(microgila_clean$group,"_",microgila_clean$anno_detail)
microgila_clean$re_names <- factor(x = microgila_clean$re_names,levels = c("Con","LPS","HYP","L/H"))
data_plot <- AverageExpression(object = microgila_clean,assays = "RNA",features = genes_use2,group.by = c("anno_detail","re_names"))$RNA


group_names1 <- rep(c("Con","LPS","HYP","L/H"),9)
group_names2 <- rep(paste0("MG",1:9),rep(4,9))
anno_group <- data.frame("Group"=group_names1,"Subtypes"=group_names2)
anno_group$Group <- factor(anno_group$Group,c("Con","LPS","HYP","L/H"))
rownames(x = anno_group) <- colnames(data_plot)
color_vec <- c("Con"="#1B9E77","LPS"="#D95F02","HYP"="#FDBF6F","L/H"="#7570B3")
color_model <- c("#ffdd00","#ae63e4","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
names(x = color_model) <- paste0("MG",1:9)

anno_color2 = list("Group" = color_vec,"Subtypes" = color_model)
data_plot <- data_ccl5
bk <- c(seq(-4,4,by=0.01))
pheatmap(mat = data_plot,display_numbers = round(data_plot,2),
         breaks = bk,color = c(colorRampPalette(c("#E6E6FA", "white"))(400), colorRampPalette(c("white","firebrick3"))(400)),
         cluster_cols = F,cluster_rows = T,scale = "row",
         #filename = "heatmap_IF_group_celltype_xlab.pdf",width = 12,height = 6,
         annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Average Expression",gaps_col = c(4,8,12,16,20,24,28,32))

pheatmap(mat = data_plot,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("#E6E6FA", "white"))(400), colorRampPalette(c("white","firebrick3"))(400)),
         cluster_cols = F,cluster_rows = T,scale = "row",labels_col = rep("",36),
         filename = "heatmap_IF_group_celltype_noxlab.pdf",width = 12,height = 5,
         annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Average Expression",gaps_col = c(4,8,12,16,20,24,28,32))

#################################3 no_MG2
##################################
genes_use2 <- intersect(Reduce(f  = union,x = genes_use[c("LPS","HYP","L/H")]),genes_use$IF)

microgila_clean$subtype_condition <- paste0(microgila_clean$group,"_",microgila_clean$anno_detail)
microgila_clean$re_names <- factor(x = microgila_clean$re_names,levels = c("Con","LPS","HYP","L/H"))
data_plot <- AverageExpression(object = microgila_clean,assays = "RNA",features = genes_use2,group.by = c("anno_detail","re_names"))$RNA
write.csv(x = data_plot,file = "heatmap_IF_group_celltype.csv")

data_plot <- data_plot[,-c(5,6,7,8)]

group_names1 <- rep(c("Con","LPS","HYP","L/H"),8)
group_names2 <- rep(paste0("MG",c(1,3:9)),rep(4,8))
anno_group <- data.frame("Group"=group_names1,"Subtypes"=group_names2)
anno_group$Group <- factor(anno_group$Group,c("Con","LPS","HYP","L/H"))
rownames(x = anno_group) <- colnames(data_plot)
color_vec <- c("Con"="#1B9E77","LPS"="#D95F02","HYP"="#FDBF6F","L/H"="#7570B3")
color_model <- c("#ffdd00","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
names(x = color_model) <- paste0("MG",c(1,3:9))

anno_color2 = list("Group" = color_vec,"Subtypes" = color_model)

bk <- c(seq(-4,4,by=0.01))
pheatmap(mat = data_plot,display_numbers = round(data_plot,2),
         breaks = bk,color = c(colorRampPalette(c("#E6E6FA", "white"))(400), colorRampPalette(c("white","firebrick3"))(400)),
         cluster_cols = F,cluster_rows = T,scale = "row",
         filename = "heatmap_IF_group_celltype_xlab_no_MG2_display.pdf",width = 12,height = 6,
         annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Average Expression",gaps_col = c(4,8,12,16,20,24,28,32))

pheatmap(mat = data_plot,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("#E6E6FA", "white"))(400), colorRampPalette(c("white","firebrick3"))(400)),
         cluster_cols = F,cluster_rows = T,scale = "row",labels_col = rep("",36),
         filename = "heatmap_IF_group_celltype_noxlab_no_MG2.pdf",width = 12,height = 5,
         annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Average Expression",gaps_col = c(4,8,12,16,20,24,28,32))



#################################
genes_inf <- read.csv(file = "inflame.csv")
genes_inf <- genes_inf[!duplicated(genes_inf$Gene),]
mouse2human <- read.csv(file = "mouse2huamn_genes.csv")
genes_inf <- merge(x = genes_inf,y = mouse2human,by.x = 1,by.y = 2,all.x = T)
genes_inf <- genes_inf[!duplicated(genes_inf$MGI.symbol),]
genes_inf <- genes_inf[!is.na(genes_inf$MGI.symbol),]
genes_use2 <- genes_inf$MGI.symbol
write.csv(x = genes_inf,file = "inflame_genes.csv")

genes_use3 <- intersect(Reduce(f  = union,x = genes_use[c("LPS","HYP","L/H")]),genes_use2)

microgila_clean$subtype_condition <- paste0(microgila_clean$group,"_",microgila_clean$anno_detail)
microgila_clean$re_names <- factor(x = microgila_clean$re_names,levels = c("Con","LPS","HYP","L/H"))
data_plot <- AverageExpression(object = microgila_clean,assays = "RNA",features = genes_use3,group.by = c("anno_detail","re_names"))$RNA


group_names1 <- rep(c("Con","LPS","HYP","L/H"),9)
group_names2 <- rep(paste0("MG",1:9),rep(4,9))
anno_group <- data.frame("Group"=group_names1,"Subtypes"=group_names2)
anno_group$Group <- factor(anno_group$Group,c("Con","LPS","HYP","L/H"))
rownames(x = anno_group) <- colnames(data_plot)
color_vec <- c("Con"="#1B9E77","LPS"="#D95F02","HYP"="#FDBF6F","L/H"="#7570B3")
color_model <- c("#ffdd00","#ae63e4","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
names(x = color_model) <- paste0("MG",1:9)

anno_color2 = list("Group" = color_vec,"Subtypes" = color_model)

bk <- c(seq(-4,4,by=0.01))
pheatmap(mat = data_plot,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("#E6E6FA", "white"))(400), colorRampPalette(c("white","firebrick3"))(400)),
         cluster_cols = F,cluster_rows = T,scale = "row",
         filename = "heatmap_IF_group_celltype_xlab.pdf",width = 12,height = 18,
         annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Average Expression",gaps_col = c(4,8,12,16,20,24,28,32))

pheatmap(mat = data_plot,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("#E6E6FA", "white"))(400), colorRampPalette(c("white","firebrick3"))(400)),
         cluster_cols = F,cluster_rows = T,scale = "row",labels_col = rep("",36),
         filename = "heatmap_IF_group_celltype_noxlab.pdf",width = 12,height = 18,
         annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Average Expression",gaps_col = c(4,8,12,16,20,24,28,32))

############################### no MG2
#################################
genes_inf <- read.csv(file = "inflame.csv")
genes_inf <- genes_inf[!duplicated(genes_inf$Gene),]
mouse2human <- read.csv(file = "mouse2huamn_genes.csv")
genes_inf <- merge(x = genes_inf,y = mouse2human,by.x = 1,by.y = 2,all.x = T)
genes_inf <- genes_inf[!duplicated(genes_inf$MGI.symbol),]
genes_inf <- genes_inf[!is.na(genes_inf$MGI.symbol),]
genes_use2 <- genes_inf$MGI.symbol
write.csv(x = genes_inf,file = "inflame_genes.csv")

genes_use3 <- intersect(Reduce(f  = union,x = genes_use[c("LPS","HYP","L/H")]),genes_use2)

microgila_clean$subtype_condition <- paste0(microgila_clean$group,"_",microgila_clean$anno_detail)
microgila_clean$re_names <- factor(x = microgila_clean$re_names,levels = c("Con","LPS","HYP","L/H"))
data_plot <- AverageExpression(object = microgila_clean,assays = "RNA",features = genes_use3,group.by = c("anno_detail","re_names"))$RNA
write.csv(x = data_plot,file = "heatmap_IF2_group_celltype.csv")
data_plot <- data_plot[,-c(5,6,7,8)]


group_names1 <- rep(c("Con","LPS","HYP","L/H"),8)
group_names2 <- rep(paste0("MG",c(1,3:9)),rep(4,8))
anno_group <- data.frame("Group"=group_names1,"Subtypes"=group_names2)
anno_group$Group <- factor(anno_group$Group,c("Con","LPS","HYP","L/H"))
rownames(x = anno_group) <- colnames(data_plot)
color_vec <- c("Con"="#1B9E77","LPS"="#D95F02","HYP"="#FDBF6F","L/H"="#7570B3")
color_model <- c("#ffdd00","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
names(x = color_model) <- paste0("MG",c(1,3:9))

anno_color2 = list("Group" = color_vec,"Subtypes" = color_model)

bk <- c(seq(-4,4,by=0.01))
pheatmap(mat = data_plot,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("#E6E6FA", "white"))(400), colorRampPalette(c("white","firebrick3"))(400)),
         cluster_cols = F,cluster_rows = T,scale = "row",
         filename = "heatmap_IF2_group_celltype_xlab_no_MG2.pdf",width = 12,height = 18,
         annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Average Expression",gaps_col = c(4,8,12,16,20,24,28,32))

pheatmap(mat = data_plot,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("#E6E6FA", "white"))(400), colorRampPalette(c("white","firebrick3"))(400)),
         cluster_cols = F,cluster_rows = T,scale = "row",labels_col = rep("",36),
         filename = "heatmap_IF2_group_celltype_noxlab_no_MG2.pdf",width = 12,height = 18,
         annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Average Expression",gaps_col = c(4,8,12,16,20,24,28,32))

#############################
library(ggrepel)
degs1 <- FindMarkers(object = microgila,ident.1 = "LPS",ident.2 = "Con",logfc.threshold = 0)
degs2 <- FindMarkers(object = microgila,ident.1 = "HYP",ident.2 = "Con",logfc.threshold = 0)
degs3 <- FindMarkers(object = microgila,ident.1 = "L/H",ident.2 = "Con",logfc.threshold = 0)
write.xlsx(x = list("LPS"=degs1,"HYP"=degs2,"L/H"=degs3),file = "degs_groups_0.xlsx",rowNames = T)
log2fc <- 0.25
group1 <- "HYP"
data_plot <- degs2
data_plot$Gene  <- rownames(data_plot)
data_plot$Group <- "notsignificant"
data_plot$Group[data_plot$avg_log2FC>=log2fc & data_plot$p_val_adj < 0.05] <- "upregulated"
data_plot$Group[data_plot$avg_log2FC<= -log2fc & data_plot$p_val_adj < 0.05] <- "downregulated"
data_plot$Group <- factor(x = data_plot$Group,levels = c("upregulated","downregulated","notsignificant"))
number0 <- sum(data_plot$p_val_adj < 10**(-300)) 
data_plot$p_val_adj[data_plot$p_val_adj < 10**(-300)] <- 10** (-300-sample(x = seq(-10,10,0.1),size = number0,replace = T))
write.csv(x = data_plot,file = paste0("degs_",group1,"_vs_","Con",".csv"))
#write.csv(x = data_plot,file = paste0("degs_","L.H","_vs_","Con",".csv"))
pdf(file = paste0("volcano_",group1,"_vs_","Con",".pdf"),width = 8,height = 6)
p <- ggplot(data = data_plot,mapping = aes(x = avg_log2FC,y = -log10(p_val_adj),color=Group))+
    geom_point(alpha=0.4, size=2) +
    scale_color_manual(values=c("#D53C50","#1F5DA7", "grey")) +
    xlim(c(-1.5, 1.5)) +
    geom_vline(xintercept=c(-log2fc,log2fc),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
    labs(x="log2(fold change)",y="-log10 (FDR)") +
    theme_bw()+
    ggtitle(label = paste0(group1,"_vs_","Con"))+
    theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
        data = subset(data_plot, data_plot$Group != "notsignificant"),
        aes(label = Gene),
        size = 3,
        box.padding = unit(0.1, "lines"),max.overlaps=10,
        point.padding = unit(0, "lines"), segment.color = "black", show.legend = FALSE)
print(p)
dev.off()


VlnPlot(object = microgila,features = "Ccl5")


###################33 
saveRDS(object = microgila,file = "./data/microgila.rds")
#features_use <- FindVariableFeatures(object = microgila_clean,nfeatures = 2000,assay = "RNA")@assays$RNA@var.features
genes_class <- read.xlsx(xlsxFile = "PAM DAM homeostatic genes.xlsx",sheet = 1)
genes_class <- read.xlsx(xlsxFile = "PAM DAM homeostatic genes.xlsx",sheet = 1)
degs_1 <- read.csv("degs_3_0.25.csv")
features_use <- intersect(x = c(genes_class$DAM,genes_class$Homeostatic),y = degs_1$gene)
features_use <- features_use[!is.na(features_use)]
features_use <- intersect(x = features_use,y = rownames(microgila_clean@assays$RNA@counts))
data_use <- t(data.frame(microgila_clean@assays$RNA@counts[features_use,]))

meta_data2 <- data.frame(row.names = gsub(pattern = "-",replacement = ".",x = rownames(microgila_clean@meta.data)),
                         "Subtype"= factor(x = microgila_clean@meta.data$anno_detail,levels = paste0("MG",c(1,3:9))))
data_use <- merge(x = data_use,y = meta_data2,by = 0)
rownames(data_use) <- data_use[,1]
data_use <- data_use[,-1]
colnames(data_use) <- gsub(pattern = "-",replacement = ".",x = colnames(data_use))
colnames(data_use) <- gsub(pattern = "^[0-9]",replacement = "A",x = colnames(data_use))

library(randomForest)
train_cells <- c()
cells_use <- paste0("MG",c(1,3:9))
for(i in cells_use){
    cell_tmp <- rownames(data_use[data_use$Subtype == i,])
    samlle_tmp <- sample(cell_tmp, length(cell_tmp)*0.75,replace = F)
    train_cells <- c(train_cells,samlle_tmp)
}

test_cells <- setdiff(x = rownames(data_use),y = train_cells)
write.csv(train_cells,file = "train_cells.csv")
write.csv(test_cells,file = "test_cells.csv")
write.csv(meta_data2,file = "metadata.csv")



mic_train <- data_use[train_cells, ]
mic_test <- data_use[test_cells, ]
mic_train.forest <- randomForest(Subtype~., data = mic_train, importance = TRUE)
mic_predict <- predict(mic_train.forest, mic_test)
# mic_predict <- predict(mic_train.forest, mic_train)

rrr <- data.frame("True"=mic_test[,"Subtype"],"Pred" = mic_predict)
ddd <- data.frame(as.matrix(table(rrr)))
cells_use <- paste0("MG",c(1,3:9))
cell_num <- table(mic_test[,"Subtype"])
cell_num2 <- data.frame("cluster"=names(cell_num),"Number"=cell_num)
percent1 <- c()
for (i in 1:nrow(ddd)){
    percent1 <- c(percent1,ddd[i,"Freq"]/(cell_num[ddd[i,"True"]]+cell_num[ddd[i,"Pred"]]))
}
ddd$Percent <- percent1
write.csv(x = )

true_cell <- c()
for (i in rownames(mic_test)){
    if (mic_test[i,"Subtype"] == mic_predict[i]){
        true_cell <- c(true_cell,i)
    }
}

cell_num
list_use <- c()
for (i in 1:7){
    for (j in (i+1):8){
        node_i <- names(cell_num)[i]
        node_j <- names(cell_num)[j]
        list_use[[paste0(i,"_",j)]] <- c(node_i,node_j,ddd[ddd$True == node_i & ddd$Pred == node_j,"Percent"]+ddd[ddd$True == node_j & ddd$Pred == node_i,"Percent"])
    }
}
edges2 <- data.frame(t(data.frame(list_use)))

edges2 <- edges[!edges$True == edges$Pred,]

nodes <- data.frame(cell_num)
edges <- ddd[,-3]
write.csv(x = nodes,file = "./cell_number.csv",quote = F,row.names = F)
write.csv(x = edges2,file = "./edges.csv",quote = F,row.names = F)


library(tidygraph)
net.tidy <- tbl_graph(
    nodes = nodes, edges = edges, directed = F
)
library(ggraph)
ggraph(net.tidy, layout = "graphopt") + 
    #geom_node_range(mapping = aes(x = Number))+
    geom_node_point(mapping = aes(width = Number)) + # 点信息
    #geom_node_point(col = 'gold',size = 2) + # 点信息
    geom_edge_link(aes(width = Percent), alpha = 0.8) +  # 边信息
    scale_edge_width(range = c(0, 0.2)) + # 控制粗细
    geom_node_text(aes(label = Node), repel = TRUE) + # 增加节点的标签，reple避免节点重叠
    labs(edge_width = "phone.call") + # 图例标签
    theme_graph()

#########################
test_cells <- read.csv("./classification/class2/test_cells.csv",row.names = 1)
meta_data2 <- read.csv(file = "./classification/class2/metadata.csv",row.names = 1)
mic_test <- data.frame(row.names = test_cells$x,"Subtype"=meta_data2[test_cells$x,])
res_list <- readRDS(file = "./classification/class2/res_list.rds")
cells_use <- paste0("MG",c(1,3:9))
cell_num <- table(mic_test[,"Subtype"])
cell_num2 <- data.frame("cluster"=names(cell_num),"Number"=cell_num)
cell_use <- lapply(X = res_list,FUN = function(x)data.frame("Cell"=names(x),"Pred"=x))
cell_all <- do.call(what = rbind,args = cell_use)

ddd <- data.frame(table(cell_all))
ddd3 <- ddd[ddd$Freq>=75,]
ddd4 <- merge(x = ddd3,y = mic_test,by.x=1,by.y=0,all.x=T)
rrr <- data.frame("True"=ddd4[,"Subtype"],"Pred" = ddd4$Pred)
ddd5 <- data.frame(as.matrix(table(rrr)))
percent1 <- c()
for (i in 1:nrow(ddd5)){
    percent1 <- c(percent1,ddd5[i,"Freq"]/(cell_num[ddd5[i,"True"]]+cell_num[ddd5[i,"Pred"]]))
}
ddd5$Percent <- percent1


list_use <- c()
for (i in 1:7){
    for (j in (i+1):8){
        node_i <- names(cell_num)[i]
        node_j <- names(cell_num)[j]
        list_use[[paste0(i,"_",j)]] <- c(node_i,node_j,ddd5[ddd5$True == node_i & ddd5$Pred == node_j,"Percent"]+ddd5[ddd5$True == node_j & ddd5$Pred == node_i,"Percent"])
    }
}
edges2 <- data.frame(t(data.frame(list_use)))

#edges2 <- edges[!edges$True == edges$Pred,]

nodes <- data.frame(cell_num)

write.csv(x = nodes,file = "./classification/class2/cell_number.csv",quote = F,row.names = F)
write.csv(x = edges2,file = "./classification/class2/edges.csv",quote = F,row.names = F)




#########################
test_cells <- read.csv("./classification/class1/test_cells.csv",row.names = 1)
meta_data2 <- read.csv(file = "./classification/class1/metadata.csv",row.names = 1)
mic_test <- data.frame(row.names = test_cells$x,"Subtype"=meta_data2[test_cells$x,])
res_list <- readRDS(file = "./classification/class1/res_list.rds")
cells_use <- paste0("MG",c(1,3:9))
cell_num <- table(mic_test[,"Subtype"])
cell_num2 <- data.frame("cluster"=names(cell_num),"Number"=cell_num)
cell_use <- lapply(X = res_list,FUN = function(x)data.frame("Cell"=names(x),"Pred"=x))
cell_all <- do.call(what = rbind,args = cell_use)

ddd <- data.frame(table(cell_all))
ddd3 <- ddd[ddd$Freq>=75,]
ddd4 <- merge(x = ddd3,y = mic_test,by.x=1,by.y=0,all.x=T)
rrr <- data.frame("True"=ddd4[,"Subtype"],"Pred" = ddd4$Pred)
ddd5 <- data.frame(as.matrix(table(rrr)))
percent1 <- c()
for (i in 1:nrow(ddd5)){
    percent1 <- c(percent1,ddd5[i,"Freq"]/(cell_num[ddd5[i,"True"]]+cell_num[ddd5[i,"Pred"]]))
}
ddd5$Percent <- percent1


list_use <- c()
for (i in 1:7){
    for (j in (i+1):8){
        node_i <- names(cell_num)[i]
        node_j <- names(cell_num)[j]
        list_use[[paste0(i,"_",j)]] <- c(node_i,node_j,ddd5[ddd5$True == node_i & ddd5$Pred == node_j,"Percent"]+ddd5[ddd5$True == node_j & ddd5$Pred == node_i,"Percent"])
    }
}
edges2 <- data.frame(t(data.frame(list_use)))

#edges2 <- edges[!edges$True == edges$Pred,]

nodes <- data.frame(cell_num)

write.csv(x = nodes,file = "./classification/class1/cell_number.csv",quote = F,row.names = F)
write.csv(x = edges2,file = "./classification/class1/edges.csv",quote = F,row.names = F)

node_label <- paste0("MG",c(1,3:9))
group_type <- c("Alternative_homeostatic","Alternative_homeostatic","Alternative_homeostatic","Classical_homeostatic","DAM-like","DAM-like","DAM-like","Classical_homeostatic")
shape_seq <- rep(x = "circle",8)
value_seq <- as.numeric(cell_num)
color_seq <- c("#ffdd00","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
node_info <- data.frame(id = node_label,
                        label = node_label,
                        group = group_type,
                        value = value_seq,
                        #shape = shape_seq,
                        color = color_seq)
edges3 <- edges2[edges2$X3>=0.01,]
edge_info <- data.frame(from = edges3$X1,
                        to = edges3$X2,
                        smooth = F,
                        color = "grey",
                        width = as.numeric(edges3$X3) * 100)
visNetwork(node_info, edge_info,physics=T, idToLabel=F)




############################
test_cells <- read.csv("./classification/class3/test_cells.csv",row.names = 1)
meta_data2 <- read.csv(file = "./classification/class3/metadata.csv",row.names = 1)
mic_test <- data.frame(row.names = test_cells$x,"Subtype"=meta_data2[test_cells$x,])
res_list <- readRDS(file = "./classification/class3/res_list.rds")
cells_use <- paste0("MG",c(1,3:9))
cell_num <- table(mic_test[,"Subtype"])
cell_num2 <- data.frame("cluster"=names(cell_num),"Number"=cell_num)
cell_use <- lapply(X = res_list,FUN = function(x)data.frame("Cell"=names(x),"Pred"=x))
cell_all <- do.call(what = rbind,args = cell_use)

ddd <- data.frame(table(cell_all))
ddd3 <- ddd[ddd$Freq>=75,]
ddd4 <- merge(x = ddd3,y = mic_test,by.x=1,by.y=0,all.x=T)
rrr <- data.frame("True"=ddd4[,"Subtype"],"Pred" = ddd4$Pred)
ddd5 <- data.frame(as.matrix(table(rrr)))
percent1 <- c()
for (i in 1:nrow(ddd5)){
    percent1 <- c(percent1,ddd5[i,"Freq"]/(cell_num[ddd5[i,"True"]]+cell_num[ddd5[i,"Pred"]]))
}
ddd5$Percent <- percent1


list_use <- c()
for (i in 1:7){
    for (j in (i+1):8){
        node_i <- names(cell_num)[i]
        node_j <- names(cell_num)[j]
        list_use[[paste0(i,"_",j)]] <- c(node_i,node_j,ddd5[ddd5$True == node_i & ddd5$Pred == node_j,"Percent"]+ddd5[ddd5$True == node_j & ddd5$Pred == node_i,"Percent"])
    }
}
edges2 <- data.frame(t(data.frame(list_use)))

#edges2 <- edges[!edges$True == edges$Pred,]

nodes <- data.frame(cell_num)

write.csv(x = nodes,file = "./classification/class3/cell_number.csv",quote = F,row.names = F)
write.csv(x = edges2,file = "./classification/class3/edges.csv",quote = F,row.names = F)
library(visNetwork)

###################################
node_label <- paste0("MG",c(1,3:9))
group_type <- c("Alternative_homeostatic","Alternative_homeostatic","Alternative_homeostatic","Classical_homeostatic","DAM-like","DAM-like","DAM-like","Classical_homeostatic")
shape_seq <- rep(x = "circle",8)
value_seq <- as.numeric(cell_num)
library(igraph)
degree(edge2)
color_seq <- c("#ffdd00","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
node_info <- data.frame(id = node_label,
                        label = node_label,
                        group = group_type,
                        value = value_seq,
                        #shape = shape_seq,
                        color = color_seq)
edges3 <- edges2[edges2$X3>=0.01,]
edge_info <- data.frame(from = edges3$X1,
                        to = edges3$X2,
                        smooth = F,
                        color = "grey",
                        width = as.numeric(edges3$X3) * 100)
visNetwork(node_info, edge_info,physics=T, idToLabel=F)




library(visNetwork)
nodes <- data.frame(id = 1:10, 
                    label = paste("Node", 1:10),                                 # labels
                    group = c("GrA", "GrB"),                                     # groups 
                    value = 1:10+0.5,                                                # size 
                    
                    title = paste0("<p><b>", 1:10,"</b><br>Node !</p>"),         # tooltip
                    color = c("darkred", "grey", "orange", "darkblue", "purple"),# color
                    shadow = c(FALSE, TRUE, FALSE, TRUE, TRUE))                  # shadow

edges <- data.frame(from = sample(1:10,8), to = sample(1:10, 8),
                    label = paste("Edge", 1:8),                                 # labels
                    length = c(100,500),                                        # length
                    arrows = c("to", "from", "middle", "middle;to"),            # arrows
                    dashes = c(TRUE, FALSE),                                    # dashes
                    title = paste("Edge", 1:8),                                 # tooltip
                    smooth = c(FALSE, TRUE),                                    # smooth
                    shadow = c(FALSE, TRUE, FALSE, TRUE))                       # shadow

visNetwork(nodes, edges, physics=T, idToLabel=T) 



################################################# classification new
###################33 
genes_class <- read.xlsx(xlsxFile = "./PAM DAM homeostatic genes.xlsx",sheet = 1)
#genes_class <- read.xlsx(xlsxFile = "/lustre/home/zhangxiaochang/project/guo/microgila/classification/class1/PAM DAM homeostatic genes.xlsx",sheet = 1)
#degs_1 <- read.csv("classification/class3/degs_3_0.25.csv")

features_use <- unique(c(genes_class$DAM,genes_class$Homeostatic))
#features_use <- features_use[!is.na(features_use)]
features_use <- intersect(x = features_use,y = rownames(microgila_clean@assays$RNA@counts))
Idents(microgila_clean) <- "anno_detail"
write.csv(x = features_use,file = "classification/class_final/feature_genes.csv")
microgila_clean2 <- subset(x = microgila_clean,idents = c("MG2"),invert = T)
data_use <- t(data.frame(microgila_clean2@assays$RNA@counts[features_use,]))

meta_data2 <- data.frame(row.names = gsub(pattern = "-",replacement = ".",x = rownames(microgila_clean2@meta.data)),
                         "Subtype"= factor(x = microgila_clean2@meta.data$anno_detail,levels = paste0("MG",c(1,3:9))))
data_use <- merge(x = data_use,y = meta_data2,by = 0)
rownames(data_use) <- data_use[,1]
data_use <- data_use[,-1]
colnames(data_use) <- gsub(pattern = "-",replacement = ".",x = colnames(data_use))
colnames(data_use) <- gsub(pattern = "^[0-9]",replacement = "A",x = colnames(data_use))

library(randomForest)
train_cells <- c()
cells_use <- paste0("MG",c(1,3:9))
for(i in cells_use){
    cell_data <- data_use[data_use[["Subtype"]] == i,]
    cell_tmp <- rownames(cell_data)
    samlle_tmp <- sample(cell_tmp, length(cell_tmp)*0.75,replace = F)
    train_cells <- c(train_cells,samlle_tmp)
}

test_cells <- setdiff(x = rownames(data_use),y = train_cells)
write.csv(train_cells,file = "classification/class_final/train_cells.csv")
write.csv(test_cells,file = "classification/class_final/test_cells.csv")
write.csv(meta_data2,file = "classification/class_final/metadata.csv")



mic_train <- data_use[train_cells, ]
mic_test <- data_use[test_cells, ]
mic_train.forest <- randomForest(Subtype~., data = mic_train, importance = TRUE)
saveRDS(mic_train.forest,file = "trained_forest_model.rds")
mic_predict <- predict(mic_train.forest, mic_test)
# mic_predict <- predict(mic_train.forest, mic_train)
#### 统计正确率
rrr <- data.frame("True"=mic_test[,"Subtype"],"Pred" = mic_predict)
ddd <- data.frame(as.matrix(table(rrr)))

cells_use <- paste0("MG",c(1,3:9))
cell_num <- table(mic_test[,"Subtype"])
cell_num2 <- data.frame("cluster"=names(cell_num),"Number"=cell_num)
percent1 <- c()
for (i in 1:nrow(ddd)){
    percent1 <- c(percent1,ddd[i,"Freq"]/(cell_num[ddd[i,"True"]]+cell_num[ddd[i,"Pred"]]))
}
ddd$Percent <- percent1
write.csv(x = ddd,file = "classification/class_final/predict_percentage.csv")

true_cell <- c()
for (i in rownames(mic_test)){
    if (mic_test[i,"Subtype"] == mic_predict[i]){
        true_cell <- c(true_cell,i)
    }
}

cell_num
list_use <- c()
for (i in 1:7){
    for (j in (i+1):8){
        node_i <- names(cell_num)[i]
        node_j <- names(cell_num)[j]
        list_use[[paste0(i,"_",j)]] <- c(node_i,node_j,ddd[ddd$True == node_i & ddd$Pred == node_j,"Percent"]+ddd[ddd$True == node_j & ddd$Pred == node_i,"Percent"])
    }
}
edges2 <- data.frame(t(data.frame(list_use)))

#edges2 <- edges[!edges2$True == edges2$Pred,]

nodes <- data.frame(cell_num)
#edges <- ddd[,-3]
write.csv(x = nodes,file = "./classification/class_final//cell_number.csv",quote = F,row.names = F)
write.csv(x = edges2,file = "./class_final//edges.csv",quote = F,row.names = F)

library(visNetwork)

###################################
node_label <- paste0("MG",c(1,3:9))
group_type <- c("Alternative_homeostatic","Alternative_homeostatic","Alternative_homeostatic","Classical_homeostatic","DAM-like","DAM-like","DAM-like","Classical_homeostatic")
shape_seq <- rep(x = "circle",8)
value_seq <- as.numeric(cell_num)

color_seq <- c("#ffdd00","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
node_info <- data.frame(id = node_label,
                        label = node_label,
                        group = group_type,
                        value = value_seq,
                        #shape = shape_seq,
                        color = color_seq)
edges3 <- edges2[edges2$X3>=0.01,] ################### 过滤 < 0.01的边
edge_info <- data.frame(from = edges3$X1,
                        to = edges3$X2,
                        smooth = F,
                        color = "grey",
                        width = as.numeric(edges3$X3) * 100)
visNetwork(node_info, edge_info,physics=T, idToLabel=F)

##################
data_imp <- data.frame(importance(mic_train.forest)) #data.frame(mic_train.forest$importance)
# data_imp1 <- data_imp[order(data_imp$MeanDecreaseGini,decreasing = T),]
data_imp2 <- data_imp[order(data_imp$MeanDecreaseAccuracy,decreasing = T),]
write.csv(x = data_imp2,file = "importance.csv")
#data_imp_SD <- data.frame(mic_train.forest[["IncNodePurity"]])
num0 <- 50
# genes_use <- unique(c(rownames(data_imp1[1:num0,]),rownames(data_imp2[1:num0,])))
# genes_use <- unique(rownames(data_imp1[1:num0,]))
genes_use <- unique(rownames(data_imp2[1:num0,]))

genes_use <- gsub(pattern = "\\.",replacement = "-",x = genes_use)
genes_use <- gsub(pattern = "A632427E13Rik",replacement = "4632427E13Rik",x = genes_use)
microgila_clean
data_plot <- AverageExpression(object = microgila_clean2,assays = "RNA",features = genes_use,group.by = c("anno_detail"))$RNA

pdf(file = "importance.pdf",width = 6,height = 12)
randomForest::varImpPlot(mic_train.forest,n.var = 50,type = 1 ,main = "Top 50 important genes in randomForest model")
dev.off()
bk <- c(seq(-2,2,by=0.01))
pheatmap(mat = data_plot,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("navy", "white"))(200), colorRampPalette(c("white","firebrick3"))(200)),
         cluster_cols = T,cluster_rows = T,scale = "row",
         filename = "heatmap_importance_genes.pdf",width = 8,height = 10,
         #annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Gene expression of top 50 important genes in randomForest model",gaps_col = c(2,5))

###################
genes_class <- read.xlsx(xlsxFile = "./PAM DAM homeostatic genes.xlsx",sheet = 1)
#genes_class <- read.xlsx(xlsxFile = "/lustre/home/zhangxiaochang/project/guo/microgila/classification/class1/PAM DAM homeostatic genes.xlsx",sheet = 1)
#degs_1 <- read.csv("classification/class3/degs_3_0.25.csv")
genes_use <- unique(c(genes_class$DAM,genes_class$Homeostatic))
genes_use1 <- unique(c(genes_class$DAM))
genes_use2 <- unique(c(genes_class$Homeostatic))

data_plot1 <- data.frame(AverageExpression(object = microgila_clean2,assays = "RNA",features = genes_use1,group.by = c("anno_detail"))$RNA)
data_plot1$Mean <- rowMeans(data_plot1)
data_plot1 <- data_plot1[order(data_plot1$Mean,decreasing = T),]
data_plot1 <- data_plot1[1:50,-9]
write.csv(x = data_plot1,file = "heatmap_highly_expressed_DAM.csv")


data_plot2 <- data.frame(AverageExpression(object = microgila_clean2,assays = "RNA",features = genes_use2,group.by = c("anno_detail"))$RNA)
data_plot2$Mean <- rowMeans(data_plot2)
data_plot2 <- data_plot2[order(data_plot2$Mean,decreasing = T),]
data_plot2 <- data_plot2[1:50,-9]
write.csv(x = data_plot2,file = "heatmap_highly_expressed_Homeostatic.csv")

pheatmap(mat = data_plot1,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("navy", "white"))(200), colorRampPalette(c("white","firebrick3"))(200)),
         cluster_cols = T,cluster_rows = T,scale = "row",
         filename = "heatmap_highly_expressed_DAM.pdf",width = 8,height = 10,
         #annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Top 50 highly expressed genes in DAM gene set")

pheatmap(mat = data_plot2,#display_numbers = round(data_p2,2),
         breaks = bk,color = c(colorRampPalette(c("navy", "white"))(200), colorRampPalette(c("white","firebrick3"))(200)),
         cluster_cols = T,cluster_rows = T,scale = "row",
         filename = "heatmap_highly_expressed_Homeostatic.pdf",width = 8,height = 10,
         #annotation_col = anno_group,annotation_colors = anno_color2,
         main = "Top 50 highly expressed genes in Homeostatic gene set")



###############################################3 classification final
#########################
test_cells <- read.csv("./classification/class_final/train_100/test_cells.csv",row.names = 1)
meta_data2 <- read.csv(file = "./classification/class_final/train_100/metadata.csv",row.names = 1)
mic_test <- data.frame(row.names = test_cells$x,"Subtype"=meta_data2[test_cells$x,])
res_list <- readRDS(file = "./classification/class_final/train_100/res_list.rds")
cells_use <- paste0("MG",c(1,3:9))
cell_num <- table(mic_test[,"Subtype"])
cell_num2 <- data.frame("cluster"=names(cell_num),"Number"=cell_num)
cell_use <- lapply(X = res_list,FUN = function(x)data.frame("Cell"=names(x),"Pred"=x))
cell_all <- do.call(what = rbind,args = cell_use)

ddd <- data.frame(table(cell_all))
ddd3 <- ddd[ddd$Freq>=75,]
ddd4 <- merge(x = ddd3,y = mic_test,by.x=1,by.y=0,all.x=T)
rrr <- data.frame("True"=ddd4[,"Subtype"],"Pred" = ddd4$Pred)
ddd5 <- data.frame(as.matrix(table(rrr)))
percent1 <- c()
for (i in 1:nrow(ddd5)){
    percent1 <- c(percent1,ddd5[i,"Freq"]/(cell_num[ddd5[i,"True"]]+cell_num[ddd5[i,"Pred"]]))
}
ddd5$Percent <- percent1


list_use <- c()
for (i in 1:7){
    for (j in (i+1):8){
        node_i <- names(cell_num)[i]
        node_j <- names(cell_num)[j]
        list_use[[paste0(i,"_",j)]] <- c(node_i,node_j,ddd5[ddd5$True == node_i & ddd5$Pred == node_j,"Percent"]+ddd5[ddd5$True == node_j & ddd5$Pred == node_i,"Percent"])
    }
}
edges2 <- data.frame(t(data.frame(list_use)))

#edges2 <- edges[!edges$True == edges$Pred,]

nodes <- data.frame(cell_num)

write.csv(x = nodes,file = "./classification/class_final/train_100/nodes.csv",quote = F,row.names = F)
write.csv(x = edges2,file = "./classification/class_final/train_100/edges.csv",quote = F,row.names = F)

node_label <- paste0("MG",c(1,3:9))
group_type <- c("Alternative_homeostatic","Alternative_homeostatic","Alternative_homeostatic","Classical_homeostatic","DAM-like","DAM-like","DAM-like","Classical_homeostatic")
shape_seq <- rep(x = "circle",8)
value_seq <- as.numeric(cell_num)
color_seq <- c("#ffdd00","#ff0000","#ff0092","#c1d82f","#6a67ce","#00a4e4","#0863b5","#f67019")
node_info <- data.frame(id = node_label,
                        label = node_label,
                        group = group_type,
                        value = value_seq,
                        #shape = shape_seq,
                        color = color_seq)
edges3 <- edges2[edges2$X3>=0.01,]
edge_info <- data.frame(from = edges3$X1,
                        to = edges3$X2,
                        smooth = F,
                        color = "grey",
                        width = as.numeric(edges3$X3) * 100)
visNetwork(node_info, edge_info,physics=T, idToLabel=F)

####################################################################
###############比例统计

seurat_anno <- readRDS(file = "./seurat_anno_final.rds")
microgila_anno <- subset(seurat_anno,idents="Microgila")
Idents(microgila_anno) <- "anno_detail"
anno_info <- data.frame(sample=microgila_anno$re_names,
                        cell_type=microgila_anno@active.ident)

cell_number <- data.frame(table(anno_info[c("sample","cell_type")]))
cell_num0 <- data.frame(table(anno_info[c("sample")]))
cell_number[["sample_num"]] <- rep(cell_num0$Freq,9)
cell_number[["percentage"]] <- round(cell_number$Freq/cell_number$sample_num*100,2)
cell_number$sample <- factor(cell_number$sample,levels = names(col_model))

cell_names <- unique(cell_number$cell_type)
write.csv(x = cell_number,file = "percentage_in_microgila.csv")

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
        ggtitle(paste0(cell_name," in microgila cells"))
}


multiplot(plotlist = list_pics,cols=3)

#############差异基因
library(ggplot2)
library(ggrepel)

Idents(object = microgila_anno) <- "anno_abbr"
deg_mg7 <- FindMarkers(object = microgila_anno,ident.1 = "MG7",assay = "RNA",logfc.threshold = 0,min.pct = 0.1)
deg_mg7 <- deg_mg7[deg_mg7$p_val_adj<0.05,]
deg_mg8 <- FindMarkers(object = microgila_anno,ident.1 = "MG8",assay = "RNA",logfc.threshold = 0,min.pct = 0.1)
deg_mg8 <- deg_mg8[deg_mg8$p_val_adj<0.05,]
write.csv(deg_mg7,"deg_mg7.csv")
write.csv(deg_mg8,"deg_mg8.csv")
deg_mg7 <- read.csv("./deg_mg7.csv",header = T,row.names = 1)
deg_mg8 <- read.csv("./deg_mg8.csv",header = T,row.names = 1)


Dat <- deg_mg7
#确定是上调还是下调，用于给图中点上色）
Dat$threshold = factor(ifelse(Dat$p_val_adj < 0.05 & abs(Dat$avg_log2FC) >= 0.25, ifelse(Dat$avg_log2FC>= 0.25 ,'Up','Down'),'NoSignifi'),levels=c('Up','Down','NoSignifi'))
Dat$Gene = rownames(Dat)
############R语言绘图超过-300次方不显示
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
        legend.title = element_blank()#不显示图例标题
    )+
    ylab('-log10 (p_val_adj)')+#修改y轴名称
    xlab('avg_log2FC')+#修改x轴名称
    geom_vline(xintercept=c(-0.25,0.25),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
    geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05


######################富集分析

library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)

deg_mg7$gene <- rownames(deg_mg7)
deg_mg8$gene <- rownames(deg_mg8)
deg_mg7$cluster <- rep(x = "MG7",nrow(deg_mg7))
deg_mg8$cluster <- rep(x = "MG8",nrow(deg_mg8))
markers <- rbind(deg_mg7,deg_mg8)

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

list_pic_count$MG7
list_pic_count$MG8

list_pic_qvalue$MG7
list_pic_qvalue$MG8

#########################
microgila_anno <- readRDS(file = "./seurat_anno_final.rds")
Idents(microgila_anno) <- "anno_detail"
########## mg7
mg7 <- subset(x = microgila_anno,idents = c("MG7"))
DefaultAssay(mg7) <- "RNA"
genes_all <- rownames(mg7@assays$RNA)
mg7 <- ScaleData(object = mg7,features = genes_all)
Idents(mg7) <- "re_names"
markers_4 <- FindAllMarkers(object = mg7,assay = "RNA",logfc.threshold = 0.25,min.pct = 0.1)
# markers_4 <- markers_4[markers_4$p_val_adj<0.05,]
# markers_4_up <- markers_4[markers_4$avg_log2FC>0,]
# markers_4_up <- markers_4_up[order(markers_4_up$avg_log2FC,decreasing = T),]
# markers_4_down <- markers_4[markers_4$avg_log2FC<0,]
# markers_4_down <- markers_4_down[order(markers_4_down$avg_log2FC,decreasing = F),]

markers_4_Con <- FindMarkers(object = mg7,ident.1 = "Con")
markers_4_Con$cluster <- "Con"
markers_4_Con$gene <- rownames(markers_4_Con)

markers_4_LPS <- FindMarkers(object = mg7,ident.1 = "LPS",ident.2 = "Con")
markers_4_LPS$cluster <- "LPS"
markers_4_LPS$gene <- rownames(markers_4_LPS)

markers_4_HYP <- FindMarkers(object = mg7,ident.1 = "HYP",ident.2 = "Con")
markers_4_HYP$cluster <- "HYP"
markers_4_HYP$gene <- rownames(markers_4_HYP)

markers_4_LH <- FindMarkers(object = mg7,ident.1 = "L/H",ident.2 = "Con")
markers_4_LH$cluster <- "L/H"
markers_4_LH$gene <- rownames(markers_4_LH)

markers_4 <- do.call(what = rbind,args = list(markers_4_Con,markers_4_LPS,markers_4_HYP,markers_4_LH))
write.csv(markers_4,"./DEGs_mg7_group_markers.csv")

markers_4 <- markers_4[markers_4$p_val_adj < 0.05,]
markers_4_up <- markers_4[markers_4$avg_log2FC>0.25,]
markers_4_up <- markers_4_up[order(markers_4_up$avg_log2FC,decreasing = T),]
markers_4_down <- markers_4[markers_4$avg_log2FC < -0.25,]
markers_4_down <- markers_4_down[order(markers_4_down$avg_log2FC,decreasing = F),]



top_5_up <- unique(unlist(lapply(X = c("LPS","HYP","L/H"),FUN = function(x){
    y <- markers_4_up[markers_4_up$cluster==x,"gene"]
    y <- y[c(1:10)]
})))
top_5_down <- unique(unlist(lapply(X = c("Con","LPS","HYP","L/H"),FUN = function(x){
    y <- markers_4_down[markers_4_down$cluster==x,"gene"]
    y <- y[c(1:5)]
})))
DoHeatmap(object = mg7,assay = "RNA",features = top_5_up,group.by = "re_names")
exp_ave <- AverageExpression(object = mg7,assays = "RNA",features = top_5_up)
DoHeatmap(object = mg7,assay = "RNA",features = top_5_down,group.by = "re_names")
pheatmap(mat = exp_ave$RNA,scale = "row",cluster_cols = F,cluster_rows = F)

markers_4_up_list <- list()
for (i in c("Con","LPS","HYP","L/H")){
    markers_4_up_list[[i]] <- markers_4_up[markers_4_up$cluster==i,"gene"]
}
markers_4_down_list <- list()
for (i in c("Con","LPS","HYP","L/H")){
    markers_4_down_list[[i]] <- markers_4_down[markers_4_down$cluster==i,"gene"]
}
library(clusterProfiler)
library(org.Mm.eg.db)
res_up <- compareCluster(geneClusters = markers_4_up_list,fun = enrichGO,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont="BP")
res_down <- compareCluster(geneClusters = markers_4_down_list,fun = enrichGO,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont="BP")

dotplot(res_up)
dotplot(res_down)
write.csv(res_up,"./res_up.csv")
write.csv(res_down,"./res_down.csv")

###########

mg8 <- subset(x = microgila_anno,idents = c("MG8"))
DefaultAssay(mg8) <- "RNA"
genes_all <- rownames(mg8@assays$RNA)
mg8 <- ScaleData(object = mg8,features = genes_all)
mg8$re_names <- factor(x = mg8$re_names,levels = c("Con","LPS","HYP","L/H"))
Idents(mg8) <- "re_names"
markers_4 <- FindAllMarkers(object = mg8,assay = "RNA",logfc.threshold = 0.25,min.pct = 0.1)
write.csv(x = markers_4,file = "./DEGs_mg8_group_markers.csv")
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
    y <- y[c(1:10)]
})))
DoHeatmap(object = mg8,assay = "RNA",features = top_5_up,group.by = "re_names")
exp_ave <- AverageExpression(object = mg8,assays = "RNA",features = top_5_up)
DoHeatmap(object = mg8,assay = "RNA",features = top_5_down,group.by = "re_names")
pheatmap(mat = exp_ave$RNA,scale = "row",cluster_cols = F,cluster_rows = F)

markers_4_up_list <- list()
for (i in c("Con","LPS","HYP","L/H")){
    markers_4_up_list[[i]] <- markers_4_up[markers_4_up$cluster==i,"gene"]
}
markers_4_down_list <- list()
for (i in c("Con","LPS","HYP","L/H")){
    markers_4_down_list[[i]] <- markers_4_down[markers_4_down$cluster==i,"gene"]
}
library(clusterProfiler)
library(org.Mm.eg.db)
res_up <- compareCluster(geneClusters = markers_4_up_list,fun = enrichGO,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont="BP")
res_down <- compareCluster(geneClusters = markers_4_down_list,fun = enrichGO,OrgDb=org.Mm.eg.db,keyType="SYMBOL",ont="BP")

dotplot(res_up,title="Up Regulation")+theme(plot.title = element_text(hjust = 0.5)) 
dotplot(res_down,title="Down Regulation")+theme(plot.title = element_text(hjust = 0.5))
write.csv(res_up,"./res_up_mg8.csv")
write.csv(res_down,"./res_down_mg8.csv")



