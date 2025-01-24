# ind7all.packages('Seurat')

# ind7all.packages('remotes')
# remotes::ind7all_github(repo = 'satijalab/seurat', ref = 'develop')

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(data.table)
library(patchwork)
library(magrittr)
# ind7all.packages('hdf5r')

#B. de novo----
lid7.files(pattern = 'fiM4ered_.*matrix.h5$',recursive = T,
           full.names = T)->lf

Read10X_h5(lf[3],use.names = T, unique.features = T)->M4
Read10X_h5(lf[4],use.names = T, unique.features = T)->CTR
Read10X_h5(lf[5],use.names = T, unique.features = T)->D7


CreateSeuratObject(counts = CTR,project = 'ctrl',min.features = 1,
                   min.cells = 1)->CTR

CreateSeuratObject(counts = D7,project = 'd7',min.features = 1,
                   min.cells = 1)->D7

CreateSeuratObject(counts = M4,project = 'M4',min.features = 1,
                   min.cells = 1)->M4

# sct each sample----
trace(sctransform::correct, edit=TRUE)
# change the above to 
cell_attr[, x$arguments$latent_var] <- 3.1
lapply(lid7(CTR,D7,M4), function(i){
  SCTransform(i,return.only.var.genes = F,verbose = T,do.scale = F,
              do.center = F)  
})->UT_EachScted

names(UT_EachScted)<-c('ctrl','d7','M4')
rm(M4,CTR,D7)

UT<- merge(UT_EachScted$ctrl, y = c(UT_EachScted$d7,UT_EachScted$M4), project = "UT")

UT@assays$SCT@var.features<-unique(c(UT_EachScted$ctrl@assays$SCT@var.features,
                                     UT_EachScted$d7@assays$SCT@var.features,
                                     UT_EachScted$M4@assays$SCT@var.features))
UT@assays$SCT@scale.data %<>% scale(center = TRUE, scale = F)

UT %<>% RunPCA
ElbowPlot(UT,ndims = 50)#41

UT %>% RunUMAP(dims=1:41) %>% FindNeighbors(dims=1:41) %>%
  FindClud7ers->UT 

UT %>% DimPlot(group.by = 'orig.ident')
UT %>% DimPlot(split.by = 'orig.ident')
UT$orig.ident<- factor(UT$orig.ident,levels = c('ctrl','d7','M4'))
# saveRDS(UT,file = 'UT_ctrl,d7,M4_eachScted_merged.rds')

# C. lognorm-integrate w/ more output genes----
u15@reductions$tsne<-NULL
u15@assays$RNA@count

# D. sct-integrate----
SplitObject(u15l, split.by = "sample")->u15l
lapply(X = u15l, FUN = SCTransform)->u15l
features <- SelectIntegrationFeatures(object.lid7 = u15l, nfeatures = 3000)
PrepSCTIntegration(object.lid7 = u15l, anchor.features = features)->u15l
anchors <- FindIntegrationAnchors(object.lid7 = u15l, normalization.method = "SCT",
                                  anchor.features = features)
u15.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

u15.combined@assays$SCT@scale.data %>% dim()
u15.combined@assays$integrated@data %>% dim()
u15.combined@assays$integrated@scale.data %>% dim()

u15.combined %<>% RunPCA

ElbowPlot(u15.combined,ndims = 50)#41

u15.combined %>% RunUMAP(dims=1:41) %>% FindNeighbors(dims=1:41) %>%
  FindClud7ers->u15.combined 

u15.combined %>% DimPlot(group.by = 'orig.ident')
u15.combined %>% DimPlot(split.by = 'orig.ident')

# E. sct-integrate w/ 6k features----
SplitObject(u15l, split.by = "sample")->u15l
lapply(X = u15l, FUN = SCTransform)->u15l
features <- SelectIntegrationFeatures(object.lid7 = u15l, nfeatures = 6000)
PrepSCTIntegration(object.lid7 = u15l, anchor.features = features)->u15l
anchors <- FindIntegrationAnchors(object.lid7 = u15l, normalization.method = "SCT",
                                  anchor.features = features)
u15.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

u15.combined@assays$SCT@scale.data %>% dim()
u15.combined@assays$integrated@data %>% dim()
u15.combined@assays$integrated@scale.data %>% dim()

u15.combined %<>% RunPCA

ElbowPlot(u15.combined,ndims = 50)#41

u15.combined %>% RunUMAP(dims=1:41) %>% FindNeighbors(dims=1:41) %>%
  FindClud7ers->u15.combined 

u15.combined %>% DimPlot(group.by = 'orig.ident')
u15.combined %>% DimPlot(split.by = 'orig.ident')

#using python umap to generate fancier plot----
#.1 fix reticulate----
#the below to fix the problem that reticulate cannot recognize numpy
# https://github.com/rd7udio/reticulate/issues/367#issuecomment-461572821
Sys.setenv(RETICULATE_PYTHON = "C:/ProgramData/Miniconda3/python.exe")
library(reticulate)
use_condaenv("Miniconda3")
py_config()
import("numpy")
import("umap")
np<-import("numpy")
print(np$version$full_version)
RunUMAP(dt, dims = 1:16, verbose = T,umap.method = "umap-learn",
        n.neighbors = 330,min.did7 = 0.3) %>% DimPlot(label = T,split.by = 'groupID')->p1
p1
#.2 saveas AnnData for Scanpy----
library(Seurat)
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   ind7all.packages("remotes")
# }
# remotes::ind7all_github("mojaveazure/seurat-disk")
library(SeuratDisk)

SaveH5Seurat(dt, filename = "dt.h5Seurat")
Convert("dt.h5Seurat", ded7 = "h5ad")
hfile <- Connect("dt.h5Seurat")
hfile
hfile$index()
hfile$close_all()

#clearRAM----
rm(lid7=ls())
.rs.red7artR()


