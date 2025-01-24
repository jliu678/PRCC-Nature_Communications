library(Seurat, lib.loc = "/opt/R/4.4.0/lib/R/library")
library('dplyr')
library('magrittr')
library('ggplot2')
library('data.table')

#dimplot----
DimPlot(
  object = kid,
  dims = c(1, 2),
  cells = NULL,
  cols = cls[levels(Idents(kid))],
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  shuffle = T,
  seed = 1,
  label = T,
  label.size = 4.7,
  label.color = "white",
  label.box = T,
  repel = T,
  cells.highlight = NULL,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  na.value = "grey50",
  ncol = NULL,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
)+theme(legend.position = 'none')->p1

p1


cls<-setNames(c('#CD0000','#2d78a1','#CD5B5B',
                '#006100','#0FD871','#7D26CD',
                'grey','#8AF6A2','#FF0000',
                '#0000FF','#FA8072','#EE81EE',
                '#FF68B3','#FF4602','#40E0CF',
                '#f5e342','#FFC3CE','#80CAEA'),
              nm = as.character(0:17))


kid@active.ident<-factor(kid@active.ident,
                         levels = as.character(0:17))


DimPlot(
  object = kid,
  dims = c(1, 2),
  cells = NULL,
  cols = cls[levels(Idents(kid))],
  pt.size = NULL,
  reduction = NULL,
  group.by = NULL,
  split.by = NULL,
  shape.by = NULL,
  order = NULL,
  shuffle = T,
  seed = 1,
  label = T,
  label.size = 4.7,
  label.color = "black",
  label.box = F,
  repel = T,
  cells.highlight = NULL,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  na.value = "grey50",
  ncol = NULL,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
)+theme(
  # legend.position = 'none'
  legend.position = 'right'
)->p1
p1
ggsave(filename = 'fig4a_dimplot.new_color.png',plot = p1,width = 5.7,height = 5)