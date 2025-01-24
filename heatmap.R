dyn.load("/apps/lib-osver/geos/3.5.0/lib/libgeos_c.so.1")
library(Seurat)
library(sctransform)

library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)

#a. get cluster markers---- 
readRDS(file = 'deg_cluster5and11AcrossTime_fig5.rds')->tt_deg
#remove those with p>0.05
tt_deg %<>% lapply(function(i){
  i[i$p_val<0.05,] %>% arrange(desc(avg_log2FC))
})

#b.plot----
readRDS('kidney.Rds')->kid
kid@meta.data<-readRDS('kidneyMetaDt.rds')
kid[[]] %>% View()
DefaultAssay(kid)<-'RNA'#default and 'integrated' show contrary trend
#and seurat suggest using 'RNA' for quantitive comparsion
#https://github.com/satijalab/seurat/issues/1717

#...b0 data prep----
# ...b0.1----
#only keep the most significant genes

names(tt_deg)[1:3] %>% sapply(function(i){
  tt_deg[[i]]$p_val_adj<0.05 & abs(tt_deg[[i]]$avg_log2FC )>0.5->tl
  # write.csv(tt_deg[[i]][tl,],
  #           file = paste0('genesShownInHeatmap_',
  #                         i,'.csv'))
  setNames(object = tt_deg[[i]][tl,"avg_log2FC",drop=T],
           nm = rownames(tt_deg[[i]])[tl])
})->tf4

vector(mode='list',length=length(c('CTR','D7','M4'))) %>%
  setNames(nm = c('CTR','D7','M4'))->tl

# ...b0.2----
names(tf4)[1] %>%
  # names(tf4)[2] %>%
  # names(tf4) %>%
  sapply(function(i){
    gsub(i,pattern = '^.*__|minus.*$',replacement = '')->ttt1
    cat(ttt1,1,'\n')
    tl[[ttt1]]<<-c(tl[[ttt1]],names(tf4[[i]][tf4[[i]]>0]))
    
    
    gsub(i,pattern = '^.*minus',replacement = '')->ttt1
    cat(ttt1,2,'\n') 
    tl[[ttt1]]<<-c(tl[[ttt1]],names(tf4[[i]][tf4[[i]]<0]))
    
    cat('done',i,'\n')
  })

sum(lengths(tl))
tl %>% unlist %>% unique()->f4
rm(tf4,tl)

#...b1 complexheat----
# devtools::install_github("jokergoo/ComplexHeatmap")
library('ComplexHeatmap')

readRDS(file = 'cells1_and_cell2_for_seuratFindmarker.rds')->cells4findmarker

#...b1.1 individual combined----

kid@assays$RNA@data[f4,c(
  cells4findmarker$`5__D7minusCTR`$js,
  cells4findmarker$`5__D7minusCTR`$bjs,
  cells4findmarker$`5__M4minusCTR`$bjs
)] %>% 
  as.matrix %>% t %>% scale %>% t -> d4heatmap



kid@assays$RNA@data[f4,c(
  cells4findmarker$`11__D7minusCTR`$js,
  cells4findmarker$`11__D7minusCTR`$bjs,
  cells4findmarker$`11__M4minusD7`$bjs
)] %>% 
  as.matrix %>% t %>% scale %>% t -> d4heatmap

c(
  length(cells4findmarker$`5__D7minusCTR`$js),
  length(cells4findmarker$`5__D7minusCTR`$bjs),
  length(cells4findmarker$`5__M4minusCTR`$bjs))->lths
grNm<-c('CTR','D7','M4')

c(
  length(cells4findmarker$`11__D7minusCTR`$js),
  length(cells4findmarker$`11__D7minusCTR`$bjs),
  length(cells4findmarker$`11__M4minusCTR`$bjs))->lths
grNm<-c('CTR','D7','M4')


#order rows
data.table(t(d4heatmap))[, lapply(.SD, mean), by = list(rep(grNm,lths))]->d4heatmap_t
d4heatmap_t[,2:ncol(d4heatmap_t)] %>%
  apply(MARGIN = 2,rank) %>% t->d4heatmap_t1
colnames(d4heatmap_t1)<-d4heatmap_t$rep

d4heatmap_t1 %<>% as.data.frame %>% arrange(CTR,D7)

d4heatmap_t1$rowNm<-rownames(d4heatmap_t1)
d4heatmap_t1 %<>% group_by(CTR) %>% sample_frac %>% as.data.frame

d4heatmap_t1<-rbind(
  d4heatmap_t1 %>% filter(CTR=='1'),
  d4heatmap_t1 %>% filter(CTR=='2') %>% arrange(desc(D7)),
  d4heatmap_t1 %>% filter(CTR=='3')
)


rownames(d4heatmap_t1)<-d4heatmap_t1$rowNm
d4heatmap_t1$rowNm<-NULL

d4heatmap_t1 %>%
  Heatmap(
    name=' ',
    column_title ='preview of heatmap',
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_dend = F,
    show_row_names = F,
    show_column_names = T)

d4heatmap<-d4heatmap[rownames(d4heatmap_t1),]

ha1<- HeatmapAnnotation(
  bar=rep(grNm[1],lths[1]),
  annotation_label='',
  col = list(
    bar=setNames(object = c("grey"),
                 nm = grNm[1])
  ),
  show_legend = F
)

ha2<- HeatmapAnnotation(
  bar=rep(grNm[2],lths[2]),
  annotation_label='',
  col = list(
    bar=setNames(object = c("red"),
                 nm = grNm[2])
  ),
  show_legend = F
)

ha3<- HeatmapAnnotation(
  bar=rep(grNm[3],lths[3]),
  annotation_label='',
  col = list(
    bar=setNames(object = c("purple"),
                 nm = grNm[3])
  ),
  show_legend = F
)

lapply(list.files(pattern = '_heatmap_highlighted\\.xlsx'),function(i){
  readxl::read_xlsx(path =i ,sheet = 'selected')
})->f1

names(f1)<-sapply(f1, function(i) {
  colnames(i)[1]
})

f1<-f1$`C5-D7/M4`$...2
f1<-f1$`C11- D7/M4`$...2

# annomark1<-rownames(d4heatmap_t1)[d4heatmap_t1$CTR=='1' | 
#                                     (d4heatmap_t1$CTR=='2' & d4heatmap_t1$D7=='1')]
# 
# annomark2<-rownames(d4heatmap_t1)[d4heatmap_t1$CTR=='3' | 
#                                     (d4heatmap_t1$CTR=='2' & d4heatmap_t1$D7=='3')]

annomark1<-rownames(d4heatmap_t1)[d4heatmap_t1$CTR=='1']

annomark2<-rownames(d4heatmap_t1)[d4heatmap_t1$CTR!='1']


sample(intersect(f1,annomark1),size = 0)->tt1


hb1<-rowAnnotation(
  foo = anno_mark(at = match(table = rownames(d4heatmap),
                             x = setdiff(intersect(f1,annomark1),tt1)), 
                  labels = setdiff(intersect(f1,annomark1),tt1),
                  side = 'left',
                  labels_gp = gpar(
                    # col = "white",
                    # fontsize = 1
                  ),
                  link_gp = gpar(alpha=0.5)),
  foo1 = anno_mark(at = match(table = rownames(d4heatmap),
                              x = c(intersect(f1,annomark2),tt1)), 
                   labels = c(intersect(f1,annomark2),tt1),
                   side = 'left',
                   labels_gp = gpar(
                     # col = "white",
                     # fontsize = 1
                   ),
                   link_gp = gpar(alpha=0.5))
)


# d4heatmap[1:3,1:lths[1]] %>%
d4heatmap[,1:lths[1]] %>%
  Heatmap(
    name=' ',
    col = col_fun,
    cluster_columns = F,
    cluster_rows = F,
    column_title='CTR',
    show_row_dend = F,
    show_column_dend = F,
    show_row_names = F,
    show_column_names = F,
    column_dend_height = unit(0.4, "cm"),
    column_dend_side ="bottom",
    left_annotation = hb1,
    top_annotation = ha1,
    use_raster = T,#turn on for ppt
    heatmap_legend_param = list(
      title = "", at = c(-2, 0, 2), 
      labels = c('-2', '0', '2')),
    width = unit(2/4, "inch")
  )->p1

# d4heatmap[1:3,(1+lths[1]):sum(lths[1:2])] %>%
d4heatmap[,(1+lths[1]):sum(lths[1:2])] %>%
  Heatmap(
    name=' ',
    col = col_fun,
    cluster_columns = F,
    column_title='D7',
    show_row_dend = F,
    show_column_dend = F,
    show_row_names = F,
    show_column_names = F,
    column_dend_height = unit(0.4, "cm"),
    column_dend_side ="bottom",
    # left_annotation = ha1,
    top_annotation = ha2,
    use_raster = T,#turn on for ppt
    width = unit(5/3, "inch")
  )->p2

# d4heatmap[1:3,(1+sum(lths[1:2])):sum(lths)] %>%
d4heatmap[,(1+sum(lths[1:2])):sum(lths)] %>%
  Heatmap(
    name=' ',
    col = col_fun,
    cluster_columns = F,
    column_title='M4',
    show_row_dend = F,
    show_column_dend = F,
    show_row_names = F,
    show_column_names = F,
    column_dend_height = unit(0.4, "cm"),
    column_dend_side ="bottom",
    # left_annotation = ha1,
    top_annotation = ha3,
    use_raster = T,#turn on for ppt
    width = unit(1, "inch")
  )->p3

png(filename = paste0('fig5e_','heatmap.png'),res=300,
    width = 6.3,height = 4.331,units = 'in')

draw(p1+p2+p3)->p4

dev.off()
