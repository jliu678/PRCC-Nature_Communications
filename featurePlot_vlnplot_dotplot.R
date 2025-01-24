library(Seurat, lib.loc = "/opt/R/4.4.0/lib/R/library")
library('dplyr')
library('magrittr')
library('ggplot2')
library('data.table')

#4. plot features----
readRDS('kidney.Rds')->kid
kid@meta.data<-readRDS('kidneyMetaDt.rds')
kid[[]] %>% View()
kid@images<-list()

f1<-c('Rida','Cdh16','Lrp2','Slc4a4','Tert',
      'Slc12a1','Hsd11b2','Kdr','Mgp',
      'Ptprc','Itgam','Mpeg1','S100a8',
      'Cd79a','Thy1','Ly6g','Ly6c1','Ly6c2')

f2<-c('Birc5',
      'Vim',
      'Krt7',
      'Pax8',
      'Mdm2',
      'Ercc6',
      'Mme',#cd10
      'Pax2',
      'Muc1',
      'Kras',
      'Ska3',
      'Mki67',
      'Tert'
)

f3<-c('Yap1',
      'Alpl',
      'Akr1c21',
      'Arg2',
      'Ass1',
      'Cldn1',
      'Ankrd1',
      'Ctgf',
      'Cyr61',
      'Lats1',
      'Lats2',
      'Hbegf',
      'Egfr',
      'Stmn1',
      'Erbb2',
      'Csf1',
      'Cxcl2',
      'Csf1r',
      'Csf2ra',
      'Lyz2',
      'Fcer1g'
)
#check what Ident is in use
Idents(kid)

FeaturePlot(kid %>% 
              subset(idents = c("5","11")),
            features = f1,
            order = T)

FeaturePlot(kid %>% 
              subset(idents = c("5","11")),
            features = f2,
            order = T)

DotPlot(kid,
        cols ='PiYG',
        # features = f1,
        features = 'Fst',
        # group.by = 'leiden.RBConfigurationVertexPartition',
        group.by = 'laipingCluster',
        split.by = 'group'
)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

kid@active.ident<-factor(kid@active.ident,
                         levels = c('1','3','4','6',
                                    '15','5','11','7',
                                    '14','9','17','0',
                                    '8','13','2','10',
                                    '16','12'))

cls<-setNames(c('#CD0000','#006100','#CD5B5B',
                '#248347','#0FD871','#7D26CD',
                '#3DFD3F','#8AF6A2','#FF0000',
                '#0000FF','#FA8072','#EE81EE',
                '#FF68B3','#FF4602','#40E0CF',
                '#7DFB25','#FFC3CE','#80CAEA'),
              nm = as.character(0:17))


VlnPlot(kid,
        features = c(f1),
        stack = T,
        flip=T,
        # group.by = 'leiden.RBConfigurationVertexPartition',
        # split.by = 'group', #ok,well
        sort = F,
        cols = cls[levels(Idents(kid))],
        fill.by = "ident"
) +theme(legend.position = 'none')->p1
