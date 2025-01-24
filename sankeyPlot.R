#plot the sankey-variant----
#A.prep data----
dyn.load("/apps/lib-osver/geos/3.5.0/lib/libgeos_c.so.1")
library(Seurat)
library(magrittr)
library(sctransform)
library(dplyr)
library(data.table)

readRDS('kidney.Rds')->kid
kid@meta.data<-readRDS('kidneyMetaDt.rds')
kid[[]] %>% View()

kid@meta.data %<>% mutate(cellType=case_when(
  laipingCluster %in% c(1,3,4,6,15) ~ 'Proximal tube epithelial cells',
  laipingCluster %in% c(5,11) ~ 'Cancer Cells',
  laipingCluster %in% c(7) ~ 'Distal tube epithelial cells',
  laipingCluster %in% c(14) ~ 'Collection tube',
  laipingCluster %in% c(9) ~ 'Endothelial cells',
  laipingCluster %in% c(17) ~ 'Fibroblast',
  laipingCluster %in% c(0,8,13,2) ~ 'MDSCs',
  laipingCluster %in% c(10) ~ 'Neutrophils',
  laipingCluster %in% c(16) ~ 'B cells',
  laipingCluster %in% c(12) ~ 'T Cells',
  T ~ 'somethingWrongIfYouSeeMe'
))

kid$cellType %>% unique()

kid@meta.data %>% 
  group_by(group,cellType,laipingCluster) %>% 
  summarise(nCellperCluster=n())->tt1

kid$group %>% table()->tt2

tt1 %<>% mutate(nCellperGroup=case_when(
  group == 'CTR' ~ tt2['CTR'],
  group == 'D7' ~ tt2['D7'],
  group == 'M4' ~ tt2['M4'],
  T ~ NA_integer_
))


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

sample(x = c25,size = 18)->c1
pie(rep(1, 18), col = c1)

tt1$cellType %>% unique()
tt1$cellType<- factor(x = tt1$cellType,
                      levels = c(
                        'Proximal tube epithelial cells',
                        'Distal tube epithelial cells',
                        'Collection tube',
                        'Cancer Cells',
                        'Endothelial cells',
                        'Fibroblast',
                        'MDSCs',
                        'Neutrophils',
                        'B cells',
                        'T Cells'
                      ))

library(ggalluvial)
ggplot(data = tt1,
       aes(x = group,
           y = percCluster2Group,
           alluvium = laipingCluster,
           stratum = laipingCluster
       )) +
  geom_alluvium(aes(fill = laipingCluster,
                    colour = laipingCluster),
                alpha = .75,
                show.legend=NA,
                decreasing = F) +
  geom_stratum(aes(fill = laipingCluster),
               decreasing = F)+
  geom_text(aes(label = laipingCluster),decreasing = F,
            stat = "stratum", size = 3,colour = "white")+
  # scale_x_continuous(breaks = seq(2003, 2013, 1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -0, hjust = 0.5)) +
  scale_fill_manual(values = c1)  +
  scale_colour_manual(values = c1) +
  facet_wrap(vars(cellType), scales = "free_y",nrow = 2L) +
  ggtitle("title to be decided")->p1

ploting(tp=p1,tw=13,th=5.7,
        tn=paste0('fig4d_sankeyPlot'))