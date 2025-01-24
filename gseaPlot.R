#code is succeeding from 'RCCpaper_visualization.R'
#"EC_MDSC_interaction.R" is another code for RCCpaper visualizaiton

dyn.load("/apps/lib-osver/geos/3.5.0/lib/libgeos_c.so.1")
library(Seurat)
library(sctransform)

library(magrittr)
library(dplyr)
library(data.table)
library(ggplot2)

readRDS(file = 'deg_cluster5and11AcrossTime_fig5.rds')->tt_deg
#a. use msgbr genesets-----
tt_deg->deg4fgsea

deg4fgsea[sapply(deg4fgsea,class)!='character'] %>% 
  lapply(function(i){
    setNames(object = i$avg_log2FC,nm = rownames(i))
  })->deg4fgsea
rm(tt_deg)

all_gene_sets = msigdbr::msigdbr(species = "Mus musculus")[1:200,]
all_gene_sets %<>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.table()

c("gs_name","entrez_gene")->colnames(all_gene_sets)


#b.fgsea----
#to show all with p<0.05 instead of padj<0.05
trace(DOSE:::GSEA_fgsea,edit = T)#comment out line44
lab1<-'msgbrGenesets_'
lab1<-'yapGenesets_'

lapply(names(deg4fgsea),
       function(i){
         # browser()
         clusterProfiler::GSEA(
           geneList = sort(deg4fgsea[[i]],
                           decreasing = T),
           minGSSize = 0,maxGSSize = 500,exponent = 1,verbose = T,by = 'fgsea',
           TERM2GENE = all_gene_sets,
           eps = 0)->t_gsea
         if (nrow(t_gsea)<1){
           return(paste0(i," has no significant term enriched"))
         }else{
           t_gsea@result %<>% arrange(desc(NES))
           saveRDS(t_gsea,file = paste0('gseaObj_',lab1,
                                        i,'.rds'))
           return(paste0(i,": gseaObj has been saved as rds"))
         }
       })->tlog


# if files are small, just combine them together
list.files(pattern = paste0('gseaObj_',lab1,
                            '.*','.rds'))->lf
saveRDS(object=lapply(lf,readRDS) %>% 
          setNames(nm=lf %>% 
                     gsub(pattern = paste0('\\.rds|^.*',lab1),replacement = '')),
        file = paste0('gseaObj_',lab1,
                      '_combined','.rds'))
rm(lf)
#check if saved RDS is good
readRDS(file = paste0('gseaObj_',lab1,
                      '_combined','.rds'))->tt1
#if above tt1 is good, remove individual gseaObj
file.remove(setdiff(list.files(pattern = paste0('gseaObj_',lab1,
                                                '.*','.rds')),
                    paste0('gseaObj_',lab1,
                           '_combined','.rds')))

names(tt1) %>% sapply(function(i){
  tt1[[i]]@result %>% write.csv(file = paste0('gseaResult_',i,'.csv'))
})

# ...b1 treeplot----
library('clusterProfiler')
library('enrichplot')
debugonce(enrichplot:::pairwise_termsim.enrichResult)
pairwise_termsim(gseaObj1$`5__D7minusCTR`,
                 showCategory=400)->tt1

# debugonce(enrichplot:::treeplot.enrichResult)
# debugonce(enrichplot:::group_tree)
# debugonce(enrichplot:::add_cladelab)
# trace(enrichplot:::add_cladelab,edit = T)

treeplot(tt1,
         showCategory = gseaXls1$`5__D7minusCTR`[[1]][['terms']],
         # color = "p.adjust",
         color = "pvalue",
         nWords = 4, 
         nCluster = gseaXls1$`5__D7minusCTR`[[1]][['nTerm']],
         cex_category = 1,
         label_format = function(x) {substr(x,start = 3,stop = 12)}, 
         fontsize = 4, 
         offset = 1, 
         offset_tiplab = 1, 
         hclust_method = "ward.D", 
         group_color = NULL, 
         extend = 0.3, 
         hilight = TRUE, 
         hexpand = 0.1, 
         align = "both")->p1


ggplot2::ggsave(file.path(#'gseaPlot_fromCustomizedGenesets',
  paste0('_CustomizedGenesets_','.jpg')),
  plot = p1,
  units = c("in", "cm", "mm", "px"),
  width = 17.15,
  height = 31.12,
  dpi = 300)

# ...d2 barlot----
library('clusterProfiler')
library('enrichplot')

lapply(names(gseaObj1),function(i){
  
  gseaObj1[[i]]@result %<>% 
    arrange(desc(NES))
  
  enrichplot:::barplot.enrichResult(
    gseaObj1[[i]],
    showCategory=50,
    color="p.adjust",
    x='NES')+
    Seurat::WhiteBackground()+
    Seurat::NoGrid()+
    cowplot::theme_cowplot(font_size = 12,
                           line_size = 1)+
    labs(caption = paste0(
      "(NES<0 for terms enriched in CTR,NES>0 for terms enriched in ",
      gsub(i,pattern = '^.*__|minus.*$',replacement = ''),
      ")"
    )
    )+
    theme(
      # plot.title = element_text(color = "red", size = 12, face = "bold"),
      # plot.subtitle = element_text(color = "blue"),
      plot.caption = element_text(face = "italic")
    )->p1
  # saveRDS(p1,file = 'fig5f_barplot1.rds')
  # readRDS(file = 'fig5f_barplot1.rds')->p1
  # readRDS(file = 'fig5f_barplot.rds')->p

  ggplot2::ggsave(file.path(#'gseaPlot_fromCustomizedGenesets',
    paste0('_CustomizedGenesets_','.jpg')),
    plot = p1,
    units = c("in", "cm", "mm", "px"),
    width = 17.15,
    height = 31.12,
    dpi = 300)
  p1
}) %>% setNames(names(gseaObj1))->bps

