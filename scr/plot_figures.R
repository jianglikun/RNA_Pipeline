# Title     : plotFigures
# Objective : plotFigures
# Created by: yuelong.chen
# Created on: 2018/4/4


args<-commandArgs(T)

if(length(args) != 3){
  print("Rscript <this>.R <input.rda> <alpha> <outdir> ")
  q()
}

library('pheatmap')
library('DESeq2')
library("RColorBrewer")

load(args[1])
rld <- rlog(dds, blind=FALSE)
res <- results(dds)
output = paste(args[3],'pheatmap.pdf',sep='/')
newres = na.omit(res)[na.omit(res)$padj<args[2],]
pheatmap(assay(rld)[rownames(newres),],
        show_rownames=TRUE,
        cellwidth = 22,
        fontsize_row = 0.5,
        cluster_cols=TRUE,
        cluster_rows=TRUE,
        annotation_col=colData,
        scale = 'row',
        filename=output,
        main='Heatmap of DEGs
        ')


