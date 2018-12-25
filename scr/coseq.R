# Title     : coseq
# Objective : 基因共表达
# Created by: yuelong.chen
# Created on: 2018/3/28
args<-commandArgs(T)

if(length(args) != 3){
  print("Rscript <this>.R <input.rda> <outdir> <alpha>")
  q()
}

library('coseq')
library('pheatmap')

infile=args[1]
outdir=args[2]
alpha=args[3]

print(paste(infile,outdir,alpha,sep='**'))
##load的是deseq导出的数据，会产生countData和condition，直接使用就可以
data=load(infile)
# print(data)
# print(countData)
#原始count进行所有基因共表达计算
#v-------------------------------------------------------------------------
runCoseq <- coseq(countData, K=2:40, transformation="arcsin",norm="TMM",
                   model="Normal",verbose=FALSE,meanFilterCutoff=200)

# 图像输出
outpdf = paste(outdir,'AllGeneCoExpress.pdf',sep='/')
pdf(outpdf)
plot(runCoseq,condition=condition)
dev.off()


# cluster输出，但是个人认为这个类别没有特别的作用

clusterResult = clusters(runCoseq)
outcsv = paste(outdir,'AllGeneCluster.csv',sep='/')
write.csv(file=outcsv,clusterResult,quote=FALSE)


#M-------------------------------------------------------------------------


#DESeq2计算后差异表达基因计算共表达
#v------------------------------------------------------------------------
runDeseq <- coseq(dds,K=2:40,verbose=FALSE,alpha=alpha)
outpdf = paste(outdir,'DESeqGeneCoExpress.pdf',sep='/')
pdf(outpdf)
plot(runDeseq,condition=condition)
dev.off()

clusterResult = clusters(runDeseq)
outcsv = paste(outdir,'DESeqGeneCluster.csv',sep='/')
write.csv(file=outcsv,clusterResult,quote=FALSE)

heatmapPdf = paste(outdir,'DESeqCoExpressHeatMap.pdf',sep='/')

clu = sort(clusterResult)
anno_row=data.frame(GeneClass=factor(clu))
cluData = tcounts(runDeseq)[names(clu),]
bk = unique(c(seq(-2,2, length=100)))
pheatmap(cluData,
        cluster_cols=FALSE,
        breaks=bk,
        annotation_row =anno_row,
        cluster_rows = FALSE,
        cellwidth = 22,
        fontsize_row = 0.5,
        scale = 'row',
        main='DEGs Coexpress Clusters',
        filename=heatmapPdf)



#M--------------------------------------------------------------------------
print('########finish#############')

