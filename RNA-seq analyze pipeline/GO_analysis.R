##设置工作环境
setwd("/Users/xuanyuanqiao/Desktop/XYQ/1_学业/major\ courses/基因组/homework2/homework2_data")
gene_data<-read.csv("raw_count.csv")
head(gene_data)
if(!requireNamespace("BiocManager",quietly=TRUE))
  + install.packages("BiocManager")
BiocManager::install("DESeq2")
##安装DESeq
library(DESeq2)
rownames(gene_data)<-gene_data[,1]
gene_count<-gene_data[,-1]
head(gene_count)
##构建表达矩阵
condition <- factor(c("trt","trt","untrt","untrt"),levels=c("trt","untrt"))
col_data<-data.frame(row.names = colnames(gene_count),condition)
dds<-DESeqDataSetFromMatrix(countData = gene_count,colData = col_data,design = ~condition)
##将所有样本表达量之和小于1的都过滤掉
dds_filter<-dds[rowSums(counts(dds))>1,]
##对基因表达数据进行差异分析
dds_out<-DESeq(dds_filter)
res<-results(dds_out)
summary(res)
##设定阈值筛选差异基因并保存结果
res <-res[order(res$padj),]
diff_gene <- subset(res,padj<0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.csv(diff_gene,file = "DEG_trt_untrt.csv")
##log格式转化
rld<-rlogTransformation(dds_out,blind=F)
#找到方差最大对前20个基因并做热图
topVarGene <- head(order(rowVars(assay(rld)),decreasing = TRUE),20)
mat <-assay(rld)[topVarGene,]
library(pheatmap)
pheatmap(mat)

##用Excel将diff_gene.csv提出gene_id并多加一列后保存为GO_data.csv
GO_data<-read.csv("GO_data.csv")
##作图前处理——提取symbol ID --> 转换为ENTREZID
DEG.gene_symbol = as.character(GO_data$gene_id)
DEG.entrez_id = mapIds(x = org.Hs.eg.db,keys = DEG.gene_symbol, keytype = "SYMBOL",column = "ENTREZID")
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(org.Hs.eg.db)
## BP（Biological process）层面上的富集分析
erich.go.BP = enrichGO(gene = DEG.entrez_id,OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pvalueCutoff = 0.5,qvalueCutoff = 0.5)
##对差异显著的基因进行GO biological process富集分析，并做出最显著对前5条通路的结果
dotplot(erich.go.BP,showCategory=5)
