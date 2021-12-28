setwd("files")

###数据拆分以及标准化
# library(edgeR)
sample_info<-read.table("gdc_sample_sheet.2021-12-20.tsv",header = T,stringsAsFactors=F,sep="\t")
mRNA_exprSet<-read.table("mRNA_exprSet.txt",header=T,sep="\t",stringsAsFactors=F)
rownames(mRNA_exprSet)<-mRNA_exprSet$gene_name
mRNA_exprSet<-mRNA_exprSet[,-1]
colnames(mRNA_exprSet)<-substr(colnames(mRNA_exprSet),1,16)
colnames(mRNA_exprSet)<-gsub("\\.","-",colnames(mRNA_exprSet))
mRNA_exprSet<-mRNA_exprSet[,sample_info$Sample.ID]
#癌症样本
cancer_exprSet<-mRNA_exprSet[,which(sample_info$Sample.Type=="Primary Tumor")]
#标准化，TPM
#计算基因转录本长度——外显子
library(tidyr)
library(rtracklayer)
gtf<-rtracklayer::import("Homo_sapiens.GRCh38.105.chr.gtf")
gtf_df<-as.data.frame(gtf)
library(data.table)
library(IRanges)
anno<-setDT(gtf_df)
anno <- anno[type=="exon",]
setnames(anno,c("seqnames","start","end","gene_name","exon_number"),c("Chr","ExonStart","ExonEnd","Gene","Exon_number"))
#mkdir bin and mean by bin
Exon_region <- unique(anno[,.(Chr,ExonStart,ExonEnd,Exon_number,Gene)])
Exon_region <- Exon_region[,{x <- IRanges(ExonStart,ExonEnd);y <- reduce(x); list(ExonStart=y@start,ExonEnd=y@start+y@width-1)},by=.(Gene,Chr)]
Exon_region[,Exon_num:=1:.N,by=Gene]
Exon_region <- Exon_region[,.(Chr,ExonStart,ExonEnd,Exon_num,Gene)]
Exon_len <- Exon_region[,.(ExonLen = ExonEnd - ExonStart + 1),by=.(Exon_num,Gene)]
gene_len <- Exon_len[,.(Length = sum(ExonLen)),by=Gene]
fwrite(Exon_region,file="All_GRCh38.105gene_exon.bed", sep = "\t", col.names = T)
fwrite(gene_len, file = "All_GRCh38.105gene_len.txt", sep = "\t", col.names = T)
#计算TPM
gene_len<-gene_len[-which(is.na(gene_len$Gene)=="TRUE"),]
transcript_len<-t(gene_len)
colnames(transcript_len)<-transcript_len[1,]
transcript_len<-transcript_len[-1,]
transcript_len<-data.frame(length=transcript_len[rownames(cancer_exprSet)])
transcript_len$length<-as.numeric(transcript_len$length)
expr1<-cancer_exprSet/transcript_len$length
TPM<-t(t(expr1)/colSums(expr1))*10^6
write.table(TPM,"mRNA_TPMnormalization.txt",row.names = T,sep = "\t",quote = F)

###对癌症样本进行分型
#一致性聚类是一种为确定数据集中可能的聚类的数量和成员提供定量证据的方法，例如作为微阵列基因表达。这种方法在癌症基因组学中得到了广泛应用，在那里发现了新的疾病分子亚类
# BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
d<-TPM
mads<-apply(d,1,mad)
d<-d[rev(order(mads))[1:5000],]
d<-sweep(d,1, apply(d,1,median,na.rm=T))
title=getwd()
results<-ConsensusClusterPlus(d,maxK=20,reps=1000,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
#consensusMatrix - the consensus matrix.  
#For .example, the top five rows and columns of results for k=2:
results[[2]][["consensusMatrix"]][1:5,1:5]
#consensusTree - hclust object 
results[[2]][["consensusTree"]]
#consensusClass - the sample classifications
results[[2]][["consensusClass"]][1:5]
#ml - consensus matrix result
#clrs - colors for cluster  
icl = calcICL(results,title=title,plot="png")
icl[["clusterConsensus"]]
icl[["itemConsensus"]][1:5,]
#选择将样本划分成两类
consensusClass<-results[[2]][["consensusClass"]]
cluster1_expr<-TPM[,which(consensusClass=="1")] #第一类310个样本
cluster2_expr<-TPM[,which(consensusClass=="2")] #第二类223个样本

###节点过滤
#保留与生物学特征有关联的基因，这里保留在亚型间，以及癌症与正常样本间的DEGs
control_counts<-mRNA_exprSet[,which(sample_info$Sample.Type=="Solid Tissue Normal")]
cluster1_counts<-cancer_exprSet[,which(consensusClass=="1")]
cluster2_counts<-cancer_exprSet[,which(consensusClass=="2")]
library(DESeq2)
#cluster1 vs control
cts<-cbind(control_counts,cluster1_counts)
colData<-data.frame(condition=factor(rep(c("Solid_Normal_Tissue","Primary_Tumor"),c(ncol(control_counts),ncol(cluster1_counts)))))
rownames(colData)<-colnames(cts)
colData$condition<-relevel(colData$condition,ref="Solid_Normal_Tissue")
dds<-DESeqDataSetFromMatrix(countData = round(cts),
                              colData = colData,
                              design = ~ condition)
dds<-DESeq(dds)
res<-results(dds)
library(apeglm)
resLFC <- lfcShrink(dds, coef="condition_Primary_Tumor_vs_Solid_Normal_Tissue", type="apeglm")
# sum(resLFC$padj<0.01,na.rm=TRUE)
# library("IHW")
# resIHW <- results(dds, filterFun=ihw)   #主要在过滤p值时使用，与Benjamini和Hochberg的方法相比，它可以为每个假设分配数据驱动的权重来进行多重测试，从而增强统计效果。IHW的输入是p值和协变量的两列表，其中协变量可以是任何连续值或者分类变量，它可以为每个假设检验的统计特性提供信息，而与零假设下的p值无关
# summary(resIHW)
# plotMA(resLFC, ylim=c(-5,5))
#绘制火山图
library('EnhancedVolcano')  #The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
library(ggrepel)
pdf("(DEGs)Solid Normal Tissue vs Class1 LUAD Tumor.pdf")
EnhancedVolcano(resLFC,
lab = rownames(resLFC),
x = 'log2FoldChange',
y = 'pvalue',
selectLab = rownames(cts)[grep("WNT",rownames(cts))],
xlab = bquote(~Log[2]~ 'fold change'),
pCutoff = 10e-2,
FCcutoff = 2.0,
pointSize = 3.0,
labSize = 4.0,
labCol = 'black',
labFace = 'bold',
boxedLabels = TRUE,
colAlpha = 4/5,
# legendPosition = 'up',
legendLabSize = 14,
legendIconSize = 4.0,
# drawConnectors = TRUE,
# widthConnectors = 1.0,
# colConnectors = 'black',
title="Solid Normal Tissue vs Class1 LUAD Tumor"
) + coord_flip()
dev.off()
node_property<-data.frame(rownames(cts),cluster1_Control=apply(resLFC,1,function(x){
    ifelse(x[4]<=10e-2,ifelse(x[2]>=2,"Up",ifelse(x[2]<=-2,"Down","NS")),"NS")
}))
#cluster1 versus cluster2
cts1<-cbind(cluster2_counts,cluster1_counts)
colData1<-data.frame(condition=factor(rep(c("Cluster2_LUAD_Tumor","Cluster1_LUAD_Tumor"),c(ncol(cluster2_counts),ncol(cluster1_counts)))))
rownames(colData1)<-colnames(cts1)
colData1$condition<-relevel(colData1$condition,ref="Cluster1_LUAD_Tumor")
dds1<-DESeqDataSetFromMatrix(countData = round(cts1),
                              colData = colData1,
                              design = ~ condition)
dds1<-DESeq(dds1)
res1<-results(dds1)
library(apeglm)
resLFC1 <- lfcShrink(dds1, coef="condition_Cluster2_LUAD_Tumor_vs_Cluster1_LUAD_Tumor", type="apeglm")
# sum(resLFC$padj<0.01,na.rm=TRUE)
# library("IHW")
# resIHW <- results(dds, filterFun=ihw)   #主要在过滤p值时使用，与Benjamini和Hochberg的方法相比，它可以为每个假设分配数据驱动的权重来进行多重测试，从而增强统计效果。IHW的输入是p值和协变量的两列表，其中协变量可以是任何连续值或者分类变量，它可以为每个假设检验的统计特性提供信息，而与零假设下的p值无关
# summary(resIHW)
# plotMA(resLFC, ylim=c(-5,5))
#绘制火山图
library('EnhancedVolcano')  #The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
library(ggrepel)
pdf("(DEGs)Class1 LUAD Tumor vs Class2 LUAD Tumor.pdf")
EnhancedVolcano(resLFC1,
lab = rownames(resLFC1),
x = 'log2FoldChange',
y = 'pvalue',
selectLab = rownames(cts)[grep("WNT",rownames(cts))],
xlab = bquote(~Log[2]~ 'fold change'),
pCutoff = 10e-2,
FCcutoff = 2.0,
pointSize = 3.0,
labSize = 4.0,
labCol = 'black',
labFace = 'bold',
boxedLabels = TRUE,
colAlpha = 4/5,
# legendPosition = 'up',
legendLabSize = 14,
legendIconSize = 4.0,
# drawConnectors = TRUE,
# widthConnectors = 1.0,
# colConnectors = 'black',
title="Class1 LUAD Tumor vs Class2 LUAD Tumor"
) + coord_flip()
dev.off()
node_property<-data.frame(node_property,cluster2_Cluster1=apply(resLFC1,1,function(x){
    ifelse(x[4]<=10e-2,ifelse(x[2]>=2,"Up",ifelse(x[2]<=-2,"Down","NS")),"NS")
}))
#cluster2 versus control
cts2<-cbind(control_counts,cluster2_counts)
colData2<-data.frame(condition=factor(rep(c("Solid_Normal_Tissue","Primary_Tumor"),c(ncol(control_counts),ncol(cluster2_counts)))))
rownames(colData2)<-colnames(cts2)
colData2$condition<-relevel(colData2$condition,ref="Solid_Normal_Tissue")
dds2<-DESeqDataSetFromMatrix(countData = round(cts2),
                              colData = colData2,
                              design = ~ condition)
dds2<-DESeq(dds2)
res2<-results(dds2)
library(apeglm)
resLFC2 <- lfcShrink(dds2, coef="condition_Primary_Tumor_vs_Solid_Normal_Tissue", type="apeglm")
# sum(resLFC$padj<0.01,na.rm=TRUE)
# library("IHW")
# resIHW <- results(dds, filterFun=ihw)   #主要在过滤p值时使用，与Benjamini和Hochberg的方法相比，它可以为每个假设分配数据驱动的权重来进行多重测试，从而增强统计效果。IHW的输入是p值和协变量的两列表，其中协变量可以是任何连续值或者分类变量，它可以为每个假设检验的统计特性提供信息，而与零假设下的p值无关
# summary(resIHW)
# plotMA(resLFC, ylim=c(-5,5))
#绘制火山图
library('EnhancedVolcano')  #The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
library(ggrepel)
pdf("(DEGs)Solid Normal Tissue vs Class2 LUAD Tumor.pdf")
EnhancedVolcano(resLFC2,
lab = rownames(resLFC2),
x = 'log2FoldChange',
y = 'pvalue',
selectLab = rownames(cts)[grep("WNT",rownames(cts))],
xlab = bquote(~Log[2]~ 'fold change'),
pCutoff = 10e-2,
FCcutoff = 2.0,
pointSize = 3.0,
labSize = 4.0,
labCol = 'black',
labFace = 'bold',
boxedLabels = TRUE,
colAlpha = 4/5,
# legendPosition = 'up',
legendLabSize = 14,
legendIconSize = 4.0,
# drawConnectors = TRUE,
# widthConnectors = 1.0,
# colConnectors = 'black',
title="Solid Normal Tissue vs Class2 LUAD Tumor"
) + coord_flip()
dev.off()
node_property<-data.frame(node_property,cluster2_Control=apply(resLFC2,1,function(x){
    ifelse(x[4]<=10e-2,ifelse(x[2]>=2,"Up",ifelse(x[2]<=-2,"Down","NS")),"NS")
}))


#把差异基因当作种子节点用随机游走传播




###蛋白活性计算
# for(i in 2:10){
#     den=density(TPM[i,])
#     lines(den$x,den$y)  #大部分基因的表达分布都呈现负二项分布
# }
#计算概率分布函数，最大似然估计

#计算活性阈值
active_threshold <- function(expr) {
   u<-mean(expr)
   var<-var(expr)
   F<-1/(1+var)
   active_threshold=u+3*sd(expr)*(1-F)
}
threshold<-data.frame(Ac.threshold=unlist(apply(TPM,1,active_threshold)))
active_protein<-data.frame(TPM[,1]>threshold[,1])
for(i in 2:dim(TPM)[2]){
    active_protein<-cbind(active_protein,TPM[,i]>threshold[,1])
}










save.image("../script/蛋白质活性表达过滤.Rdata")
