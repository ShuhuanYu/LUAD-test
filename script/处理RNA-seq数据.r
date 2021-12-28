setwd("RNA-seq-counts")

###数据整合与预处理
dir.create("00000_extracted_counts")
dir_counts<-dir()[nchar(dir())==36]
#decompress
for(data in dir_counts){
    filename=list.files(data,pattern = ".*counts")
    R.utils::gunzip(paste0(data,"/",filename))
    file_extracted=list.files(data,pattern = ".*counts$")
    file.copy(paste0(data,"/",file_extracted),"00000_extracted_counts")
}
#合并
setwd("00000_extracted_counts")
library(dplyr)
file_extracted_counts<-list.files()
first_file<-read.table(file_extracted_counts[1],header=F,sep="\t",stringsAsFactors=F)
names(first_file)<-c("gene_id",substr(file_extracted_counts[1],1,nchar(file_extracted_counts[1])-13))
for(extracted_counts in file_extracted_counts[2:length(file_extracted_counts)]) {
    file_appended=read.table(extracted_counts,header = T,sep = "\t",stringsAsFactors=F)
    names(file_appended)=c("gene_id",substr(extracted_counts,1,nchar(extracted_counts)-13))
    first_file=dplyr::inner_join(first_file,file_appended,by="gene_id")
}
expr_counts<-first_file
#处理样本id
library(jsonlite)
metadata<-jsonlite::fromJSON(txt="../../files/metadata.cart.2021-12-20.json")
metadata_id<-dplyr::select(.data=metadata,c(file_name,associated_entities))
sub_TCGA_ID_func <- function(x) {
   x$entity_submitter_id
}
ID_conversion_table<-data.frame(file_names=substr(metadata_id$file_name,1,36),
                                TCGA_IDs=as.character(lapply(metadata_id$associated_entities,FUN=sub_TCGA_ID_func)),
                                stringsAsFactors = F)
rownames(ID_conversion_table)<-ID_conversion_table$file_names
colnames(expr_counts)<-append("gene_id",ID_conversion_table[colnames(expr_counts)[2:length(expr_counts)],"TCGA_IDs"])
#更换ENSEMBL id为gene symbol
library(tidyr)
library(rtracklayer)
expr_counts<-tidyr::separate(expr_counts,gene_id,into=c("gene_id","junk"),sep="\\.")
expr_counts<-expr_counts[,-2]
gtf<-rtracklayer::import("../../files/Homo_sapiens.GRCh38.105.chr.gtf")
gtf_df<-as.data.frame(gtf)
#select protein-coding genes
mRNA_exprSet<- gtf_df %>%
    dplyr::filter(type=="gene",gene_biotype=="protein_coding") %>%
    dplyr::select(c(gene_name,gene_id,gene_biotype)) %>%
    dplyr::inner_join(expr_counts,by="gene_id")
#重复探针取均值
mRNA_exprSet1<-aggregate(mRNA_exprSet[,-c(1:3)],by=list(mRNA_exprSet$gene_name),mean)
colnames(mRNA_exprSet1)[1]<-"gene_name"
#select lncRNA genes
lncRNA_exprSet<- gtf_df %>%
    dplyr::filter(type=="gene",gene_biotype=="lncRNA") %>%
    dplyr::select(c(gene_name,gene_id,gene_biotype)) %>%
    dplyr::inner_join(expr_counts,by="gene_id")
lncRNA_exprSet1<-aggregate(lncRNA_exprSet[,-c(1:3)],by=list(lncRNA_exprSet$gene_name),mean)
colnames(lncRNA_exprSet1)[1]<-"gene_name"
#select miRNA genes
miRNA_exprSet<- gtf_df %>%
    dplyr::filter(type=="gene",gene_biotype=="miRNA") %>%
    dplyr::select(c(gene_name,gene_id,gene_biotype)) %>%
    dplyr::inner_join(expr_counts,by="gene_id")
miRNA_exprSet1<-aggregate(miRNA_exprSet[,-c(1:3)],by=list(miRNA_exprSet$gene_name),mean)
colnames(miRNA_exprSet1)[1]<-"gene_name"
#select snRNA genes
snRNA_exprSet<- gtf_df %>%
    dplyr::filter(type=="gene",gene_biotype=="snRNA") %>%
    dplyr::select(c(gene_name,gene_id,gene_biotype)) %>%
    dplyr::inner_join(expr_counts,by="gene_id")
snRNA_exprSet1<-aggregate(snRNA_exprSet[,-c(1:3)],by=list(snRNA_exprSet$gene_name),mean)
colnames(snRNA_exprSet1)[1]<-"gene_name"

#写出
write.table(mRNA_exprSet1,"../../files/mRNA_exprSet.txt",row.names = F,sep = "\t",quote = F)
write.table(lncRNA_exprSet1,"../../files/lncRNA_exprSet.txt",row.names = F,sep = "\t",quote = F)
write.table(miRNA_exprSet1,"../../files/miRNA_exprSet.txt",row.names = F,sep = "\t",quote = F)
write.table(snRNA_exprSet1,"../../files/snRNA_exprSet.txt",row.names = F,sep = "\t",quote = F)
save.image("../../script/数据清洗.Rdata")
