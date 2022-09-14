# load the data
# load the expression data
edata<-read.table("M_PR_raw_count.tsv",sep = "",header = TRUE)

rownames(edata)<-edata$ID

# add the gene symbols and length
table_human<-read.table("table_human_index.csv",header = TRUE,sep = ",",row.names = 1)

id<-match(rownames(edata),table_human$ensembl_gene_id)

edata$Length<-table_human$transcript_length[id]

edata$symbol<-table_human$external_gene_name[id]


# DEG analysis use EdgeR package
# Put the data into a DGEList object
library(edgeR)

genelist<-rownames(edata)

y<-DGEList(counts=edata[,2:97],genes=genelist)

# add the gene length information
y$genes$Length<-edata$Length

RPKM<-rpkm(y)

RPKM<-as.data.frame(RPKM)

ID2<-match(rownames(RPKM),table_human$ensembl_gene_id)

RPKM$symbol<-table_human$external_gene_name[ID2]

write.csv(RPKM,file = "M_PR_RPKM_annotated.csv")


