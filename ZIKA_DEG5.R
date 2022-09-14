# load the data
#load the phenotype table
pdata=read.table("phenotype_data_ZIKA5.txt",sep = "",row.names = 1)
# load the expression data
edata=read.table("expression_data_ZIKA5.txt",sep = "",row.names = 1)

# table for factor/character variables
pdata$group=as.factor(pdata$group)

table(pdata$group)


# DEG analysis use EdgeR package
# Put the data into a DGEList object
library(edgeR)

genelist<-rownames(edata)

y<-DGEList(counts=edata,genes=genelist)

# add transcript length
# load the match table which have the gene_name, gene_id,and entrez_gene_id
table_human<-read.table("table_human_index.csv",header = TRUE,sep = ",",row.names = 1)

id<-match(rownames(y$genes),table_human$ensembl_gene_id)

y$genes$Length<-table_human$transcript_length[id]

RPKM<-rpkm(y)

RPKM<-as.data.frame(RPKM)

id<-match(rownames(RPKM),table_human$ensembl_gene_id)

RPKM$symbol<-table_human$external_gene_name[id]

write.csv(RPKM,file = "compare5_RPKM.csv")


# Filtering
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) > 1)
y <- y[keep, ]


# Normalization
y <- calcNormFactors(y, method="TMM")

y$samples$group <- pdata$group

group<-pdata$group

design <- model.matrix(~0+group)

colnames(design) <- levels(group)

design

# the DEA result for all the genes
# dea <- lrt$table

y <- estimateDisp(y, design, robust = TRUE)

fit<-glmQLFit(y,design,robust = TRUE)

MOCK2_vs_MOCK1<-makeContrasts(mock2-mock1,levels = design)

res1<-glmQLFTest(fit,contrast = MOCK2_vs_MOCK1)

toptag1 <- topTags(res1, n = nrow(y$genes), p.value = 1)

dea1 <- toptag1$table 

dea1 <- dea1[order(dea1$FDR, -abs(dea1$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

id1<-match(dea1$genes,table_human$ensembl_gene_id)

dea1$symbol<-table_human$external_gene_name[id1]

write.csv(dea1,file = "compare5_mock_DEG.csv")



IN2_vs_IN1<-makeContrasts(IN2-IN1,levels = design)

res2<-glmQLFTest(fit,contrast = IN2_vs_IN1)

toptag2 <- topTags(res2, n = nrow(y$genes), p.value = 1)

dea2 <- toptag2$table 

dea2 <- dea2[order(dea2$FDR, -abs(dea2$logFC), decreasing = FALSE), ]  # sort the table: ascending of FDR then descending of absolute valued of logFC

id2<-match(dea2$genes,table_human$ensembl_gene_id)

dea2$symbol<-table_human$external_gene_name[id2]

write.csv(dea2,file = "compare5_infect_DEG.csv")


