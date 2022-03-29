install.packages("ggrepel")
BiocManager::install("biomaRT")
library(ggplot2)
library(ggrepel)
library('biomaRt')
Deseq_Res = read.csv(choose.files())
Deseq_Res = na.omit(Deseq_Res)
Deseq_Res$threshold = as.numeric(Deseq_Res$pvalue < 0.01)
Deseq_Res$threshold2 = as.numeric(Deseq_Res$log2FoldChange < -1)
Deseq_Res$threshold3 = as.numeric(Deseq_Res$log2FoldChange > 1)
Deseq_Res$threshold4 = Deseq_Res$threshold2 + Deseq_Res$threshold + Deseq_Res$threshold3
Deseq_Res2 = Deseq_Res[Deseq_Res$threshold4 == 2,]
Deseq_Res2 = Deseq_Res2[order(Deseq_Res2$pvalue),]
Deseq_Res2 = Deseq_Res2[1:20,]
genes <- Deseq_Res2$X
ensembl <- useMart("ensembl", host= "www.ensembl.org", dataset= "drerio_gene_ensembl")
genes_ensembl_org <- getBM(attributes <- c("entrezgene_id", "ensembl_gene_id", "external_gene_name", "description"), filters = "ensembl_gene_id", values = genes, mart = ensembl, uniqueRows=T)
pmatch_table <- pmatch(genes, genes_ensembl_org[,2], duplicates.ok=T)
ensembl_table <- as.data.frame(matrix(NA, nrow=length(genes), ncol=8))
ensembl_table[which(!is.na(pmatch_table)),] <- genes_ensembl_org[pmatch_table[(!is.na(pmatch_table))], ];
rownames(ensembl_table) <- genes;
colnames(ensembl_table) <- colnames(genes_ensembl_org);
results <- cbind(ensembl_table[,3:4], Deseq_Res2);
colnames(results) <- c("symbol", "description",  colnames(Deseq_Res2));
rownames(results) <- rownames(Deseq_Res2)
g = ggplot(data = Deseq_Res)
g + geom_point(aes(x = log2FoldChange, y = -log10(pvalue), col = as.factor(threshold4))) +
geom_hline(yintercept= -log10(0.01)) +
geom_vline(xintercept = c(1,-1)) +
scale_colour_manual(values = c("#666666","#0066CC","#CC0000")) +
geom_text_repel(data = results, 
                aes(x = log2FoldChange,
                        y = -log10(pvalue)),
                    label = results$symbol,
                    min.segment.length = 0,
                    size = 3,
                box.padding = 0.5,
                nudge_x = .4,
                nudge_y = .4,) +
  theme_light() + 
  theme(legend.position = "None")
