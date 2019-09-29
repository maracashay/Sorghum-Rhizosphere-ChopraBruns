esv.time <- phyloseq_to_deseq2(ps.pruned, ~ Frost)
esv.time <- DESeq(esv.time, test="Wald", fitType="parametric")
res.esv.time = results(esv.time, cooksCutoff = FALSE)
alpha=0.05
res.esv.time = res.esv.time[which(res.esv.time$padj < alpha), ]
res.esv.time.df = data.frame(res.esv.time)
res.esv.time.df = cbind(fam_name=row.names(res.esv.time.df), res.esv.time.df)
res.esv.time.df = res.esv.time.df[order(res.esv.time.df$log2FoldChange),]
res.esv.time.df
