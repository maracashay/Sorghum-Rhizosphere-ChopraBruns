esvs.line <- phyloseq_to_deseq2(ps, ~ Line)
#use ps.pruned if running ITS data (ps.pruned excludes two samples that had 0 or close to 0 sequences)
esvs.line.ds <- DESeq(esvs.line, test="Wald", fitType="parametric")
resultsNames(esvs.line.ds)
alpha=0.05

# Line N321 v 326 #
res.321.326<- results(esvs.line.ds, cooksCutoff = FALSE, contrast = c("Line", "N326", "N321"))
sigtab.321.326.esv = res.321.326[which(res.321.326$padj < alpha), ]
sigtab.321.326.esv.df = data.frame(sigtab.321.326.esv) # convers the output to a data frame
sigtab.321.326.esv.df = cbind(fam_name=row.names(sigtab.321.326.esv.df), sigtab.321.326.esv.df)
sigtab.df.1.2.esv = sigtab.321.326.esv.df[order(sigtab.321.326.esv.df$log2FoldChange),]
# orders the esvs by the log 2 fold change
sigtab.df.1.2.esv

#re-run above and change line numbers for all possible combinations

# for genus/family/order/class level #

ps.genera.perc<- taxa_level(ps.pruned, "Genus") #taxa_level from the microbiome seq package
diagdds.flav <- phyloseq_to_deseq2(ps.genera.perc, ~ Line)
diagdds.flav <- DESeq(diagdds.flav, test="Wald", fitType="parametric")
resultsNames(diagdds.flav)
alpha=0.05

# Line N321 v 326 #
res<- results(diagdds.flav, cooksCutoff = FALSE, contrast = c("Line", "N321", "N326"))
sigtab.1.2 = res[which(res$padj < alpha), ]
sigtab.df.1.2 = data.frame(sigtab.1.2) # convers the output to a data frame
sigtab.df.1.2 = cbind(fam_name=row.names(sigtab.df.1.2), sigtab.df.1.2)
sigtab.df.1.2 = sigtab.df.1.2[order(sigtab.df$log2FoldChange),]
# orders the esvs by the log 2 fold change
sigtab.df.1.2

