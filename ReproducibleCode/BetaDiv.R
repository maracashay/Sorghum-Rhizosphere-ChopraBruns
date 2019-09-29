# Bray curtis # 
bc_rt_highcover_ns_abun <- phyloseq::distance(otu_table(ps), "bray")
adonis(bc_rt_highcover_ns_abun ~ Line * Frost, data = meta.df)


# unifrac weighted
w.unifrac <- UniFrac(ps, weighted = TRUE, normalized = TRUE)
w.u.pcoa<- pcoa(w.unifrac)

# global permanova 
meta.df <- data.frame(sample_data(ps))
w.u.perm <- adonis(w.unifrac ~ Line * Frost, data = meta.df)
w.u.perm

# pairwise permanova
pair.perm.line <- pairwise.perm.manova(bc_rt_highcover_ns_abun, meta.df$Line ,nperm=500)
pair.perm.line
