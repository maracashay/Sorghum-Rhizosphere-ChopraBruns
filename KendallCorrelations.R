ps.filtered <- filter_taxa(ps.pruned, function(x) sum(x > 1) > (0.75*length(x)), TRUE)
esv.filtered <- data.frame(otu_table(ps.filtered))
metadata.df <- data.frame(sample_data(ps.filtered))
metadata.df <- metadata.df[,-c(3:11)]
cor.kendall <- cor(esv.filtered, metadata.df, method = "kendall")
write.table(data.frame(cor.kendall), file = "ASV~ContinuousVariable-Correlations-Kendall.csv")


#### Assess Correlations by "Kendall" Genera ####
ps.gen <- taxa_level(ps.pruned, "Genus")
ps.gen.filtered <- filter_taxa(ps.gen, function(x) sum(x > 1) > (0.75*length(x)), TRUE)
gen.filtered <- data.frame(otu_table(ps.gen.filtered))
metadata.df <- data.frame(sample_data(ps.gen.filtered))
metadata.df <- metadata.df[,-c(3:11)]
cor.kendall <- cor(gen.filtered, metadata.df, method = "kendall")
write.table(data.frame(cor.kendall), file = "Gen~ContinuousVariable-Correlations-Kendall.csv")

#### Correlations with Non-producers ####
ps.nprod <- subset_samples(ps.pruned, Line == "N321" | Line == "N326")

ps.np.filtered <- filter_taxa(ps.nprod, function(x) sum(x > 1) > (0.75*length(x)), TRUE)
esv.np.filtered <- data.frame(otu_table(ps.np.filtered))
metadata.np.df <- data.frame(sample_data(ps.np.filtered))
metadata.np.df <- metadata.np.df[,-c(3:11)]
cor.kendall.np <- cor(esv.np.filtered, metadata.np.df, method = "kendall")
write.table(data.frame(cor.kendall.np), file = "ASV~ContinuousVariable-Correlations-Kendall-NonProducers1.csv")


#### Correlations with Producers ####
ps.prod <- subset_samples(ps.pruned, Line == "N334" | Line == "N336")

ps.gen.p <- taxa_level(ps.prod, "Genus")
ps.gen.p.filtered <- filter_taxa(ps.gen.p, function(x) sum(x > 1) > (0.75*length(x)), TRUE)
gen.p.filtered <- data.frame(otu_table(ps.gen.p.filtered))

cor.kendall <- cor(gen.p.filtered, metadata.p.df, method = "kendall")
write.table(data.frame(cor.kendall), file = "Gen~ContinuousVariable-Correlations-Kendall-Producer.csv")
