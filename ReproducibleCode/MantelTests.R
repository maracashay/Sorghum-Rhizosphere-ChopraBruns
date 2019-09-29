xdist <- vegdist(asv.table, method = "bray")
metadata.df <- metadata.df[,-c(1:7)]
ydist <- vegdist(metadata.df, method = "euclid")

vegan::mantel(xdist, ydist, method = "spear")

lm.all <- lm(xdist~ydist)
jpeg(file="BC-Dist~AllDist.jpg", width=10, height=8, units="in", res=800)
plot(xdist~ydist, ylab = "Bray-Curtis Dissimilarity", xlab = "Variable Dissimilarity", pch = 19, frame = FALSE, cex = 1, font = 2)
abline(lm.all, col = "blue", lty =1, lwd=4 )
dev.off()

flav.dist <- vegdist(metadata.df$Total_Flavonoids, method = "euclid")
vegan::mantel(xdist, flav.dist, method = "spear")

lm.flav <- lm(xdist~flav.dist)
jpeg(file="BC-Dist~FlavDist-BothTimePoints.jpg", width=10, height=8, units="in", res=800)
plot(xdist~flav.dist, ylab = "Bray-Curtis Dissimilarity", xlab = "Flavonoid Dissimilarity", pch = 19, frame = FALSE, cex = 1)
abline(lm.flav, col = "blue", lty =1, lwd=4 )
dev.off()

# now with phenolics
phenol.dist <- vegdist(metadata.df$Total_Phenolics, method = "euclid")
vegan::mantel(xdist, phenol.dist, method = "spear")

# now with DPPH
DPPH.dist <- vegdist(metadata.df$DPPH, method = "euclid")
vegan::mantel(xdist, DPPH.dist, method = "spear")

# now with Api
Api.dist <- vegdist(metadata.df$X3DA_Luteo, method = "euclid")
vegan::mantel(xdist, Api.dist, method = "spear")

# to do these with producer and non-producer lines separately
# subset phyloseq object accordingly
ps.prod <- subset_samples(ps.pruned, Line == "N334" | Line == "N336")
ps.nprod <- subset_samples(ps.pruned, Line == "N321" | Line == "N326")
