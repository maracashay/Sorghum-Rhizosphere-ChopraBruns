library(dada2)
library(ShortRead)
library(Biostrings)
library(devtools)
library(phyloseq)
library(DESeq2)
library(biomformat)
library(microbiome)
library(microbiomeSeq)
library(lme4)
library(lmerTest)
library(Hmisc)
library(ggplot2)
library(tibble)
library(vegan)
library(igraph)
library(visNetwork)
library(scales)
library(adespatial)
library(gtools)

#### 1. set working directory to set 1 ####
setwd("~/Downloads/Surinder_16S_Aug2018_raw_data")
path<- setwd("~/Downloads/Surinder_16S_Aug2018_raw_data")

list.files(path)

#### 2. Sort the forward and reverse reads #####
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

#### 3. Examine the quality profiles ####


quartz()
plotQualityProfile(fnFs[1:2])

quartz()
plotQualityProfile(fnRs[1:2])

#### 4. Assign Filtering samples ####

# Place filtered files in filtered/ subdirectory
filt_path <- file.path(path, "filtered") #you should see this in your global environment
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out.1 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 


head(out) 

#### 5. Train Dada2 to assess errors ####

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

quartz()
plotErrors(errF, nominalQ=TRUE)


#### 6. Dereplicate the filtered fastq files ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### 7. Infer the sequence variants in each sample #####

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)


dadaFs[[1]]
dadaRs[[1]]

#### 8.  Merge the denoised forward and reverse reads: ####

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)

dim(seqtab)


table(nchar(getSequences(seqtab)))

seqtab.1 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)] #this removes sequences shorter than 250 or longer than 256

saveRDS(seqtab.1, "~/Downloads/Surinders_ITS-seqs/SetII/ITS-Sample-Set-II/Pre-filtered/seqtab1.rds")

# clear workspace and rerun pipeline with set II

# then read in sequence table from first set
seqtab.1 <- readRDS("~/Downloads/Surinders_ITS-seqs/SetII/ITS-Sample-Set-II/Pre-filtered/seqtab1.rds")

# merge sequence tables
seqtab.full <- mergeSequenceTables(seqtab.1, seqtab.2)

#### 9.  Remove chimeric sequences ####



seqtab.nochim <- removeBimeraDenovo(seqtab.full, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab.full)


####10. Make phylogenetic tree ####
#Now, we can make a phylogenetic tree####
#based on paper by Callahan et al Bioconductor workflow for microbiome data analysis: from raw reads to community analyses

seqs<- getSequences(seqtab.nochim)

names(seqs) <- seqs

# BiocManager::install("msa")
mult <- msa(seqs, method = "ClustalW", type = "dna", class = "input")

#library(phangorn)
phang.aling <- as.phyDat(mult, type = "DNA", names = getSequences(seqtab.nochim))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data = phang.align)

#Phangorn package used to make tree


fitGTR <- update(fit, k=4, inv=0.2) #fit to a generalized time-reversible with gamma rate variation maximum liklihood tree 
fitGTR <- optim.pml(fitGTR, model= "GTR", optInv=TRUE, optGamma = TRUE, rearrangement = "stochastic", control = pml.control(trace=0))
detach("package:phangorn", unload=TRUE)

#### 11. Assign taxonomy ####
# taxonomic assignment with silva database
taxa.full<- assignTaxonomy(seqtab.nochim, "C:/Users/muc345/Downloads/silva_nr_v132_train_set.fa", multithread=TRUE)
#you should see assignments at taxonomic levels

taxa.full.sp <- assignSpecies(taxa.full, "C:/Users/muc345/Downloads/silva_species_assignment_v132.fa")


#### 12. Format data for analysis ####

Map<-import_qiime_sample_data("16S-metadata.txt")


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxa.full.sp), Map, phy_tree(fitGTR$tree))


# Now let's exchange the sequences for ASvs
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


##### Remove non bacteria/archaeal taxa ####
ps = subset_taxa(ps, Kingdom %in% c("Archaea", "Bacteria"))
