# the ITS seqs prepared by the Wright lab had a large amount of seqs that had low complexity (many repetitive seqs)
# seqs had to be pre-filtered using fastP- using the following commands
# conda install -c bioconda fastp
# fastp -y low_complexity_filter -i RI-4_S4_R1_001.fastq -I RI-4_S4_R2_001.fastq -o RI4_R1_001.fastq -O RI4_R2_001.fastq
# Resulting in the removal of low complexity reads

source("https://bioconductor.org/biocLite.R")
#biocLite("dada2")
#biocLite("ShortRead")
#biocLite("Biostrings")
#biocLite("phyloseq")
#biocLite("microbiome")
#biocLite("Deseq2")
#biocLite("biomformat")

library(dada2)
library(ShortRead)
library(Biostrings)
library(phyloseq)
library(microbiome)
library(DESeq2)
library(biomformat)

#devtools::install_github("benjjneb/dada2", force=TRUE)

library(devtools)

#install_github("umerijaz/microbiomeSeq")  # Install the package
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

# Set Working Directory for Sequence Set I ####
setwd("D:/QualityFiltered")
path<- setwd("D:/QualityFiltered")

list.files(path)


# First we read in the names of the fastq files, and perform some string manipulation to get lists of the 
# forward and reverse fastq files in matched order


#### 2. Sort the forward and reverse reads #####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
#You should see in your global env, fnFs chr[1:20] "and the path directory"
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
#this is depedent on how you samples are named

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
#should see a list of the sample names


##### 3. Examine the quality profiles ####

## forward reads ##

# plot A in the Dada2 Outputs  # 
quartz()
#jpeg(file="ForwardRead_QualityProfile.jpeg", width=10, height=8, units="in", res=800)
plotQualityProfile(fnFs[1:4])
#dev.off()
## reverse reads ##
# plot B in the Dada2 Outputs  # 
quartz()
#jpeg(file="ReverseRead_QualityProfile.jpeg", width=10, height=8, units="in", res=800)
plotQualityProfile(fnRs[1:4])
#dev.off()

#Check complexities of the forward and reverse reads
cmp1 <- seqComplexity(getSequences(filtFs[[1]])); 
#jpeg(file="Forward_cmp1.jpeg", width=10, height=8, units="in", res=800)
plot(cmp1)
#dev.off()

cmpR <- seqComplexity(getSequences(filtRs[[1]]))
#jpeg(file="Reverse_cmp1.jpeg", width=10, height=8, units="in", res=800)
plot(cmpR)
#dev.off()

#### 4. Assign Filtering samples ####
# Assign the filenames for the filtered fastq.gz files

# Place filtered files in filtered/ subdirectory
filt_path <- file.path(path, "filtered") #you should see this in your global environment
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(120,120),
                               maxN=0, maxEE=c(3,3), truncQ=2,
                               compress=TRUE, multithread=FALSE) 


#### 5. Train Dada2 to assess errors ####
errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

## Visualize errors
plotErrors(errF, nominalQ=TRUE)

#### 6. Dereplicate the filtered fastq files ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### 7. Infer the sequence variants in each sample #####

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]


#### 8.  Merge the denoised forward and reverse reads: ####

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.1 <- seqtab[,nchar(colnames(seqtab)) %in% seq(187,203)] #this removes sequences shorter than 250 or longer than 256

#export table 
saveRDS(seqtab.1, "D:/Surinders_ITS_BothSets/seqtab-setI.rds")

# save workspace and then clear workspace
# rerun pipeline on second set of ITS sequences 
# read in sequence table from second set
seqtab.1 <- readRDS("~/Downloads/Surinders_ITS-seqs/SetI/seqtab-setI.rds")

# after running pipeline to here:
#export table 
saveRDS(seqtab2.setII, "~/Downloads/Surinders_ITS-seqs/SetII/ITS-Sample-Set-II/Pre-filtered/seqtab-setII.rds")

# now combine sequence tables from both sequencing runs
seqtab.full <- mergeSequenceTables(seqtab.setI, seqtab.setII)
#### 9.  Remove chimeric sequences ####

seqtab.nochim.full <- removeBimeraDenovo(seqtab.full, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim.full)
sum(seqtab.nochim.full)/sum(seqtab.full)

getN <- function(x) sum(getUniques(x))

track <- cbind(rowSums(seqtab.full), rowSums(seqtab.nochim.full))

colnames(track) <- c("tabled", "nonchim")
head(track)
track

#### 10. Assign taxonomy ####

# taxonomy assigned using the UNITE database
# go to the Unite website https://unite.ut.ee/repository.php and the "General FASTA Release" and download the most recent file 

unite.ref<- "C:/Users/muc345/Downloads/Unite/sh_general_release_dynamic_02.02.2019.fasta"

taxonomy <- assignTaxonomy(seqtab.nochim.full, unite.ref, multithread = TRUE, tryRC = TRUE)

# check your taxonomic classifications #
taxa.print<- taxa
rownames(taxa.print)
head(taxa.print)

#### 11. Format data for analysis ####

#now we need to import our mapping file or sometimes referred to as our metadata
Map<-import_qiime_sample_data("DirectorytoFungalMappingFile.txt")


esv.table<-otu_table(seqtab.nochim.full.mod, taxa_are_rows=FALSE)

#Now we can make the phyloseq object
ps <- phyloseq(esv.table, tax_table(taxonomy), Map)

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#### 12. general processing and analysis using silva phyloseq object####
##### First thing we do is count the number of seqs per sample and remove samples with low read counts
sample_sums(ps)
#there are two samples with 0 reads, let's remove them by setting a low sequence count threshold
ps.pruned <- prune_samples(sample_sums(ps)>=2000, ps)
sample_sums(ps.pruned)


ps.pruned = subset_taxa(ps.pruned, Kingdom =="k__Fungi")

# Since we removed two samples from the phyloseq object, we need to create a new mapping file
Map.1 <- data.frame(sample_data(ps.pruned))

#### Rarefy dataset for Alpha-diversity analyses ####
ordered(sampleSums(ps.pruned)) 
set.seed(500)
ps.rare<-rarefy_even_depth(ps, sample_size = 38000)
