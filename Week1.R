source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("GenomicRanges")
biocLite("IRanges")
biocLite("AnnotationHub")

library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(AnnotationHub)
library(AnnotationDbi)

#Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.
#Question 1: How many islands exists on the autosomes?
ah = AnnotationHub()
ah = subset(ah, species == "Homo sapiens")
ah_cpg = query(ah, "CpG Islands")
ah_cpg
ah_cpg$genome
cpgdata = ah_cpg[[1]]
head(cpgdata)
cpgdata = dropSeqlevels(cpgdata, "chrX", pruning.mode = "coarse")
cpgdata = dropSeqlevels(cpgdata, "chrY", pruning.mode = "coarse")
cpgdata = keepStandardChromosomes(cpgdata, pruning.mode = "coarse")
cpgdata

#Question 2: How many CpG Islands exists on chromosome 4.
cpgdata4 = cpgdata
seqlevels(cpgdata4, pruning.mode = "coarse") = "chr4"
cpgdata4

#Obtain the data for the H3K4me3 histone modification for the H1 cell line from Epigenomics Roadmap, using AnnotationHub. Subset these regions to only keep regions mapped to the autosomes (chromosomes 1 to 22).
#Question 3: How many bases does these regions cover?
ah_k4 = query(ah, c("H3K4me3", "E003"))
ah_k4
ah_k4[2]
dat1 = ah_k4[["AH29884"]]
head(dat1)
dat1 = dropSeqlevels(dat1, "chrX", pruning.mode = "coarse")
dat1 = dropSeqlevels(dat1, "chrY", pruning.mode = "coarse")
dat1 = keepStandardChromosomes(dat1, pruning.mode = "coarse")
dat1
sum(width(dat1))

#Obtain the data for the H3K27me3 histone modification for the H1 cell line from Epigenomics Roadmap, using the AnnotationHub package. Subset these regions to only keep regions mapped to the autosomes. In the return data, each region has an associated "signalValue".
#Question 4: What is the mean signalValue across all regions on the standard chromosomes?
ah_k27 = query(ah, c("H3K27me3", "E003"))
ah_k27
ah_k27[2]
dat2 = ah_k27[["AH29892"]]
head(dat2)
dat2 = dropSeqlevels(dat2, "chrX", pruning.mode = "coarse")
dat2 = dropSeqlevels(dat2, "chrY", pruning.mode = "coarse")
dat2 = keepStandardChromosomes(dat2, pruning.mode = "coarse")
dat2
mean(dat2$signalValue)

#Bivalent regions are bound by both H3K4me3 and H3K27me3.
#Question 5: Using the regions we have obtained above, how many bases on the standard chromosomes are bivalently marked?
bi_valent = intersect(dat1, dat2)
sum(width(bi_valent))

#We will examine the extent to which bivalent regions overlap CpG Islands.
#Question 6: how big a fraction (expressed as a number between 0 and 1) of the bivalent regions, overlap one or more CpG Islands?
overlaps_bi_cpg <- findOverlaps(bi_valent, cpgdata)
length(unique(queryHits(overlaps_bi_cpg)))/length(bi_valent)

#Question 7: How big a fraction (expressed as a number between 0 and 1) of the bases which are part of CpG Islands, are also bivalent marked.
inter_bi_cpg = intersect(bi_valent, cpgdata)
sum(width(inter_bi_cpg))/sum(width(cpgdata))

#Question 8: How many bases are bivalently marked within 10kb of CpG Islands?
#Tip: consider using the "resize()"" function.
resized_cpgdata = resize(cpgdata, width=20000+width(cpgdata), fix="center")
inter_bi_resizedcpgdata = intersect(bi_valent, resized_cpgdata)
sum(width(inter_bi_resizedcpgdata))

#Question 9: How big a fraction (expressed as a number between 0 and 1) of the human genome is contained in a CpG Island?
#Tip 1: the object returned by AnnotationHub contains "seqlengths".
#Tip 2: you may encounter an integer overflow. As described in the session on R Basic Types, you can address this by converting integers to numeric before summing them, "as.numeric()".
ah_genome = query(ah, c("hg19","Assembly"))
genome = ah_genome[[1]]
genome = dropSeqlevels(genome, "chrX", pruning.mode = "coarse")
genome = dropSeqlevels(genome, "chrY", pruning.mode = "coarse")
genome = keepStandardChromosomes(genome, pruning.mode = "coarse")
genome_size = sum(as.numeric(seqlengths(genome)))
sum(as.numeric(width(cpgdata)))/genome_size

#Question 10: Compute an odds-ratio for the overlap of bivalent marks with CpG islands.
inOut = matrix(0, ncol = 2, nrow = 2)
colnames(inOut) = c("in", "out")
rownames(inOut) = c("in", "out")
inOut
inOut[1,1] = sum(width(intersect(bi_valent, cpgdata)))
inOut[1,2] = sum(width(setdiff(bi_valent, cpgdata)))
inOut[2,1] = sum(width(setdiff(cpgdata, bi_valent)))
inOut[2,2] = genome_size - sum(inOut)
inOut
odd_ratio <- inOut[1,1]*inOut[2,2]/(inOut[1,2]*inOut[2,1])
odd_ratio
