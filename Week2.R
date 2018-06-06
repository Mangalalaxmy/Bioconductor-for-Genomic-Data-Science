source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("Biostrings")
biocLite("BSgenome")
biocLite("GenomicRanges")
biocLite("IRanges")
biocLite("AnnotationHub")
biocLite("AnnotationDbi")
biocLite("GenomicFeatures")
biocLite("shiny")
biocLite("rtracklayer")
biocLite("stringi")
biocLite("bit")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(AnnotationHub)
library(AnnotationDbi)
library(Biostrings)
library(BSgenome)
library(GenomicFeatures)
library(shiny)
library(rtracklayer)
library(stringi)
library(bit)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Question 1: What is the GC content of "chr22" in the "hg19" build of the human genome?
#Tip: The reference genome includes "N" bases; you will need to exclude those.
biocLite("BSgenome.Hsapiens.UCSC.hg19")
available.genomes()
library("BSgenome.Hsapiens.UCSC.hg19")
A = letterFrequency(Hsapiens$chr22, "A")
t = letterFrequency(Hsapiens$chr22, "T")
g = letterFrequency(Hsapiens$chr22, "G")
c = letterFrequency(Hsapiens$chr22, "C")
Total = A + t + g + c
GC = letterFrequency(Hsapiens$chr22, "GC")
g + c == GC
GCcontent = GC/Total
GCcontent

#Background: In the previous assessment we studied H3K27me3 "narrowPeak" regions from the H1 cell line (recall that the Roadmap ID for this cell line is "E003"). We want to examine whether the GC content of the regions influence the signal; in other words wether the reported results appear biased by GC content.
#Question 2: What is mean GC content of H3K27me3 "narrowPeak" regions from Epigenomics Roadmap from the H1 stem cell line on chr 22.
#Clarification: Compute the GC content for each peak region as a percentage and then average those percentages to compute a number between 0 and 1.
ah = AnnotationHub()
ah = subset(ah, species == "Homo sapiens")
ah_k27 = query(ah, c("H3K27me3", "E003"))
ah_k27
ah_k27[2]
k27narrowdata = ah_k27[["AH29892"]]
head(k27narrowdata)
k27data_22 = keepSeqlevels(k27narrowdata, "chr22", pruning.mode = "coarse")
k27data_22_seq = Views(Hsapiens, k27data_22)
gcContents = letterFrequency(k27data_22_seq, "GC", as.prob = TRUE)
meanGcContents = mean(gcContents)
meanGcContents

#The "narrowPeak" regions includes information on a value they call "signalValue".
#Question 3: What is the correlation between GC content and "signalValue" of these regions (on chr22)?
sigval = mcols(k27data_22_seq)$signalValue
cor(sigval, gcContents)

#The "narrowPeak" regions are presumably reflective of a ChIP signal in these regions. To confirm this, we want to obtain the "fc.signal" data from AnnotationHub package on the same cell line and histone modification. This data represents a vector of fold-change enrichment of ChIP signal over input.
#Question 4: what is the correlation between the "signalValue" of the "narrowPeak" regions and the average "fc.signal" across the same regions?
#Clarification: First compute the average "fc.signal" for across each region, for example using "Views"; this yields a single number of each region. Next correlate these numbers with the "signalValue" of the "narrowPeaks".
ah_k27_fc = ah_k27[["AH32033"]]
ah_k27_fc
k27_fc_rle = import(ah_k27_fc, which=GRanges("chr22", ranges = IRanges(start =1, end = 10^8)), as="Rle")
k27_fc_rle22 = k27_fc_rle$chr22
fc_signal_22 = Views(k27_fc_rle22, start = start(k27data_22), end = end(k27data_22))
fc_signal_mean = mean(fc_signal_22)
cor(fc_signal_mean, sigval)

#Referring to the objects made and defined in the previous question.
#Question 5: How many bases on chr22 have an fc.signal greater than or equal to 1?
sum(k27_fc_rle22 >= 1)

#The H1 stem cell line is an embryonic stem cell line, a so-called pluripotent cell. Many epigenetic marks change upon differentiation. We will examine this. We choose the cell type with Roadmap ID "E055" which is foreskin fibroblast primary cells.
#We will use the "fc.signal" for this cell type for the H3K27me3 mark, on chr22. We now have a signal track for E003 and a signal track for E055. We want to identify regions of the genome which gain H3K27me3 upon differentiation. These are regions which have a higher signal in E055 than in E003. To do this properly, we would need to standardize (normalize) the signal across the two samples; we will ignore this for now.
#Question 6: Identify the regions of the genome where the signal in E003 is 0.5 or lower and the signal in E055 is 2 or higher.
#Tip: If you end up with having to intersect two different Views, note that you will need to convert the Views to IRanges or GRanges first with \verb|ir <- as(vi, "IRanges")|ir<-as(vi,"IRanges").
ah_k27_55 = query(ah, c("H3K27me3", "E055"))
ah_k27_55
ah_k27_55fc = ah_k27_55[["AH32470"]]
ah_k27_55fc
k27_55fc_rle = import(ah_k27_55fc, which=GRanges("chr22", ranges = IRanges(start =1, end = 10^8)), as="Rle")
k27_55fc_rle22 = k27_55fc_rle$chr22
reg_3 = slice(k27_fc_rle22, upper = 0.5)
reg_55 = slice(k27_55fc_rle22, lower = 2)
reg_3 = as(reg_3, "IRanges")
reg_55 = as(reg_55, "IRanges")
inter_region = intersect(reg_3, reg_55)
sum(width(inter_region))

#CpG Islands are dense clusters of CpGs. The classic definition of a CpG Island compares the observed to the expected frequencies of CpG dinucleotides as well as the GC content.
#Specifically, the observed CpG frequency is just the number of "CG" dinucleotides in a region. The expected CpG frequency is defined as the frequency of C multiplied by the frequency of G divided by the length of the region.
#Question 7: What is the average observed-to-expected ratio of CpG dinucleotides for CpG Islands on chromosome 22?
ah_cpg = query(ah, c("hg19", "CpG"))
ah_cpg
cpg_data = ah_cpg[["AH5086"]]
cpg_22 = keepSeqlevels(cpg_data, "chr22", pruning.mode = "coarse")
cpg_22views = Views(Hsapiens, cpg_22)
region_length = width(cpg_22views)
obs_freq = dinucleotideFrequency(cpg_22views)[,7]/region_length
freq_C = letterFrequency(cpg_22views, "C")
freq_G = letterFrequency(cpg_22views, "G")
exp_freq = (freq_C/region_length)*(freq_G/region_length)
mean(obs_freq/exp_freq)

#A TATA box is a DNA element of the form "TATAAA". Around 25% of genes should have a TATA box in their promoter. We will examine this statement.
#Question 8: How many TATA boxes are there on chr 22 of build hg19 of the human genome?
#Clarification: You need to remember to search both forward and reverse strands.
TATA_boxes = countPattern("TATAAA", Hsapiens$chr22) + countPattern("TATAAA", reverseComplement(Hsapiens$chr22))
TATA_boxes

#Question 9: How many promoters of transcripts on chromosome 22 containing a coding sequence, contains a TATA box on the same strand as the transcript?
#Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence. Here, we defined a promoter to be 900bp upstream and 100bp downstream of the transcription start site.
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
gr = GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 52330658))
gr_trans22 = subsetByOverlaps(transcripts(txdb), gr, ignore.strand = TRUE)
proms = promoters(gr_trans22, upstream = 900, downstream = 100)
cdseq = subsetByOverlaps(genes(txdb), gr, ignore.strand = TRUE)
proms_cds = findOverlaps(proms, cdseq)
unique(queryHits(proms_cds))
count = 0
for (i in unique(queryHits(proms_cds))){
  proms_cds_vi <- Views(Hsapiens, proms[i])
  count = count + vcountPattern("TATAAA", DNAStringSet(proms_cds_vi))
}
count

#It is possible for two promoters from different transcripts to overlap, in which case the regulatory features inside the overlap might affect both transcripts. This happens frequently in bacteria.
#Question 10: How many bases on chr22 are part of more than one promoter of a coding sequence?
#Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to define transcripts and coding sequence. Here, we define a promoter to be 900bp upstream and 100bp downstream of the transcription start site. In this case, ignore strand in the analysis.
gr = GRanges(seqnames = "chr22", ranges = IRanges(start = 1, end = 52330658))
gr_trans22 = subsetByOverlaps(transcripts(txdb), gr, ignore.strand = TRUE)
length(gr_trans22) 
proms = promoters(gr_trans22, upstream = 900, downstream = 100)
tl_chr22 = transcriptLengths(txdb, with.cds_len = TRUE) #rtn df
tl_chr22 = tl_chr22[tl_chr22$cds_len > 0,]
trans_eval = proms[mcols(proms)$tx_id %in% tl_chr22$tx_id]
sum(coverage(trans_eval) > 1)
