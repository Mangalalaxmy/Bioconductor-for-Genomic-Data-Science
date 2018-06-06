source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("ShortRead")
library(ShortRead) 
biocLite("yeastRNASeq")
library(yeastRNASeq)
fastqFilePath <-  system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")
fqFile <- FastqFile(fastqFilePath)
reads <- readFastq(fqFile)
reads_set <- sread(reads)
sum(DNAStringSet(reads_set,5,5) == "A") / length(reads_set)

qm <- as(quality(reads), "matrix")
mean(qm[,5:5]) # 28.93

biocLite("leeBamViews")
library(leeBamViews)
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")
# In this interval, how many reads are duplicated by position?
bamFile <- BamFile(bamFilePath)
seqinfo(bamFile)
aln <- scanBam(bamFile)
aln <- aln[[1]]  
names(aln)
lapply(aln, function(xx) xx[1])
unique(aln$rname)
gr <- GRanges(seqnames = "Scchr13", ranges = IRanges(start = 800000, end = 801000))
params <- ScanBamParam(which = gr, what = scanBamWhat())
aln <- scanBam(bamFile, param = params)
aln <- aln[[1]]  
aln$pos # 327 total 
duplicatedValues = unique(aln$pos[duplicated(aln$pos)]) # duplicated positions and their names
sum(aln$pos %in% duplicatedValues)

library(GenomicRanges)
bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)
gr <- GRanges(seqnames = "Scchr13", ranges = IRanges(start = 807762, end = 808068))
# Scchr13:807762-808068
## ----BamViews------------------------------------------------------------
bamView <- BamViews(bpaths)
bamRanges(bamView) <- gr
aln <- scanBam(bamView)
names(aln)
names(aln[[1]])
lens <- list()
for(i in 1:length(aln)) {
  lens[i] <- length(aln[[i]][[1]]$seq)
}

mean(unlist(lens))

biocLite("bit")
biocLite("RSQLite")
library(bit)
library(RSQLite)
biocLite("oligo")
library(oligo)
biocLite("pd.hugene.1.0.st.v1") # affy and nimblegen arrays GE and snp arrays
library(pd.hugene.1.0.st.v1)
library(GEOquery)
geoMat <- getGEO("GSE38792")
pD.all <- pData(geoMat[[1]])
getGEOSuppFiles("GSE38792") # remember, raw data is in supplementary
list.files("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")
list.files("GSE38792/CEL") # samples from control to samples with sleep apnea
celfiles <- list.files("GSE38792/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)
rawData  # pd.hugene.1.0.st.v1  human gene vs 1 based on random priming
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),"OSA", "Control")
pData(rawData)
normData <- rma(rawData)
expr <- exprs(normData)
mean(expr["8149273",1:8])

library(limma)
design <- model.matrix(~ normData$group)
fit <- lmFit(normData, design)
fit <- eBayes(fit)
topTable(fit)
abs(topTable(fit, n=1)$logFC) # 0.7126

topTable(fit, p.value = 0.05)

library(minfi)
require(minfiData)
data(RGsetEx)
p <- preprocessFunnorm(RGsetEx)
b <- getBeta(p)
is <- getIslandStatus(p)
pData(p)$status
norm <- b[,c(1,2,5)]
can <- b[,c(3,4,6)]
norm_os <- norm[is == "OpenSea",]
can_os <- can[is == "OpenSea",]
mean(norm_os) - mean(can_os)

library(AnnotationHub)
ah <- AnnotationHub()
qah_h1 <- query(ah, c("Caco2", "AWG"))
h <- qah_h1[["AH22442"]]
sum(countOverlaps(p,h))  
ah_s <- subset(ah, genome == "hg19")
ah_s <- subset(ah, dataprovider == "UCSC")
# write.csv(mcols(ah), "ah.csv")
# g <- ah[["AH5018"]] # assembly
cpg <- ah_s[["AH5086"]] # CpG islands
h_cpg <- subsetByOverlaps(cpg, h)
ov <- subsetByOverlaps(h_cpg, p)

biocLite("DESeq2")
library(DESeq2)
biocLite("zebrafishRNASeq")
library(zebrafishRNASeq)
data(zfGenes)
#exclude spike-in controls
tx <- zfGenes[grep("^ERCC", rownames(zfGenes), invert = T),]
counts_mat <- as.matrix(tx)
colData <- DataFrame(sampleID=colnames(tx), group=as.factor(c("control", "control", "control", "treatment", "treatment", "treatment")))
ddsMat <- DESeqDataSetFromMatrix(counts_mat, colData, design = ~ group)
ddsMat <- DESeq(ddsMat)
res <- results(ddsMat)
res <- res[order(res$padj),]
sigRes <- subset(res, padj <= 0.05)
dim(sigRes)