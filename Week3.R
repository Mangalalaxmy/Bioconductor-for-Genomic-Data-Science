source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("ALL")
library(ALL)

#Question 1: What is the mean expression across all features for sample 5 in the ALL dataset (from the ALL package)?
data(ALL)
ALL
mean(exprs(ALL)[,5])

#We will use the biomaRt package to annotate an Affymetrix microarray. We want our results in the hg19 build of the human genome and we therefore need to connect to Ensembl 75 which is the latest release on this genome version. How to connect to older versions of Ensembl is described in the biomaRt package vignette; it can be achived with the command \verb|mart <- useMart(host=&#x27;feb2014.archive.ensembl.org&#x27;, biomart = "ENSEMBL_MART_ENSEMBL")|mart<-useMart(host=&#x27;feb2014.archive.ensembl.org&#x27;,biomart="ENSEMBL_MART_ENSEMBL").
#Question 2: Using this version of Ensembl, annotate each feature of the ALL dataset with the Ensembl gene id. How many probesets (features) are annotated with more than one Ensembl gene id?
biocLite("biomaRt")
install.packages("stringi")
library(stringi)
library(biomaRt)
mart = useMart(host='feb2014.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL")
listDatasets(mart)
ensembl_data = useDataset("hsapiens_gene_ensembl", mart)
ids = featureNames(ALL)
attributePages(ensembl_data)
att = listAttributes(ensembl_data, page="feature_page")
result = getBM(attributes = c("affy_hg_u95av2", "ensembl_gene_id", "chromosome_name"),
                filters = "affy_hg_u95av2", values = ids, mart = ensembl_data)
install.packages("dplyr")
library(dplyr)
set = result %>%
  group_by (affy_hg_u95av2) %>% summarise(count = n())
sum(set$count > 1)

#Question 3: How many probesets (Affymetrix IDs) are annotated with one or more genes on the autosomes (chromosomes 1 to 22).
result_auto = subset(result, chromosome_name < 23)
set_auto = result_auto %>%
  group_by (affy_hg_u95av2) %>% summarise(count_auto = n())
sum(set_auto$count_auto > 0)

result_autosome <- subset(result, chromosome_name < 23)
prob_set_autosome <- result_autosome %>%
  group_by (affy_hg_u95av2) %>%
  summarise(
    prob_count = n()
  )
tail(result_autosome)
sum(prob_set_autosome$prob_count > 0)

#Use the MsetEx dataset from the minfiData package. Part of this question is to use the help system to figure out how to address the question.
#Question 4: What is the mean value of the Methylation channel across the features for sample "5723646052_R04C01"?
biocLite("minfiData")
library(minfiData)
data(MsetEx)
head(MsetEx)
pData(MsetEx)
sample = MsetEx[,2]
?`MethylSet-class`
mean(getMeth(sample))

#Question 5: Access the processed data from NCBI GEO Accession number GSE788. What is the mean expression level of sample GSM9024?
biocLite("GEOquery")
library(GEOquery)
elist = getGEO("GSE788")
length(elist)
names(elist)
edata = elist[[1]]
sampleNames(phenoData(edata))
selected = exprs(edata)[,2]
mean(selected)

#We are using the airway dataset from the airway package.
#Question 6: What is the average of the average length across the samples in the expriment?
biocLite("airway")
library(airway)
data(airway)
airway
colData(airway)
colnames(airway)
mean(airway$avgLength)

#We are using the airway dataset from the airway package. The features in this dataset are Ensembl genes.
#Question 7: What is the number of Ensembl genes which have a count of 1 read or more in sample SRR1039512?
SRR1039512 = airway[,3]
Counts = assay(SRR1039512, "counts")
sum(Counts>=1)

#Question 8: The airway dataset contains more than 64k features. How many of these features overlaps with transcripts on the autosomes (chromosomes 1-22) as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?
#Clarification: A feature has to overlap the actual transcript, not the intron of a transcript. So you will need to make sure that the transcript representation does not contain introns.
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
Exons = exons(txdb)
auto = paste0("chr", c(1:22))
Exons = keepSeqlevels(Exons, auto, pruning.mode = "coarse")
ncbiStyleLevels = mapSeqlevels(seqlevels(Exons),"NCBI")
Exons = renameSeqlevels(Exons, ncbiStyleLevels)
ov = subsetByOverlaps(airway,Exons)
ov

#The expression measures of the airway dataset are the number of reads mapping to each feature. In the previous question we have established that many of these features do not overlap autosomal transcripts from the TxDb.Hsapiens.UCSC.hg19.knownGene. But how many reads map to features which overlaps these transcripts?
#Question 9: For sample SRR1039508, how big a percentage (expressed as a number between 0 and 1) of the total reads in the airway dataset for that sample, are part of a feature which overlaps an autosomal TxDb.Hsapiens.UCSC.hg19.knownGene transcript?
SRR1039508 = airway[,1]
subset_SRR1039508 = subsetByOverlaps(SRR1039508, Exons)
counts = assay(SRR1039508, "counts")
subset_counts = assay(subset_SRR1039508, "counts")
sum(subset_counts)/sum(counts)

#Consider sample SRR1039508 and only consider features which overlap autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene. We should be able to very roughly divide these transcripts into expressed and non expressed transcript. Expressed transcripts should be marked by H3K4me3 at their promoter. The airway dataset have assayed "airway smooth muscle cells". In the Roadmap Epigenomics data set, the E096 is supposed to be "lung". Obtain the H3K4me3 narrowPeaks from the E096 sample using the AnnotationHub package.
#Question 10: What is the median number of counts per feature (for sample SRR1039508) containing a H3K4me narrowPeak in their promoter (only features which overlap autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene are considered)?
#Clarification: We are using the standard 2.2kb default Bioconductor promotor setting.
#Conclusion Compare this to the median number of counts for features without a H3K4me3 peak. Note that this short analysis has not taken transcript lengths into account and it compares different genomic regions to each other; this is highly suscepticle to bias such as sequence bias.
library(AnnotationHub)
ah = AnnotationHub()
qah_h1 = query(ah, c("E096", "H3K4me3"))
qah_h1
h1 = qah_h1[[2]]
h1 = keepSeqlevels(h1, auto, pruning.mode = "coarse")
h1 = renameSeqlevels(h1, ncbiStyleLevels)

t = range(rowRanges(subset_SRR1039508))
auto_ncbi = extractSeqlevelsByGroup(species="Homo sapiens", style="NCBI", group="auto")
t = keepSeqlevels(t, auto_ncbi)
p = promoters(t)

ov = subsetByOverlaps(p, h1)
t2 = subsetByOverlaps(subset_SRR1039508, ov)
counts = assay(t2, "counts")
median(counts)
