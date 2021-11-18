setwd("~/chipseq_exercise")


library(ChIPpeakAnno)


case1 <- toGRanges("~/datasets_class/chipseq/peaks/class_chip_peaks.narrowPeak", format="narrowPeak", header = FALSE)
case2 <- toGRanges("~/datasets_class/chipseq/peaks/class_chip2_peaks.narrowPeak", format="narrowPeak", header = FALSE)

overlaps <- findOverlapsOfPeaks(case1, case2)
overlaps <- addMetadata(overlaps, colNames="score", FUN=mean)

makeVennDiagram(overlaps, fill=c("#009E73", "#F0E442"), # circle fill color
                col=c("#D55E00", "#0072B2"), #circle border color
                cat.col=c("#D55E00", "#0072B2")) # label color, keep same as circle border color


library(EnsDb.Hsapiens.v86) ##(hg38)
## create annotation file from EnsDb or TxDb
annoData <- readRDS("/home/rstudio/datasets_class/chipseq/peaks/annoData.RData")

overlappingPeaks <- overlaps$peaklist[["case1///case2"]]

binOverFeature(overlappingPeaks, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               xlab="distance from TSS (bp)", ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")




library(TxDb.Hsapiens.UCSC.hg38.knownGene)
peaks <- GRangesList(rep1=case1,
                     rep2=case2)

genomicElementDistribution(peaks, 
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))



genomicElementDistribution(overlappingPeaks, 
                           TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))



library(UpSetR)

upsetData <- genomicElementUpSetR(overlappingPeaks, 
                                  TxDb.Hsapiens.UCSC.hg38.knownGene)
upset(upsetData$plotData, nsets=13, nintersects=NA)





overlappingPeaks <- annotatePeakInBatch(overlappingPeaks, 
                                        AnnotationData=annoData, 
                                        output="nearestBiDirectionalPromoters",
                                        bindingRegion=c(-2000, 500))
library(org.Hs.eg.db)
overlappingPeaks.anno <- addGeneIDs(overlappingPeaks,
                                    "org.Hs.eg.db",
                                    IDs2Add = c("symbol"))
head(overlappingPeaks.anno)
write.csv(as.data.frame(unname(overlappingPeaks.anno)), "peakes_annotation.csv")

## fix necessaria per un bug nel pacchetto
library(DBI)

over <- getEnrichedGO(overlappingPeaks.anno, 
                      orgAnn="org.Hs.eg.db",
                      maxP = 0.05,
                      minGOterm = 10,
                      condense=TRUE)


## for some reasons go.term which is used to plot is left blank
## so we need to fill it in manually

library(GO.db)
xx <- as.list(GOTERM)
over$bp$go.term <- unlist(lapply(which(over$bp$go.id %in% names(xx)), function(z) Term(xx[[z]])))
over$mf$go.term <- unlist(lapply(which(over$mf$go.id %in% names(xx)), function(z) Term(xx[[z]])))
over$cc$go.term <- unlist(lapply(which(over$cc$go.id %in% names(xx)), function(z) Term(xx[[z]])))

enrichmentPlot(over)
