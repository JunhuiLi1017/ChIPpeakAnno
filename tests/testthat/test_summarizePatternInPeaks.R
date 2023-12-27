test_that("summarizePatternInPeaks works not correct", {
filepath =system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
peaks = GRanges(seqnames=c("chr17","chr3","chr12","chr8"),
                 IRanges(start=c(41275784,10076141,4654135,31024288),
                         end=c(41276382,10076732,4654728,31024996),
                         names=paste0("peak",1:4)))
result <- summarizePatternInPeaks(patternFilePath=filepath,
peaks=peaks,
BSgenomeName=Hsapiens)

expect_equal(result$motif_enrichment[1:2,1], c(2,9))
expect_equal(result$motif_enrichment[1:2,2], c(2474,2474))
expect_equal(result$motif_enrichment[1:2,3], c(0.0003712243,0.0003455440),
             tolerance=1e-6)
expect_equal(result$motif_enrichment[1:2,4], c(0.2342546,3.086031e-07),
             tolerance=1e-5)

expect_equal(as.character(result$motif_occurence[,1]), c("chr17","chr8","chr17",
                                                         "chr17","chr3","chr12",
                                                         "chr12","chr12","chr12"
                                                         ,"chr8","chr8"))
expect_equal(result$motif_occurence[,2], c(41275801,31024627,41276253,41276355,
                                           10076636,4654256,4654355,4654604,
                                           4654653,31024477,31024934))

expect_equal(result$motif_occurence[,3], c(41275806,31024632,41276258,41276360,
                                           10076641,4654261,4654360,4654609,
                                           4654658,31024482,31024939))
expect_equal(result$motif_occurence[,8], c("CGGACC","GGTCCT","AACCAA","AACCCC",
                                           "GTGGTT","AACCCA","GTGGTT","GGGGTT",
                                           "GTGGTT","TAGGTT","GAGGTT"))
})

