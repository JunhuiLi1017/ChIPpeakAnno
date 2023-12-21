#' Oligonucleotide enrichment analysis
#'
#' Test if given oligonucleotide patterns are enriched in peak regions or not.
#'
#' @param filepath The path to the file that contains the oligonucleotide pattern. 
#'
#' @param format The format of file containing the oligonucleotide pattern, either "fasta" (default) or "fastq".
#'
#' @param peaks A GRanges object containing the peaks.
#' 
#' @param upstream The number of base pairs to expand the peak regions (upstream) when searching for oligonucleotide patterns. Default to 0.
#' 
#' @param downstream The number of base pairs to expand the peak regions (downstream) when searching for oligonucleotide patterns. Default to 0.
#'
#' @param genome BSgenome object or mart object. Please refer to available.genomes in BSgenome package and useMart in bioMaRt package for details
#'
#' @param methodBackground The method to get the background of compared oligonucleotide. "select.chr.randomly" (default) is used to select background chromosomes from all chromosomesor, and "shuffle" will  shuffle the letters within input sequences with any k-let size.
#'
#' @param chromosome Specify which chromosome will be selected to randomly pick back ground sequences. Default is the chromosome in peaks. Note that this parameter is valid for 'select.chr.randomly' method.
#'
#' @param ... could be parameters of function \code{\link[universalmotif]{shuffle_sequences}}
#'
#' @param times The times of permutation test, default is 1000
#'
#' @param alpha The significant level for permutation test, default is 0.05
#'
#' @return
#'
#' A data frame with 5 columns as x (number of match of the pattern), n (total number of oligonucleotide with the same length of pattern in the input) and prop.background (the proportions of pattern in background sequence), binom.pvalue (p value for the null that probabilities of success equal certain given values ) and threshold(p value threshold with given times of sample or shuffle).
#'
#' @details
#'
#' Please see \link[universalmotif]{shuffle_sequences} for the more information about 'shuffle' method.
#'
#' @author Junhui Li
#'
#' @importFrom  Biostrings oligonucleotideFrequency
#' 
#' @importFrom  Biostrings readDNAStringSet
#' 
#' @importFrom  Biostrings reverseComplement
#' 
#' @importFrom stats binom.test
#' 
#' @importFrom universalmotif shuffle_sequences
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' filepath =system.file("extdata", "examplePattern.fa", package="ChIPpeakAnno")
#' peaks = GRanges(seqnames=c("chr17","chr3","chr12","chr8"),
#'                  IRanges(start=c(41275784,10076141,4654135,31024288),
#'                          end=c(41276382,10076732,4654728,31024996),
#'                          names=paste0("peak",1:4)))
#' result <- oligoNucleotideEnrichment(filepath=filepath,
#' peaks=peaks,
#' genome=Hsapiens,
#' method="binom.test")
#' @export

oligoNucleotideEnrichment <- function(filepath,
                                      format = "fasta",
                                      revcomp = TRUE,
                                      peaks,
                                      upstream = 0,
                                      downstream = 0,
                                      genome,
                                      method = c("binom.test","permutation.test"),
                                      backgroud = c("select.chr.randomly","shuffle"),
                                      chromosome = NULL,
                                      ...,
                                      times = 1000,
                                      alpha = 0.05){
  stopifnot("File doesn't exist!"=file.exists(filepath))
  oligoNset <- readDNAStringSet(filepath=filepath, format=format, use.names = TRUE)
  patternName <- names(oligoNset)
  forwardpattern <- as.character(oligoNset)
  
  method <- match.arg(method)
  backgroud <- match.arg(backgroud)
  stopifnot("genome is required parameter, 
           please pass in either a BSgenome object or a Mart object."=
              is(genome, "BSgenome") | is(genome, "Mart"))
  stopifnot("The 'times' or 'alpha' parameter should be increased to a sufficient extent."=round(times*alpha,0) != 0)
  nPeak <- length(peaks)
  stopifnot("There is no peak, please check your data."=nPeak != 0)
  peakSeq <- getAllPeakSequence(peaks, upstream = upstream, downstream = downstream, genome = genome)
  if(method == "binom.test"){
    peakPvalue <- binomEnrichment(forwardpattern, revcomp=revcomp, seqs = peakSeq)
  }else if(method == "permutation.test"){
    peakPvalue <- permutationEnrichment(forwardpattern, revcomp=revcomp, seqs = peakSeq, genome=genome, backgroud = backgroud, times=times)
  }
  return(peakPvalue)
}

binomEnrichment <- function(forwardpattern, revcomp=TRUE, seqs){
  ## get frequency of single nuleotide in seqs
  ACGTcount <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),width=1))
  ACGTfreq <- ACGTcount/sum(ACGTcount)
  allBaseFreq <- ACGTfreq
  allBaseFreq['D'] <- sum(ACGTfreq[c("G","A","T")])
  allBaseFreq['H'] <- sum(ACGTfreq[c("C","A","T")])
  allBaseFreq['B'] <- sum(ACGTfreq[c("G","C","T")])
  allBaseFreq['V'] <- sum(ACGTfreq[c("G","A","C")])
  allBaseFreq['W'] <- sum(ACGTfreq[c("A","T")])
  allBaseFreq['S'] <- sum(ACGTfreq[c("G","C")])
  allBaseFreq['K'] <- sum(ACGTfreq[c("G","T")])
  allBaseFreq['M'] <- sum(ACGTfreq[c("C","A")])
  allBaseFreq['Y'] <- sum(ACGTfreq[c("C","T")])
  allBaseFreq['R'] <- sum(ACGTfreq[c("G","A")])
  allBaseFreq['N'] <- 1
  
  ## get expected frequency of fasta
  forwardpattern.t <- translatePattern(forwardpattern)
  backwardpattern.t <- NULL
  backwardpattern <- NULL
  if(revcomp == TRUE){
    backwardpattern <- reverseComplement(DNAStringSet(forwardpattern))
    backwardpattern <- as.character(backwardpattern)
    backwardpattern.t <- translatePattern(backwardpattern)
  }
  
  ## get expected frequency of forward and backward fasta 
  expFreqSet <- sapply(c(forwardpattern,backwardpattern),function(i){
    patternTab <- table(strsplit(i,"")[[1]])
    expFre <- prod(allBaseFreq[names(patternTab)]^patternTab)
    expFre
  })
  expFreqSet <- colSums(matrix(expFreqSet,ncol=length(forwardpattern),byrow = TRUE))
  names(expFreqSet) <- names(forwardpattern)
  
  ## get real frequency of input seq in peak region
  peakFreq <- matrix(0,length(forwardpattern),3)
  rownames(peakFreq) <- names(forwardpattern)
  colnames(peakFreq) <- c("motif_num","all_motif_num","motif_background_rate")
  peakFreq[,3] <- expFreqSet
  
  peakFreq[,1:2] <- t(sapply(seq.int(forwardpattern),function(i){
    oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$sequence),width=nchar(forwardpattern[i])))
    posPlus = gregexpr(forwardpattern.t[i], names(oligoNVec), perl = TRUE)
    posPlusFreq <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
    posMinusFreq <- 0
    if(revcomp == TRUE){
      posMinus = gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
      posMinusFreq <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
    }
    c(posPlusFreq + posMinusFreq, sum(oligoNVec))
  }))

  cut_off_binom_test <- sapply(seq.int(nrow(peakFreq)),function(i){
    stats::binom.test(peakFreq[i,1],peakFreq[i,2],peakFreq[i,3],alternative="greater")$p.value
  })
  peakFreq <- data.frame(peakFreq,cut_off_binom_test)
  return(peakFreq)
}

permutationEnrichment <- function(forwardpattern, revcomp = TRUE, seqs, upstream = 0, downstream = 0, genome, backgroud = c("select.chr.randomly","shuffle"), chromosome=NULL, times = 1000, alpha = 0.05){
  forwardpattern.t <- translatePattern(forwardpattern)
  backwardpattern.t <- NULL
  backwardpattern <- NULL
  if(revcomp == TRUE){
    backwardpattern <- reverseComplement(DNAStringSet(forwardpattern))
    backwardpattern <- as.character(backwardpattern)
    backwardpattern.t <- translatePattern(backwardpattern)
  }
  allseqLen <- seqlengths(genome)
  seqWidth <- width(seqs)
  nPeak <- length(seqs)
  colnames(mcols(seqs))[3] <- "orig.sequence"
  names(seqs$orig.sequence) <- do.call(paste, c(as.data.frame(seqs)[,1:3], sep="_"))
  inputSeq <- DNAStringSet(seqs$orig.sequence)
  candiChrom <- names(allseqLen >= max(seqWidth))
  stopifnot("The length of peak sequence should be less than the length of chromosomes"=
              all(seqWidth < allseqLen[as.character(seqnames(seqs)@values)]))

  peakFreq <- matrix(0,length(forwardpattern),4)
  rownames(peakFreq) <- names(forwardpattern)
  colnames(peakFreq) <- c("motif_num","all_motif_num","motif_rate","cut_off_permutation_test")
  peakFreq[,1:2] <- t(sapply(seq.int(forwardpattern),function(i){
    oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(seqs$orig.sequence),width=nchar(forwardpattern[i])))
    posPlus = gregexpr(forwardpattern.t[i], names(oligoNVec), perl = TRUE)
    posPlusFreq <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
    posMinusFreq <- 0
    if(revcomp == TRUE){
      posMinus = gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
      posMinusFreq <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
    }
    c(posPlusFreq + posMinusFreq, sum(oligoNVec))
  }))
  peakFreq[,3] <- peakFreq[,1]/peakFreq[,2]
  
  permutationFreq <- do.call(rbind,lapply(seq.int(times),function(n) {
    if(backgroud == "select.chr.randomly"){
      if(is.null(chromosome)){
        chrs <- as.character(seqnames(seqs))
        givenSeqLen <- allseqLen[chrs]
      }else{
        chrs <- sample(candiChrom,nPeak)
        givenSeqLen <- allseqLen[chrs]
      }
      startPos <- vapply(givenSeqLen-seqWidth,function(x){sample(seq.int(x),1)},numeric(1))
      endPos <- startPos + seqWidth - 1
      backgroudPeak <- GRanges(seqnames = chrs,
                                IRanges(start = startPos,
                                        end = endPos,
                                        names = paste0("peak",seq.int(nPeak))))
      backgroudPeakseq <- getAllPeakSequence(backgroudPeak, upstream = upstream, downstream = downstream, genome = genome)
    }else if(backgroud == "shuffle"){
      sequences.shuffled <- universalmotif::shuffle_sequences(inputSeq, ...)
      #sequences.shuffled <- universalmotif::shuffle_sequences(inputSeq)
      peakSeq$sequence <- as.vector(sequences.shuffled)
      backgroudPeakseq <- peakSeq
    }
    ## get frequency from permutation
    sapply(seq.int(forwardpattern),function(i){
      oligoNVec <- colSums(oligonucleotideFrequency(DNAStringSet(backgroudPeakseq$sequence),width=nchar(forwardpattern[i])))
      posPlus <- gregexpr(forwardpattern.t[i], names(oligoNVec), perl = TRUE)
      expFre <- sum(oligoNVec[which(sapply(posPlus, "[[", 1) > 0)])
      if(revcomp == TRUE){
        posMinus <- gregexpr(backwardpattern.t[i], names(oligoNVec), perl = TRUE)
        expFreMinus <- sum(oligoNVec[which(sapply(posMinus, "[[", 1) > 0)])
        expFre <- expFre + expFreMinus
      }
      c(expFre/sum(oligoNVec))
    })
  }))
  colnames(permutationFreq) <- names(forwardpattern)
  permutationFreqSorted <- apply(permutationFreq, 2, sort, decreasing=TRUE)
  peakFreq[,4] <- as.vector(permutationFreqSorted[round(times*alpha,0),])
  return(peakFreq)
}


