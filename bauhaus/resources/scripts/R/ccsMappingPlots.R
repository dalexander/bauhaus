library(pbbamr)
library(dplyr)
library(ggplot2)
library(xml2)
library(stringr)
library(feather)

toPhred <- function(acc, maximum=60) {
    err = pmax(1-acc, 10^(-maximum/10))
    -10*log10(err)
}


getConditionTable <- function(wfOutputRoot)
{
    read.csv(file.path(wfOutputRoot, "condition-table.csv"))
}


## This is not a good idea at all really---the dataset could contain
## "filter" operations that we are ignoring via this mechanism.  We
## need a better solution in pbbamr, to somehow provide access to a
## virtual pbi.
listDatasetContents <- function(datasetXmlFile)
{
    x <- read_xml(datasetXmlFile)
    ns <- xml_ns(x)
    allResourceFiles <- sapply(xml_find_all(x, ".//pbbase:ExternalResource/@ResourceId", ns), xml_text)
    isBam <- str_detect(allResourceFiles, ".*.bam$")
    bams <- unique(allResourceFiles[isBam])
    bams
}

makeCCSDataFrame1 <- function(datasetXmlFile, conditionName, sampleFraction=1.0)
{
    print(datasetXmlFile)
    ## Do subsampling at the BAM level
    allBams <- listDatasetContents(datasetXmlFile)
    if (sampleFraction < 1) {
        set.seed(42) # Do we want to do this globally instead?
        n <- max(1, floor(length(allBams)*sampleFraction))
        sampledBams <- as.character(sample_n(data.frame(fname=allBams), n)$fname)
    } else {
        sampledBams <- allBams
    }
    pbis <- lapply(sampledBams, pbbamr::loadPBI,
                   loadSNR = TRUE, loadNumPasses = TRUE, loadRQ = TRUE)
    ## This would be more efficient, but it crashes!
    ##do.call(bind_rows, sampledBams)
    combinedPbi <- do.call(rbind, pbis)

    ## TODO: moviename??

    ## TODO: readlength not yet available, unfortunately, due to the
    ## qstart/qend convention for CCS reads.
    with(combinedPbi,
         tbl_df(data.frame(
             Condition=conditionName,
             NumPasses = np,
             HoleNumber = hole,
             ReadQuality = qual,
             ReadQualityPhred = toPhred(qual),
             Identity = 1. - (mismatches + inserts + dels)/(aend-astart),
             IdentityPhred = toPhred(1. - (mismatches + inserts + dels)/(aend-astart)),
             NumErrors=(mismatches+inserts+dels),
             TemplateSpan=(tend-tstart),
             ReadLength=(aend-astart),          ## <-- this is a lie, see above!
             SnrA = snrA,
             SnrC = snrC,
             SnrG = snrG,
             SnrT = snrT)))
}

makeCCSDataFrame <- function(wfOutputRoot, sampleFraction=1.0)
{
    ct <- getConditionTable(wfOutputRoot)
    conditions <- unique(ct$Condition)
    dsetXmls <- sapply(conditions, function(condition) file.path(wfOutputRoot, condition, "ccs_mapping/all_movies.consensusalignments.xml"))
    dfs <- mapply(makeCCSDataFrame1, dsetXmls, conditions, sampleFraction=sampleFraction, SIMPLIFY=F)
    tbl_df(do.call(rbind, dfs))
}


doCCSCumulativeYieldPlots <- function(ccsDf)
{
  cumByCut <- function(x) {
    qvOrder <- order(x$IdentityPhred, decreasing=TRUE)
    xo <- x[qvOrder,]
    xo$NumReads <- seq(1, nrow(xo))
    xo$YieldFraction <- cumsum(xo$ReadLength) / sum(xo$ReadLength)
    xo[seq(1,nrow(xo), by=10),]
  }

  ## yield <- ddply(ccsDf, "Condition", cumByCut)
  yield <- ccsDf %>% group_by(Condition) %>% do(cumByCut(.))

  ## NumReads on y-axis
  p <- qplot(IdentityPhred, NumReads, colour=Condition, data=yield, main="Yield of reads by CCS accuracy")
  print(p)

  ## Fraction of reads on y-axis
  p <- qplot(IdentityPhred, YieldFraction, colour=Condition, data=yield, main="Fractional yield by CCS accuracy")
  print(p)
}

doCCSNumPassesHistogram <- function(ccsDf)
{
    p <- qplot(NumPasses, data=ccsDf, geom="density", color=Condition,
               main="NumPasses distribution (density)")
    print(p)
}

doCCSNumPassesCDF <- function(ccsDf)
{
    p <- (ggplot(aes(x=NumPasses, color=Condition), data=ccsDf) +
          stat_ecdf(geom="step") +
          ggtitle("NumPasses distribution (ECDF)"))
    print(p)
}


## calibration plot...

doCCSReadQualityCalibrationPlots <- function(ccsDf)
{
    ccsDf <- sample_n(ccsDf, 5000)

    p <- qplot(ReadQuality, Identity, alpha=I(0.1), data=ccsDf) + facet_grid(.~Condition) +
        geom_abline(slope=1, color="red") +
        ggtitle("Read quality versus empirical accuracy")
    print(p)

    p <- qplot(ReadQualityPhred, IdentityPhred, alpha=I(0.1), data=ccsDf) + facet_grid(.~Condition) +
        geom_abline(slope=1, color="red") +
        ggtitle("Read quality versus empirical accuracy (Phred scale)")
    print(p)
}


doCCSTitrationPlots <- function(ccsDf)
{
     accVsNp <- ccsDf %>% group_by(Condition, NumPasses) %>% summarize(
       MeanIdentity=1-(max(1, sum(NumErrors))/sum(ReadLength)),
       TotalBases=sum(ReadLength)) %>% mutate(
       MeanIdentityPhred=toPhred(MeanIdentity))

     p <- qplot(NumPasses, MeanIdentityPhred, size=TotalBases, weight=TotalBases, data=filter(accVsNp, NumPasses<20)) +
         facet_grid(.~Condition) + geom_smooth()
     print(p)
}


doAllCCSPlots <- function(ccsDf)
{
    doCCSTitrationPlots(ccsDf)
    doCCSNumPassesHistogram(ccsDf)
    doCCSNumPassesCDF(ccsDf)
    doCCSReadQualityCalibrationPlots(ccsDf)
    doCCSCumulativeYieldPlots(ccsDf)
}


## Main, when run as a script.
if (!interactive())
{
    args <- commandArgs(TRUE)
    wfRootDir <- args[1]
    ccsDf <- makeCCSDataFrame(wfRootDir)
    ##write.csv(ccsDf, "ccs-mapping.csv")  ## TOO BIG, TOO SLOW
    write_feather(ccsDf, "ccs-mapping.feather")
    pdf("ccs-mapping.pdf", 11, 8.5)
    doAllCCSPlots(ccsDf)
    dev.off()
}


if (0) {
    ##wfRoot = "/home/UNIXHOME/dalexander/Projects/Analysis/EchidnaConsensus/2kLambda_4hr_postTrain_CCS/"
    wfRoot <- "/home/UNIXHOME/ayang/projects/bauhaus/Echidna_PerfVer/EchidnaVer_CCS_postTrain"
    df <- makeCCSDataFrame(wfRoot, 1.0)
}
