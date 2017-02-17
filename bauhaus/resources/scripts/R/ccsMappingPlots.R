library(pbbamr)
library(dplyr)
library(ggplot2)
library(xml2)
library(stringr)
# library(feather)

## FIXME: make a real package
myDir = "./scripts/R"
source(file.path(myDir, "Bauhaus2.R"))

toPhred <- function(acc, maximum=60) {
    err = pmax(1-acc, 10^(-maximum/10))
    -10*log10(err)
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

makeCCSDataFrame <- function(report, wfOutputRoot, sampleFraction=1.0)
{
    ct <- report$condition.table
    conditions <- unique(ct$Condition)
    dsetXmls <- sapply(conditions, function(condition) file.path(wfOutputRoot, condition, "ccs_mapping/all_movies.consensusalignments.xml"))
    dfs <- mapply(makeCCSDataFrame1, dsetXmls, conditions, sampleFraction=sampleFraction, SIMPLIFY=F)
    tbl_df(do.call(rbind, dfs))
}


doCCSCumulativeYieldPlots <- function(report, ccsDf)
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
  report$ggsave(
    "yield_reads_ccs_accuracy.png",
    p,
    id = "yield_reads_ccs_accuracy",
    title = "Yield of reads by CCS accuracy",
    caption = "Yield of reads by CCS accuracy"
  )

  ## Fraction of reads on y-axis
  p <- qplot(IdentityPhred, YieldFraction, colour=Condition, data=yield, main="Fractional yield by CCS accuracy")
  report$ggsave(
    "fractional_yield_ccs_accuracy.png",
    p,
    id = "fractional_yield_ccs_accuracy",
    title = "Fractional yield by CCS accuracy",
    caption = "Fractional yield by CCS accuracy"
  )
}

doCCSNumPassesHistogram <- function(report, ccsDf)
{
    p <- qplot(NumPasses, data=ccsDf, geom="density", color=Condition,
               main="NumPasses distribution (density)")
    report$ggsave(
      "numpasses_dist_density.png",
      p,
      id = "numpasses_dist_density",
      title = "NumPasses distribution (density)",
      caption = "NumPasses distribution (density)"
    )
}

doCCSNumPassesCDF <- function(report, ccsDf)
{
    p <- (ggplot(aes(x=NumPasses, color=Condition), data=ccsDf) +
          stat_ecdf(geom="step") +
          ggtitle("NumPasses distribution (ECDF)"))
    report$ggsave(
      "numpasses_dist_ecdf.png",
      p,
      id = "numpasses_dist_ecdf",
      title = "NumPasses distribution (ECDF)",
      caption = "NumPasses distribution (ECDF)"
    )
}


## calibration plot...

doCCSReadQualityCalibrationPlots <- function(report, ccsDf)
{
    ccsDf <- sample_n(ccsDf, min(5000, nrow(ccsDf)))

    p <- qplot(ReadQuality, Identity, alpha=I(0.1), data=ccsDf) + facet_grid(.~Condition) +
        geom_abline(slope=1, color="red") +
        ggtitle("Read quality versus empirical accuracy")
    report$ggsave(
      "read_quality_vs_empirical_accuracy.png",
      p,
      id = "read_quality_vs_empirical_accuracy",
      title = "Read quality versus empirical accuracy",
      caption = "Read quality versus empirical accuracy"
    )

    p <- qplot(ReadQualityPhred, IdentityPhred, alpha=I(0.1), data=ccsDf) + facet_grid(.~Condition) +
        geom_abline(slope=1, color="red") +
        ggtitle("Read quality versus empirical accuracy (Phred scale)")
    report$ggsave(
      "read_quality_vs_empirical_accuracy_phred.png",
      p,
      id = "read_quality_vs_empirical_accuracy_phred",
      title = "Read quality versus empirical accuracy (Phred scale)",
      caption = "Read quality versus empirical accuracy (Phred scale)"
    )
}


doCCSTitrationPlots <- function(report, ccsDf)
{
     accVsNp <- ccsDf %>% group_by(Condition, NumPasses) %>% summarize(
       MeanIdentity=1-(max(1, sum(NumErrors))/sum(ReadLength)),
       TotalBases=sum(ReadLength)) %>% mutate(
       MeanIdentityPhred=toPhred(MeanIdentity))

     p <- qplot(NumPasses, MeanIdentityPhred, size=TotalBases, weight=TotalBases, data=filter(accVsNp, NumPasses<20)) +
         facet_grid(.~Condition) + geom_smooth()
     report$ggsave(
       "ccs_titration.png",
       p,
       id = "ccs_titration",
       title = "CCS Titration Plots",
       caption = "CCS Titration Plots"
     )
}


doAllCCSPlots <- function(report, ccsDf)
{
    doCCSTitrationPlots(report, ccsDf)
    doCCSNumPassesHistogram(report, ccsDf)
    doCCSNumPassesCDF(report, ccsDf)
    doCCSReadQualityCalibrationPlots(report, ccsDf)
    doCCSCumulativeYieldPlots(report, ccsDf)
}

makeReport <- function(report) {
  if (!interactive()){
    args <- commandArgs(TRUE)
    wfRootDir <- args[1]
    ccsDf <- makeCCSDataFrame(report, wfRootDir)  
    # write_feather(ccsDf, "ccs-mapping.feather")
    doAllCCSPlots(report, ccsDf)
    
    # Save the report object for later debugging
    save(report, file = file.path(report$outputDir, "report.Rd"))
    
    # At the end of this function we need to call this last, it outputs the report
    report$write.report()
  }
  if (0) {
    ##wfRoot = "/home/UNIXHOME/dalexander/Projects/Analysis/EchidnaConsensus/2kLambda_4hr_postTrain_CCS/"
    wfRoot <- "/home/UNIXHOME/ayang/projects/bauhaus/Echidna_PerfVer/EchidnaVer_CCS_postTrain"
    df <- makeCCSDataFrame(report, wfRoot, 1.0)
  }
}

main <- function()
{
  report <- bh2Reporter(
    "condition-table.csv",
    "reports/CCSMappingReports/report.json",
    "CCS Mapping Reports")
  makeReport(report)
}

## Leave this as the last line in the file.
logging::basicConfig()
main()