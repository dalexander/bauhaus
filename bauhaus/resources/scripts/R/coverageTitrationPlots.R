library(ggplot2)
library(plyr)
library(stringr)
library(pbutils)
library(gtools) # for "permutations"

ignore <- function(x) {}

.squishWhitespace <- function(s) {
  str_replace_all(s, "[[:space:]]+", " ")
}

## Abbrev: wf == "workflow"

## TODO: functions like these should move to pbexperiment or whatnot, or at least a
##       common.R or somesuch

getConditionTable <- function(wfOutputRoot)
{
    read.csv(file.path(wfOutputRoot, "condition-table.csv"))
}

variableNames <- function(ct)
{
    nms <- names(ct)
    matches <- str_detect(nms, "(Algorithm)|(Genome)|(p_.*)")
    nms[matches]
}

makeCoverageTitrationTable1 <- function(maskedVariantsDir)
{
    ##
    ## Get the precomputed masked files, which we expect to find as
    ##   {worfklowOutputRoot}/{condition}/variant_calling/{algo}/masked-variants-{coverage}.gff
    ##
    alnSummary <- file.path(maskedVariantsDir, "../alignments_summary.gff")
    mvs <- Sys.glob(file.path(maskedVariantsDir, "masked-variants-*.gff"))

    algo <- basename(maskedVariantsDir)
    condition <- basename(dirname(dirname(maskedVariantsDir)))
    coverageFromPath   <- function(path) as.integer(str_match(path, ".*/masked-variants-(.*)\\.gff")[,2])

    variantsCount.fast <- function(fname) {
        CMD <- sprintf("egrep -v '^#' %s | wc -l", fname)
        result <- as.integer(.squishWhitespace(system(CMD, intern=T)))
        result
    }

    empiricalQV <- function(numErrors, genomeSize) {
        err <- (numErrors + 1)/(genomeSize + 1)
        -10*log10(err)
    }

    computeGenomeSize <- function(alignmentSummaryFile)
    {
        lines <- readLines(alignmentSummaryFile)
        sequenceRegionLines <- lines[grep("##sequence-region", lines)]
        cols <- str_split(sequenceRegionLines, " ")
        chromosomeSizes <- as.integer(lapply(cols, tail, 1))
        sum(chromosomeSizes)
    }

    computeQ01Coverage <- function(alignmentSummaryFile) {
        cov2 <- readGFF(alignmentSummaryFile)$cov2
        cov <- as.numeric(str_extract(cov2, "[^,]*"))
        as.numeric(quantile(cov, 0.01))
    }

    q01Coverage <- computeQ01Coverage(alnSummary)
    numVariants <- sapply(mvs, variantsCount.fast)
    genomeSize <- computeGenomeSize(alnSummary)
    concordanceQV <- empiricalQV(numVariants, genomeSize)

    ctt <- data.frame(
        MaskedVariantsFile  = mvs,
        Coverage            = coverageFromPath(mvs),
        AvailableCoverage   = q01Coverage,
        Condition           = condition,
        NumVariants         = sapply(mvs, variantsCount.fast),
        ConcordanceQV       = concordanceQV,
        Algorithm           = algo)

    ctt$ShouldCensor = (ctt$Coverage > ctt$AvailableCoverage)
    ctt
}

makeCoverageTitrationTable <- function(wfOutputRoot)
{
    ct <- getConditionTable(wfOutputRoot)

    ## dig around the wf output to find masked-variants files.
    mvDirs <- Sys.glob(file.path(wfOutputRoot, "*/variant_calling/*/"))

    rawCtt <- ldply(mvDirs,  function(d) { makeCoverageTitrationTable1(d) })

    ## Get the table of just the conditions and associated variables---no runcodes/other inputs
    ## Is this a concept that is more broadly useful?
    keepColumns = append("Condition", variableNames(ct))
    condensedCt <- unique(ct[,keepColumns])

    merge(rawCtt, condensedCt, by="Condition")
}


doTitrationPlots <- function(tbl)
{
    tbl <- tbl[!tbl$ShouldCensor,]

    ## Implicit variables: Genome, Algorithm
    ## Explicit variables: have the "p_" prefix
    variables <- names(tbl)[grep("^p_|^Genome$|^Algorithm$", names(tbl))]
    nvals <- lapply(variables, function(x) length(unique(tbl[,x])))
    variables <- variables[nvals > 1]
    stopifnot(length(variables) %in% c(1,2,3,4))

    print(ggplot(tbl, aes(x=Coverage, y=ConcordanceQV, color=Condition)) +
          geom_line() +
          scale_y_continuous("Concordance (QV)") +
          #scale_color_manual(values = extras$colors) +
          ggtitle("Consensus Performance, all conditions"))

    # Facet on individual variables
    for (v in variables) {
        print(
            ggplot(tbl, aes(x=Coverage, y=ConcordanceQV, color=Condition)) +
            geom_line() +
            facet_grid(paste(".~", v)) +
            scale_y_continuous("Concordance (QV)") +
            #scale_color_manual(values = extras$colors) +
            ggtitle(paste("Consensus Performance By", v)))
    }

    if (length(variables) >= 2)
    {
      # Take pairs of variables, facet on first and color by second
      apply(permutations(n=length(variables), r=2, v=variables), 1,
            function(twoVars)
            {
              ignore <- print(
                ggplot(tbl, aes_string(x="Coverage", y="ConcordanceQV", color=twoVars[2], group="Condition")) +
                  geom_line() +
                  facet_grid(paste(".~", twoVars[1])) +
                  scale_y_continuous("Concordance (QV)") +
                  #scale_color_manual(values = extras$colors) +
                  ggtitle(paste("Consensus Performance By", twoVars[1])))
            })
    }

    # Take triples of variables, facet on first two and color by third
    if (length(variables) >= 3)
    {
      # Take pairs of variables, facet on first and color by second
      apply(permutations(n=length(variables), r=3, v=variables), 1,
            function(twoVars)
            {
              print(
                ggplot(tbl, aes_string(x="Coverage", y="ConcordanceQV", color=twoVars[3], group="Condition")) +
                  geom_line() +
                  facet_grid(paste(twoVars[1], "~", twoVars[2])) +
                  scale_y_continuous("Concordance (QV)") +
                  #scale_color_manual(values = extras$colors) +
                  ggtitle(paste0("Consensus Performance By ", twoVars[1], ", ", twoVars[2])))
            })
    }
}


doCoverageDiagnosticsPlot <- function(tbl)
{
    ## Diagnostic plot of coverage
    print(qplot(Condition, AvailableCoverage, data=tbl, color=Condition) +
          theme(axis.text.x = element_text(angle=90, hjust=1)) +
          geom_hline(yintercept=80, col="red") +
          ggtitle("1-percentile coverage level by Condition"))
}



makeResidualsTable <- function(ccsDf, variables)
{
    ## Plots of residual error modes
    tbl <- tbl[!tbl$ShouldCensor,]

    summarizedResidualErrors <- function(variantsGffFile) {
        print(variantsGffFile)
        gff <- readGFF(variantsGffFile)
        if (nrow(gff) > 0) {
            base <- ifelse(gff$reference == ".", gff$variantSeq, gff$reference)
            base <- str_sub(base, 1, 1)
            ##count(interaction(gff$type, base))
            df <- as.data.frame(table(interaction(base, gff$type)))
            df$Base      <- str_sub(df$Var1, 1, 1)
            df$ErrorMode <- str_sub(df$Var1, 3)
            df
        }
    }

    #MIN.COVERAGE <- 40
    ROUND.COVERAGE <- c(40, 60, 100)
    varTypeCounts <- ddply(
        subset(tbl, Coverage %in% ROUND.COVERAGE),
        append(c("Condition", "Algorithm", "Coverage", "MaskedVariantsFile", "Genome"), variables),
        function(df) {
            res <- summarizedResidualErrors(as.character(df$MaskedVariantsFile))
        })
}

doResidualErrorsPlot <- function(varTypeCounts, variables)
{
    for (coverage in unique(varTypeCounts$Coverage)) {
        plt <- (ggplot(data=varTypeCounts[varTypeCounts$Coverage==coverage,],
                      mapping=aes(x=ErrorMode, y=Freq, fill=Base)) + geom_bar(stat="identity") +
                theme(axis.text.x = element_text(angle = 75, hjust = 1)))

        facet.formula <- as.formula(".~Condition")
        plt.fix_y  <- (plt + facet_grid(facet.formula)
                       +  ggtitle(sprintf("Residual errors in %dx consensus sequence (fixed y-axis)", coverage)))
        plt.free_y <- (plt + facet_grid(facet.formula, scale="free_y")
                       + ggtitle(sprintf("Residual errors in %dx consensus sequence (free y-axis)", coverage)))

        print(plt.fix_y)
        print(plt.free_y)
    }
}


if (!interactive())
{
    args <- commandArgs(TRUE)
    wfRootDir <- args[1]
    tbl <- makeCoverageTitrationTable(wfRootDir)
    variables <- names(tbl)[grep("^p_", names(tbl))]
    residualsTable <- makeResidualsTable(tbl, variables)

    ## Dump the ctt table
    write.csv(tbl, file="coverage-titration.csv")

    ## Generate plot
    pdf("coverage-titration.pdf", 11, 8.5)
    doTitrationPlots(tbl)
    doCoverageDiagnosticsPlot(tbl)
    doResidualErrorsPlot(residualsTable, variables)
    dev.off()
}




if(0) {

    wfRootDir <- "/home/UNIXHOME/ayang/projects/bauhaus/Echidna_PerfVer/6kecoli_2hrImmob_postTrain_CoverageTitration"
    tbl <- makeCoverageTitrationTable(wfRootDir)

    pdf("/tmp/coverage-titration.pdf", 11, 8.5)
    variables <- names(tbl)[grep("^p_", names(tbl))]
    residualsTable <- makeResidualsTable(tbl)
    doResidualErrorsPlot(residualsTable, variables)
    dev.off()


}
