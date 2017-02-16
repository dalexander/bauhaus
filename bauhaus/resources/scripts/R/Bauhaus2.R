library(data.table)
library(dplyr)
library(ggplot2)
library(jsonlite)

##
## Core functions for operation in a bauhaus2/zia environment.
##

# Random utils

chkClass <- function(var, classname, msg) {
  if (!(classname %in% class(var))) {
    stop(msg)
  }
}

#' Check the file ends as a png
chkPng <- function(fname) {
  substr(fname, nchar(fname) - 3, nchar(fname)) == '.png'
}

#' Check the file ends as a csv
chkCsv <- function(fname) {
  substr(fname, nchar(fname) - 3, nchar(fname)) == '.csv'
}



#' Returns a default there for making plots.
#' @export
getPBTheme <- function() {
  return(ggplot2::theme_bw(base_size = 14))
}


defaultPalette = "Set1"
defaultLargePaletteGetter = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

#' Get the default color scheme, adding more colors if they are available.
#' @export
getPBColorScale <- function(numLevels = 9) {
  if (numLevels <= 9) {
    return(ggplot2::scale_colour_brewer(palette = defaultPalette))
  } else {
    return(ggplot2::scale_color_manual(values = defaultLargePaletteGetter(numLevels)))
  }
}


#' Get the default color scheme for fills, adding more colors if the current
#' palette is maxed out
#' @export
getPBFillScale <- function(numLevels = 9) {
  if (numLevels <= 9) {
    return(ggplot2::scale_fill_brewer(palette = "Set1"))
  } else {
    return(ggplot2::scale_fill_manual(values = defaultLargePaletteGetter(numLevels)))
  }
}



##
## Condition table grokking.
##
loadConditionTable <- function(conditionTableCSV,
                               wfOutRoot=dirname(conditionTableCSV))
{
    ## We load the Bauhaus2 CSV and augment it with extra columns R analysis is
    ## going to need:
    ##  - Reference      (path to reference FASTA)
    ##  - MappedSubreads (path to the mapped subreads)
    ##  - more??
    ## (These files should all found in the wfOutRoot)
    ct <- read.csv(conditionTableCSV)
    data.table(ct)
    mappedSubreads <- file.path(wfOutRoot, "conditions", ct$Condition, "mapped/mapped.alignmentset.xml")
    ct$MappedSubreads <- mappedSubreads
    referenceFasta <- file.path(wfOutRoot, "conditions", ct$Condition, "reference.fasta")
    ct$Reference <- referenceFasta
    ct
}


##
## Generating reports about plots, tables
##
bh2Reporter <- function(conditionTableCSV, outputFile, reportTitle) {

    reportOutputDir <- dirname(outputFile)
    reportOutputFile = outputFile
    version <- version
    plotsToOutput <- NULL
    tablesToOutput <- NULL
    cond_table <- loadConditionTable(conditionTableCSV)

    ## Save a ggplot in the report.
    .ggsave <- function(img_file_name, plot, id="Default ID", title="Default Title",
                       caption="No caption specified", tags=c(), ...)
    {
        if (!chkPng(img_file_name)) {
            img_file_name = paste0(img_file_name, ".png")
        }
        img_path = file.path(reportOutputDir, img_file_name)
        ggplot2::ggsave(img_path, plot = plot, ...)
        logging::loginfo(paste("Wrote img to: ", img_path))

        thisPlot <- list(id=unbox(id),
                         image=unbox(img_file_name),
                         title=unbox(title),
                         caption=unbox(caption),
                         tags=as.vector(tags, mode="character"))
        plotsToOutput <<- rbind(plotsToOutput, thisPlot)
    }

    ## Add a table to the report.
    .write.table<- function(tbl_file_name, tbl, id, title = "Default Title", tags=c())
    {
        if (!chkCsv(tbl_file_name)) {
            tbl_file_name = paste(tbl_file_name, ".csv")
        }
        tbl_path = file.path(reportOutputDir, tbl_file_name)
        write.csv(tbl, file=tbl_path)
        logging::loginfo(paste("Wrote table to: ", tbl_path))

        thisTbl <- list(id=unbox(id),
                        csv=unbox(tbl_file_name),
                        title=unbox(title),
                        tags=as.vector(tags, mode="character"))

        tablesToOutput <<- rbind(tablesToOutput, thisTbl)
    }

    ## Output the report file as json.
    .write.report <- function()
    {
        pp <- as.data.frame(plotsToOutput)
        row.names(pp) <- NULL
        tt <- as.data.frame(tablesToOutput)
        row.names(tt) <- NULL
        write_json(list(plots=pp, tables=tt), reportOutputFile, pretty=T)
        logging::loginfo(paste("Wrote report to: ", reportOutputFile))
    }


    list(condition.table = cond_table,
         ggsave = .ggsave,
         write.table = .write.table,
         write.report = .write.report,
         outputDir = reportOutputDir,
         outputJSON = reportOutputFile)
}




if (0) {

    r <- bh2Reporter("~dalexander/Projects/rsync/bauhaus/test/data/two-tiny-movies.csv", "/tmp/report.json")

    p <- qplot()
    r$ggsave("1.png", p, "Id1", "Title1", "Caption1", tags=c())
    r$ggsave("2.png", p, "Id2", "Title2", "Caption2", tags=c("A"))

    t <- data.frame(a=1:3, b=4:6)
    r$write.table("foo.csv", t, "T1", "TTitle1", tags=c("A"))
    r$write.table("foo2.csv", t, "T2", "TTitle2", tags=c())

    r$write.report()
}
