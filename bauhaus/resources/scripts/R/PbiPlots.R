#!/usr/bin/env Rscript
# Test Simple Comparison Script
# This script takes a set of conditions and produces a png for display

library(argparser, quietly = TRUE)
library(data.table, quietly = TRUE)
library(jsonlite, quietly = TRUE)
library(logging)
library(ggplot2)
library(pbbamr)
library(uuid, quietly = TRUE)
library(gridExtra)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(survival)
library(ggfortify)

## FIXME: make a real package
myDir = "./scripts/R"
source(file.path(myDir, "Bauhaus2.R"))

#' Define a basic addition to all plots
plTheme <- theme_bw(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))
clScale <- NULL #scale_colour_brewer(palette = "Set1")
clFillScale <- NULL# scale_fill_brewer(palette = "Set1")
themeTilt = theme(axis.text.x = element_text(angle = 45, hjust = 1))
sampleSize = 5000
plotwidth = 7.2
plotheight = 4.2

makeReadLengthSurvivalPlots <- function(report, cd) {
  loginfo("Making Template Span Survival Plots")
  
  cd2 = cd %>% dplyr::group_by(hole, Condition) %>% dplyr::summarise(tlen = sum(tlen))
  cd2 <- as.data.frame(cd2)
  
  cd2$SurvObj <- with(cd2, Surv(tlen))
  cd2.by.con <- survfit(SurvObj ~ Condition, data = cd2)
  p1 <-
    autoplot(cd2.by.con) + labs(x = "Template Span", title = "Template Span Survival") + plTheme
  p2 <-
    autoplot(cd2.by.con) + scale_x_log10() + labs(x = "Template Span", title = "Template Span Survival (Log-scale)") + plTheme
  #  tp <- arrangeGrob(p1, p2, nrow = 2)
  
  report$ggsave(
    "template_span_survival.png",
    p1,
    width = plotwidth,
    height = plotheight,
    id = "template_span_survival",
    title = "Template Span Survival",
    caption = "Template Span Survival"
  )
  report$ggsave(
    "template_span_survival (Log-scale).png",
    p2,
    width = plotwidth,
    height = plotheight,
    id = "template_span_survival(log)",
    title = "Template Span Survival (Log-scale)",
    caption = "Template Span Survival (Log-scale)"
  )
  
  loginfo("Making Aligned Read Length Survival Plots")
  
  cd2 = cd %>% dplyr::group_by(hole, Condition) %>% dplyr::summarise(alen = sum(alen))
  cd2 <- as.data.frame(cd2)
  
  cd2$SurvObj <- with(cd2, Surv(alen))
  cd2.by.con <- survfit(SurvObj ~ Condition, data = cd2)
  p1 <-
    autoplot(cd2.by.con) + labs(x = "Aligned Read Length", title = "Aligned Read Length Survival") + plTheme
  p2 <-
    autoplot(cd2.by.con) + scale_x_log10() + labs(x = "Aligned Read Length", title = "Aligned Read Length Survival (Log-scale)") + plTheme
  #  tp <- arrangeGrob(p1, p2, nrow = 2)
  
  report$ggsave(
    "aligned_read_length_survival.png",
    p1,
    width = plotwidth,
    height = plotheight,
    id = "aligned_read_length_survival",
    title = "Aligned Read Length Survival",
    caption = "Aligned Read Length Survival"
  )
  report$ggsave(
    "aligned_read_length_survival (Log-scale).png",
    p2,
    width = plotwidth,
    height = plotheight,
    id = "aligned_read_length_survival(log)",
    title = "Aligned Read Length Survival (Log-scale)",
    caption = "Aligned Read Length Survival (Log-scale)"
  )
}

makeAccuracyDensityPlots <- function(report, cd) {
  loginfo("Making Accuracy Density Plots")
  tp = ggplot(cd, aes(x = Accuracy, colour = Condition)) + geom_density(alpha = .5) + plTheme + clScale +
    labs(x = "Accuracy (1 - Mean Errors Per Template Position)", title = "Accuracy by Condition")
  report$ggsave(
    "acc_density.png",
    tp,
    width = plotwidth,
    height = plotheight,
    id = "acc_density",
    title = "Accuracy Distribution",
    caption = "Accuracy Distribution"
  )
  
  loginfo("Making Accuracy Violin Plots")
  tp = ggplot(cd, aes(x = Condition, y = Accuracy, fill = Condition)) + geom_violin() +
    plTheme + clFillScale + themeTilt +
    labs(y = "Accuracy (1 - Mean Errors Per Template Position)", title = "Accuracy by Condition", x = "Condition") +
    geom_boxplot(width = 0.1, fill = "white")
  report$ggsave(
    "acc_violin.png",
    tp,
    width = plotwidth,
    height = plotheight,
    id = "acc_violin",
    title = "Accuracy Violin",
    caption = "Accuracy Violin"
  )
  
  loginfo("Making Template Length vs. Accuracy")
  samps_per_group = sampleSize / length(levels(cd$Condition))
  sample_nigel <-
    function(tbl,
             size,
             replace = FALSE,
             weight = NULL,
             .env = parent.frame())
    {
      #assert_that(is.numeric(size), length(size) == 1, size >= 0)
      weight <- substitute(weight)
      index <- attr(tbl, "indices")
      sizes = sapply(index, function(z)
        min(length(z), size)) # here's my contribution
      sampled <-
        lapply(1:length(index), function(i)
          dplyr:::sample_group(
            index[[i]],
            frac = FALSE,
            tbl = tbl,
            size = sizes[i],
            replace = replace,
            weight = weight,
            .env = .env
          ))
      idx <- unlist(sampled) + 1
      grouped_df(tbl[idx, , drop = FALSE], vars = groups(tbl))
    }
  
  cd2 = cd %>% group_by(Condition) %>% sample_nigel(size = samps_per_group) %>% ungroup()
  tp = ggplot(cd2, aes(x = tlen, y = Accuracy, color = Condition)) + geom_point(alpha = .2) +
    plTheme  + geom_smooth(fill = NA) + clScale + labs(
      y = "Accuracy (1 - Mean Errors Per Template Position)",
      title = paste("Accuracy vs. Template Length\n(Sampled to <= ", sampleSize, ")"),
      x = "Template Length"
    ) + facet_wrap( ~ Condition, nrow = length(levels(cd2$Condition)))
  img_height = min(49.5, 3.6 * length(levels(cd2$Condition)))
  report$ggsave(
    "acc_accvtl.png",
    tp,
    id = "acc_accvtl",
    width = plotwidth,
    height = img_height,
    title = "Template Length v. Accuracy",
    caption = "Template Length v. Accuracy"
  )
  
  tp = ggplot(cd2, aes(x = alen, y = Accuracy, color = Condition)) + geom_point(alpha = .2) +
    plTheme  + geom_smooth(fill = NA) + clScale + labs(
      y = "Accuracy (1 - Mean Errors Per Template Position)",
      title = paste(
        "Accuracy vs. Aligned Read Length (aend - astart)\n(Sampled to <= ",
        sampleSize,
        ")"
      ),
      x = "Aligned Read Length"
    ) + facet_wrap( ~ Condition, nrow = length(levels(cd2$Condition)))
  report$ggsave(
    "acc_accvrl.png",
    tp,
    id = "acc_accvrl",
    width = plotwidth,
    height = img_height,
    title = "Aligned Read Length v. Accuracy",
    caption = "Aligned Read Length v. Accuracy"
  )
  
  tp = ggplot(cd2, aes(x = qrlen, y = alen, color = Condition)) + geom_point(alpha = .2) +
    plTheme  + geom_abline(intercept = 0,
                           slope = 1,
                           color = "red") + clScale + labs(
                             y = "Aligned Read Length (aend - astart)",
                             title = paste("Unaligned Read Length v. Aligned Read Length\n(Sampled to <= ", sampleSize, ")"),
                             x = "Unaligned Read Length (qend - qstart)"
                           ) + facet_wrap( ~ Condition, nrow = length(levels(cd2$Condition)))
  report$ggsave(
    "alen_v_qlen.png",
    tp,
    id = "alen_v_qlen",
    width = plotwidth,
    height = img_height,
    title = "Unaligned Read Length v. Aligned Read Length",
    caption = "Unaligned Read Length v. Aligned Read Length"
  )
  
  loginfo("Making Template Span Violin Plot")
  tp = ggplot(cd, aes(x = Condition, y = tlen, fill = Condition)) + geom_violin() + 
    plTheme + clFillScale + themeTilt + 
    labs(y = "Template Span (tend - tstart)", title = "Template Span Violin Plot", x = "Condition")
  report$ggsave(
    "tlen_violin.png",
    tp,
    id = "tlen_violin",
    width = plotwidth,
    height = plotheight,
    title = "Template Span Violin Plot",
    caption = "Template Span Violin Plot"
  )
  
  loginfo("Making Template Span Density Plot")
  tp = ggplot(cd, aes(x = tlen, colour = Condition)) + geom_density() + 
    plTheme + clScale + themeTilt + 
    labs(y = "Density", title = "Template Span Density Plot", x = "Template Span (tend - tstart)")
  report$ggsave(
    "tlen_density.png",
    tp,
    id = "tlen_density",
    width = plotwidth,
    height = plotheight,
    title = "Template Span Density Plot",
    caption = "Template Span Density Plot"
  )
  
  loginfo("Making Aligned Read Length Density Plot")
  tp = ggplot(cd, aes(x = alen, colour = Condition)) + geom_density() + 
    plTheme + clScale + themeTilt + 
    labs(y = "Density", title = "Aligned Read Length Density Plot", x = "Aligned Read Length (aend - astart)")
  report$ggsave(
    "alen_density.png",
    tp,
    id = "alen_density",
    width = plotwidth,
    height = plotheight,
    title = "Aligned Read Length Density Plot",
    caption = "Aligned Read Length Density Plot"
  )
  
  loginfo("Making Template Span Boxplot")
  tp = ggplot(cd, aes(x = Condition, y = tlen, fill = Condition)) + geom_boxplot() + 
    plTheme + clFillScale + themeTilt + 
    labs(y = "Template Span (tend - tstart)", title = "Template Span Boxplot", x = "Condition")
  report$ggsave(
    "tlen_box.png",
    tp,
    id = "tlen_box",
    width = plotwidth,
    height = plotheight,
    title = "Template Span Boxplot",
    caption = "Template Span Boxplot"
  )
}

makeErateViolinPlots <- function(report, cd) {
  loginfo("Making Error Rate Violin Plots")
  vnames = c("mmrate", "irate", "drate")
  labels = c("Mismatch", "Insertion", "Deletion")
  mkErate <- function(i) {
    vname = vnames[i]
    label = labels[i]
    tp = ggplot(cd, aes_string(x = "Condition", y = vname, fill = "Condition")) + geom_violin() +
      plTheme + clFillScale + geom_boxplot(width = 0.1, fill = "white") + themeTilt +
      labs(
        y = label,
        title = paste(label, " Rate by Condition - Violin Plot"),
        x = "Condition"
      )
    report$ggsave(
      paste("etype_", "_", vname, "_violin.png", sep = ""),
      tp,
      width = plotwidth,
      height = plotheight,
      id = paste("etype_", "_", vname, "_violin", sep = ""),
      title = paste(label, "Rate - Violin Plot"),
      caption = paste(label, "Rate - Violin Plot")
    )
  }
  pv = lapply(1:3, mkErate)
  return(pv)
}

makeErateBoxPlots <- function(report, cd) {
  loginfo("Making Error Rate Boxplots")
  vnames = c("mmrate", "irate", "drate")
  labels = c("Mismatch", "Insertion", "Deletion")
  mkErateBox <- function(i) {
    vname = vnames[i]
    label = labels[i]
    tp = ggplot(cd, aes_string(x = "Condition", y = vname, fill = "Condition")) + geom_boxplot() +
      plTheme + clFillScale + themeTilt + stat_summary(
        fun.y = median,
        colour = "black",
        geom = "text",
        show.legend = FALSE,
        vjust = -0.8,
        aes(label = round(..y.., digits = 4))
      ) +
      labs(
        y = label,
        title = paste(label, " Rate by Condition - Boxplot"),
        x = "Condition"
      )
    report$ggsave(
      paste("etype_", "_", vname, "_boxplot.png", sep = ""),
      tp,
      width = plotwidth,
      height = plotheight,
      id = paste("etype_", "_", vname, "_boxplot", sep = ""),
      title = paste(label, "Rate - Boxplot"),
      caption = paste(label, "Rate - Boxplot")
    )
  }
  pv = lapply(1:3, mkErateBox)
  return(pv)
}

makeBasesDistribution <- function(report, cd) {
  loginfo("Making Bases Data Distribution")
  res = cd %>% group_by(Condition) %>% summarise(Template = sum(tlen), Read = sum(alen)) %>% gather(BaseType, Bases, Template:Read)
  tp = ggplot(res,
              aes(
                x = BaseType,
                y = Bases,
                fill = Condition,
                group = Condition
              )) + geom_bar(stat = "identity", position = "dodge") + 
    plTheme + clFillScale + 
    labs(y = "Total Bases\nsum(end - start)", title = "Total Bases by Condition", x = "Base Type")
  report$ggsave(
    "base_count_bar.png",
    tp,
    width = plotwidth,
    height = plotheight,
    id = "base_count_bar",
    title = "Total Bases",
    caption = "Total Bases"
  )
}

makeYieldHistogram <- function(report, cd) {
  loginfo("Making Yield Histogram")
  tp = ggplot(cd, aes(Condition, fill = Condition)) + geom_bar() + 
    plTheme + themeTilt  + clFillScale + 
    labs(x = "Condition", y = "nReads", title = "nReads by Condition")
  report$ggsave(
    "nreads_hist.png",
    tp,
    width = plotwidth,
    height = plotheight,
    id = "nreads_histogram",
    title = "nReads Histogram",
    caption = "nReads Histogram"
  )
}

# The core function, change the implementation in this to add new features.
makeReport <- function(report) {
  # Make fake data - for debugging
  #condFile = "/pbi/dept/secondary/siv/smrtlink/smrtlink-internal/userdata/jobs-root/005/005578/tasks/pbcommandR.tasks.pbiplot_reseq_condition-0/resolved-tool-contract.json"
  #condjson = jsonlite::fromJSON(condFile)
  #input = condjson$resolved_tool_contract$input_files
  #decoded <- loadReseqConditionsFromPath(input)
  #conds = decoded@conditions
  #tmp = lapply(conds, function(z) data.frame(condition = z@condId,subreadset = z@subreadset, alignmentset = z@alignmentset, referenceset = z@referenceset))
  #conditions = do.call(rbind, tmp)
  
  conditions = report$condition.table
  # Load the pbi index for each data frame
  dfs = lapply(as.character(unique(conditions$MappedSubreads)), function(s) {
    loginfo(paste("Loading alignment set:", s))
    loadPBI2(s)
  })
  # Now combine into one large data frame
  cd = combineConditions(dfs, as.character(conditions$Condition))
  
  ## Let's set the graphic defaults
  n = length(levels(conditions$Condition))
  clFillScale <<- getPBFillScale(n)
  clScale <<- getPBColorScale(n)
  
  cd$tlen = as.numeric(cd$tend - cd$tstart)
  cd$alen = as.numeric(cd$aend - cd$astart)
  cd$errors = as.numeric(cd$mismatches + cd$inserts + cd$dels)
  cd$Accuracy = 1 - cd$errors / cd$tlen
  cd$mmrate = cd$mismatches / cd$tlen
  cd$irate  = cd$inserts / cd$tlen
  cd$drate  = cd$dels / cd$tlen
  cd$qrlen = as.numeric(cd$qend - cd$qstart)
  
  summaries = cd[, .(
    AccuracyRate.Median = median(Accuracy),
    AlnLength = median(tlen),
    ReadLength = median(alen),
    QReadLength = median(qrlen),
    InsertRate = median(irate),
    DeletionRate = median(drate),
    MismatchRate = median(mmrate),
    NumberZMWs = nrow(distinct(cd, hole)),
    NumberAlns = length(hole),
    TotalAlignedBases = sum(tlen),
    TotalReadBases = sum(alen),
    BAMFiles = length(unique(file))
  ),
  by = Condition]
  colnames(summaries) <-
    c(
      "Condition",
      "Accuracy",
      "Aln Length",
      "Aln Read Length (a)",
      "Aln Read Length (q)",
      "Insert Rate",
      "Deletion Rate",
      "Mismatch Rate",
      "# ZMWs",
      "# Alignments",
      "Total Template Bases",
      "Total Read Bases",
      "Total BAM Files"
    )
  report$write.table("sumtable.csv",
                     summaries,
                     id = "sumtable",
                     title = "Summary Statistics (Median Values)")
  
  # Make Plots
  makeReadLengthSurvivalPlots(report, cd)
  makeAccuracyDensityPlots(report, cd)
  makeErateViolinPlots(report, cd)
  makeErateBoxPlots(report, cd)
  makeBasesDistribution(report, cd)
  makeYieldHistogram(report, cd)
  
  # Save the report object for later debugging
  save(report, file = file.path(report$outputDir, "report.Rd"))
  
  # At the end of this function we need to call this last, it outputs the report
  report$write.report()
}


main <- function()
{
    report <- bh2Reporter(
        "condition-table.csv",
        "reports/PbiPlots/report.json",
        "Sampled ZMW metrics")
    makeReport(report)
    0
}

## Leave this as the last line in the file.
logging::basicConfig()
main()
