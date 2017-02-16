#!/usr/bin/env Rscript
# This is file is for plotting metrics we can obtain by loading the indexes with optional parsing of the original BAM
# files as well as plots that can be made by taking a sample of those.
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
library(stats)
library(IRanges)

## FIXME: make a real package
myDir = "./scripts/R"
source(file.path(myDir, "Bauhaus2.R"))

# Define a basic addition to all plots
plTheme <- theme_bw(base_size = 14) + theme(plot.title = element_text(hjust = 0.5))
clScale <- scale_colour_brewer(palette = "Set1")
clFillScale <- scale_fill_brewer(palette = "Set1")
themeTilt = theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Custom sampler function to sample min(data, sample) which can't be done with dplyr
### it's a modified copy of sample_n.grouped_df
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

# Subsample the data set to load SNR values
loadSNRforSubset <- function(cd, SNRsampleSize = 5000) {
  cd = cd %>% group_by(file) %>% mutate(framePerSecond = as.numeric(as.character(loadHeader(as.character(file[1]))$readgroups$framerate))) %>% ungroup()
  cd2 = cd %>% group_by(Condition, framePerSecond) %>% sample_nigel(size = SNRsampleSize) %>% ungroup()
  cd2snr = loadExtras(cd2, loadSNR = TRUE)
  cd2 = cbind(cd2, cd2snr)
  cd2
}

makeSamplingPlots <-
  function(report,
           cd,
           conditions,
           sampleSize = 1000) {
    loginfo("Making Sampling plots")
    
    load_alns <- function(tbl) {
      names(tbl)[names(tbl) == "ref"] <- "refName"
      curCondition = tbl$Condition[1]
      rsfname = as.character(conditions$Reference[conditions$Condition == curCondition][1])
      fasta = pbbamr::getReferencePath(rsfname)
      loginfo(paste("Loading fasta file:", fasta))
      alns = loadAlnsFromIndex(tbl, fasta)
      sampleSize = min(nrow(tbl), sampleSize)
      for (i in (1:sampleSize)) {
        alns[[i]]$hole = as.factor(as.character(tbl$hole))[i]
        alns[[i]]$refName = as.factor(as.character(tbl$refName))[i]
      }
      alnsTotal = data.table::rbindlist(alns)
      alnsTotal$Condition = curCondition
      grouped_df(alnsTotal, vars = groups(tbl))
    }
    cd2 = cd %>% group_by(Condition, framePerSecond) %>% sample_nigel(size = sampleSize) %>% do(load_alns(.)) %>% ungroup()
    cd2$AccuBases <- "Inaccurate"
    cd2$AccuBases[cd2$read == cd2$ref] = "Accurate"
    cd2$DC <- cd2$ipd + cd2$pw
    
    # Set a boolean variable to see if all the conditions are internal mode
    # Only when all the conditions are internal mode, the variable is set to TRUE
    internalBAM = TRUE
    if (!("sf" %in% colnames(cd2))) {
      internalBAM = FALSE
    } else {
      cd2internal = cd2 %>% group_by(Condition) %>% summarise(sf = all(unique(sf) %in% NA))
      if (any(cd2internal$sf) == TRUE) {
        internalBAM = FALSE
      }
    }
    
    # For internal mode and non-internal mode produce different dataframes
    # For internal mode, the dataframs contain time information
    if (internalBAM) {
      cd2 = cd2 %>% group_by(Condition, framePerSecond) %>% mutate(
        snrCfac = cut(snrC, breaks = c(0, seq(3, 20), 50)),
        time = cut(sf / framePerSecond, breaks = c(seq(
          0, ceiling(max(cd2$sf) / min(framePerSecond) / 600)
        ) * 600))
      ) %>% ungroup()
      cd3 = cd2[!cd2$read == "-", ] %>% group_by(Condition, framePerSecond, hole, refName) %>% summarise(
        medianpw = median(pw),
        medianipd = median(ipd),
        PolRate = mean(ipd + pw),
        startTime = min(sf) / unique(framePerSecond),
        # Here we group by condition, frame rate and hole, so there should be only one unique frame rate
        endTime = max(sf) / unique(framePerSecond)
      )
    } else {
      # Write a tabel for the missing plots due to non-internal BAM files
      missingPlots = c("Active ZMW - Normalized",
                       "pkMid Box Plot - all reference reads",
                       "pkMid Box Plot - accurate reference reads",
                       "pkMid Box Plot - inaccurate reference reads",
                       "pkMid Violin Plot - all reference reads",
                       "pkMid Violin Plot - accurate reference reads",
                       "pkMid Violin Plot - inaccurate reference reads",
                       "pkMid Density Plot - all reference reads",
                       "pkMid Density Plot - accurate reference reads",
                       "pkMid Density Plot - inaccurate reference reads",
                       "pkMid CDF - all reference reads",
                       "pkMid CDF - accurate reference reads",
                       "pkMid CDF - inaccurate reference reads",
                       "pkMid Histogram - all reference reads",
                       "pkMid Histogram - accurate reference reads",
                       "pkMid Histogram - inaccurate reference reads",
                       "pkMid Density Plot - Accurate vs Inaccurate bases",
                       "PW Trend by Time",
                       "Mean Pulse Width by Time",
                       "Mean Pkmid by Time",
                       "Median Pkmid by Time",
                       "Median Pkmid by Time (Normalized)")
      noninternalBAM = as.data.frame(missingPlots)
      report$write.table("noninternalBAM.csv",
                         noninternalBAM,
                         id = "noninternalBAM",
                         title = "Missing plots that require internal BAM files")
      cd2 = cd2 %>% group_by(Condition, framePerSecond) %>% mutate(snrCfac = cut(snrC, breaks = c(0, seq(3, 20), 50))) %>% ungroup()
      cd3 = cd2[!cd2$read == "-", ] %>% group_by(Condition, framePerSecond, hole, refName) %>% summarise(
        medianpw = median(pw),
        medianipd = median(ipd),
        PolRate = mean(ipd + pw)
      )
    }
    
    cd3$DC <- cd3$medianipd + cd3$medianpw
    cd3$DutyCycle  = cd3$medianpw / (cd3$medianpw + cd3$medianipd)
    
    # Plots based on startFrame (Only produced when "sf" is loaded)
    if (internalBAM) {
      # ActiveZMWs
      activeZMW <- function(rngs) {
        dta <- as.data.frame(rngs)
        dta$time <- 1:dim(dta)[1]
        dtm <-
          melt(as.data.frame(dta),
               id.vars = "time",
               variable.name = "condition")
        dtm
      }
      m <- ceiling(max(cd3$endTime))
      # rngs_unnorm <-
      #   do.call(cbind, tapply(seq.int(1, nrow(cd3)), factor(cd3$Condition), function(idxs) {
      #     x <-
      #       as(coverage(IRanges(cd3$startTime[idxs], cd3$endTime[idxs]), width = m), "vector")
      #   }))
      rngs_norm <-
        do.call(cbind, tapply(seq.int(1, nrow(cd3)), factor(cd3$Condition), function(idxs) {
          x <-
            as(coverage(IRanges(cd3$startTime[idxs], cd3$endTime[idxs]), width = m), "vector")
          x / length(idxs)
        }))
      #    dtm_unnorm = activeZMW(rngs_unnorm)
      dtm_norm = activeZMW(rngs_norm)
      
      # tp = ggplot(dtm_unnorm, aes(x=time,y=value,color=condition,group=condition)) +
      #   geom_line(lty=1,lwd=1) + xlab("Seconds") + ylab("Percentage of Alignments") +
      #   labs(title = "Alignment Percentage by Time: Unnormalized")
      # report$ggsave("active_zmw_unnormalized.png", tp, id = "active_zmw_unnormalized.png", title = "Active ZMW - Unnormalized", caption = "Active ZMW - Unnormalized")
      
      tp = ggplot(dtm_norm,
                  aes(
                    x = time,
                    y = value,
                    color = condition,
                    group = condition
                  )) +
        geom_line(lty = 1, lwd = 1) + xlab("Seconds") + ylab("Percentage of Alignments") +
        labs(title = "Alignment Percentage by Time: Normalized") + plTheme
      report$ggsave(
        "active_zmw_normalized.png",
        tp,
        id = "active_zmw_normalized.png",
        title = "Active ZMW - Normalized",
        caption = "Active ZMW - Normalized"
      )
      
      # pkMid for complete data set, accurate bases, and inaccurate bases
      
      cd2.1 <- cd2[cd2$AccuBases == "Accurate", ]
      cd2.2 <- cd2[cd2$AccuBases == "Inaccurate", ]
      reads <- list(cd2, cd2.1, cd2.2)
      variableTitle <-
        c("all reference reads",
          "accurate reference reads",
          "inaccurate reference reads")
      
      for (i in 1:3) {
        tp = ggplot(reads[[i]], aes(x = Condition, y = pkmid, fill = Condition)) +
          geom_boxplot() + plTheme + themeTilt  + clFillScale +
          facet_wrap( ~ ref)
        report$ggsave(
          paste("pkMid_Box_", variableTitle[i], ".png", sep = ""),
          tp,
          id = paste("pkMid_boxplot_", variableTitle[i], sep = ""),
          title = paste("pkMid Box Plot - ", variableTitle[i], sep = ""),
          caption = paste(
            "Distribution of pkMid for ",
            variableTitle[i],
            " (Boxplot)",
            sep = ""
          )
        )
        
        tp = ggplot(reads[[i]], aes(x = Condition, y = pkmid, fill = Condition)) + geom_violin() +
          geom_boxplot(width = 0.1, fill = "white") + plTheme + themeTilt  + clFillScale +
          facet_wrap( ~ ref)
        report$ggsave(
          paste("pkMid_Violin_", variableTitle[i], ".png", sep = ""),
          tp,
          id = paste("pkMid_violinplot_", variableTitle[i], sep = ""),
          title = paste("pkMid Violin Plot - ", variableTitle[i], sep = ""),
          caption = paste(
            "Distribution of pkMid for ",
            variableTitle[i],
            " (Violin plot)",
            sep = ""
          )
        )
        
        tp = ggplot(reads[[i]], aes(x = pkmid, colour = Condition)) + geom_density(alpha = .5) +
          plTheme + themeTilt  + clScale + facet_wrap( ~ ref) +
          labs(x = "pkMid (after normalization)", title = "pkMid by Condition")
        report$ggsave(
          paste("pkMid_Dens_", variableTitle[i], ".png", sep = ""),
          tp,
          id = paste("pkMid_densityplot_", variableTitle[i], sep = ""),
          title = paste("pkMid Density Plot - ", variableTitle[i], sep = ""),
          caption = paste(
            "Distribution of pkMid for ",
            variableTitle[i],
            " (Density plot)",
            sep = ""
          )
        )
        
        tp = ggplot(reads[[i]], aes(x = pkmid, colour = Condition)) + stat_ecdf() +
          plTheme + themeTilt  + clScale + facet_wrap( ~ ref) +
          labs(x = "pkMid", y = "CDF", title = "pkMid by Condition (CDF)")
        report$ggsave(
          paste("pkMid_CDF_", variableTitle[i], ".png", sep = ""),
          tp,
          id = paste("pkMid_cdf_", variableTitle[i], sep = ""),
          title = paste("pkMid CDF - ", variableTitle[i], sep = ""),
          caption = paste("Distribution of pkMid for ", variableTitle[i],  " (CDF)", sep = "")
        )
        
        tp = ggplot(reads[[i]], aes(x = pkmid, fill = Condition)) + geom_histogram() +
          plTheme + themeTilt  + clFillScale + facet_wrap( ~ ref) +
          labs(x = "pkMid", title = "pkMid by Condition")
        report$ggsave(
          paste("pkMid_Hist_", variableTitle[i], ".png", sep = ""),
          tp,
          id = paste("pkMid_histogram_", variableTitle[i], sep = ""),
          title = paste("pkMid Histogram - ", variableTitle[i], sep = ""),
          caption = paste(
            "Distribution of pkMid for ",
            variableTitle[i],
            " (Histogram)",
            sep = ""
          )
        )
      }
      
      # Density plots to compare pkMid for accurate bases and inaccurate bases
      
      tp = ggplot(cd2, aes(x = pkmid, colour = AccuBases)) + geom_density(alpha = .5) +
        plTheme + themeTilt  + clScale + facet_wrap( ~ Condition + ref, ncol = 5) +
        labs(x = "pkMid", title = "pkMid for accurate bases and inaccurate bases")
      report$ggsave(
        "pkMid_Accu_vs_Inaccu_Dens.png",
        tp,
        id = "pkMid_Accu_Inaccu_densityplot",
        title = "pkMid Density Plot - Accurate vs Inaccurate bases",
        caption = "Distribution of pkMid for Accurate vs inaccurate bases (Density plot)"
      )
      
      # Make Pkmid / PW / PolRate by time plot
      cd2time = cd2 %>% group_by(Condition) %>% mutate(PKMID.Median.Con = median(pkmid)) %>% ungroup() %>% group_by(Condition, time) %>% summarise(
        PW.Mean = mean(pw),
        PKMID.Median = median(pkmid),
        PKMID.Mean = mean(pkmid),
        PolRate = mean(ipd + pw),
        PKMID.Median.Con = median(PKMID.Median.Con)
      )
      cd2time$time = as.numeric(cd2time$time) * 10
      
      tp = ggplot(cd2time,
                  aes(
                    x = time,
                    y = PW.Mean,
                    color = Condition,
                    group = Condition
                  )) + geom_point() +
        geom_line()  + clScale + plTheme + themeTilt + labs(y = "Mean Pulse Width",
                                                            x = "Minute",
                                                            title = "PW Trend by Time")
      report$ggsave(
        "pw_mean_by_time.png",
        tp,
        id = "pw_mean_by_time",
        title = "Mean Pulse Width by Time",
        caption = "Mean Pulse Width by Time"
      )
      
      tp = ggplot(cd2time,
                  aes(
                    x = time,
                    y = PKMID.Mean,
                    color = Condition,
                    group = Condition
                  )) + geom_point() +
        geom_line()  + clScale + plTheme + themeTilt + labs(y = "Mean Pkmid",
                                                            x = "Minute",
                                                            title = "Pkmid Trend by Time")
      report$ggsave(
        "pkmid_mean_by_time.png",
        tp,
        id = "pkmid_mean_by_time",
        title = "Mean Pkmid by Time",
        caption = "Mean Pkmid by Time"
      )
      
      tp = ggplot(cd2time,
                  aes(
                    x = time,
                    y = PKMID.Median,
                    color = Condition,
                    group = Condition
                  )) + geom_point() +
        geom_line()  + clScale + plTheme + themeTilt + labs(y = "Median Pkmid",
                                                            x = "Minute",
                                                            title = "Pkmid Trend by Time")
      report$ggsave(
        "pkmid_median_by_time.png",
        tp,
        id = "pkmid_median_by_time",
        title = "Median Pkmid by Time",
        caption = "Median Pkmid by Time"
      )
      
      tp = ggplot(cd2time,
                  aes(
                    x = time,
                    y = PKMID.Median / PKMID.Median.Con,
                    color = Condition,
                    group = Condition
                  )) + geom_point() +
        geom_line()  + clScale + plTheme + themeTilt + labs(y = "Median Pkmid (Normalized)",
                                                            x = "Minute",
                                                            title = "Pkmid Trend by Time (Normalized)")
      report$ggsave(
        "pkmid_median_by_time_normalized.png",
        tp,
        id = "pkmid_median_by_time_normalized",
        title = "Median Pkmid by Time (Normalized)",
        caption = "Median Pkmid by Time (Normalized)"
      )
    }
    
    # Polymerization Rate by Reference
    tp = ggplot(cd3, aes(x = refName, y = PolRate, fill = Condition)) + geom_boxplot(position = "dodge") + 
      plTheme + themeTilt  + clFillScale + 
      labs(x = "Reference", y = "Polymerization Rate", title = "Polymerization Rate by Reference")
    
    report$ggsave(
      "polrate_ref_box.png",
      tp,
      id = "polrate_ref_box",
      title = "Polymerization Rate by Reference",
      caption = "Polymerization Rate by Reference"
    )
    
    # PW by Template Base
    tp = ggplot(cd2, aes(x = pw, colour = Condition)) + geom_density(alpha = .5) + xlim(0, 50) +
      labs(
        y = "frequency",
        x = "frames",
        title = paste("Pulse Width\n(From ", sampleSize, "Sampled Alignments)")
      ) + plTheme + themeTilt + clScale + facet_wrap( ~ ref)
    report$ggsave(
      "pw_by_template.png",
      tp,
      id = "pw_by_template.png",
      title = "Pulse Width by Template Base",
      caption = "Pulse Width by Template Base"
    )
    
    tp = ggplot(cd2, aes(x = pw, colour = Condition)) + stat_ecdf() + xlim(0, 50) +
      labs(
        y = "CDF",
        x = "frames",
        title = paste("Pulse Width_CDF\n(From ", sampleSize, "Sampled Alignments)")
      ) + plTheme + themeTilt + clScale + facet_wrap( ~ ref)
    report$ggsave(
      "pw_by_template_cdf.png",
      tp,
      id = "pw_by_template_cdf.png",
      title = "Pulse Width by Template Base (CDF)",
      caption = "Pulse Width by Template Base (CDF)"
    )
    
    # # Local Polymerization Rate
    # # Note that thsi plot only works for the unrolled data set (one read per ZMW)
    # if (length(cd$hole) = length(unique(cd$hole))) {
    #   cd2Unrolled = cd2 %>% group_by(hole, Condition) %>% mutate(UnrolledTemplateLocation = seq_len(n())) %>% ungroup() %>% mutate(UnrolledTemplateGroup = cut(UnrolledTemplateLocation, breaks = c((seq(0, 1000) * 50), max(UnrolledTemplateLocation)))) %>% group_by(UnrolledTemplateGroup, Condition) %>% summarise(mdPolRate = median(1/(ipd + pw)))
    #   pd <- position_dodge(0.2) # move them .05 to the left and right
    #   tp = ggplot(cd2Unrolled, aes(x = UnrolledTemplateGroup, y = mdPolRate, color = Condition, group = Condition)) + geom_point() + geom_line() +
    #     # geom_errorbar(aes(ymin = mdPolRate - madPolRate, ymax = mdPolRate + madPolRate), width = .1, position = pd) +
    #     plTheme + themeTilt + clScale + labs(x = "Position Bin of Unrolled Template Span", y = "Median Polymerization Rate (50bp bins)")
    #   report$ggsave("localpolrate.png", tp,
    #                 id = "local_polrate",
    #                 title = "Local Polymerization Rate", caption = "Local Polymerization Rate")
    # }
    
    # IPD Plots
    maxIPD = 100
    tp = ggplot(cd2[cd2$ipd < maxIPD,], aes(x = Condition, y = ipd, fill = Condition)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white") +
      labs(
        y = paste("IPD (Truncated < ", maxIPD, ")", sep = ""),
        title = paste("IPD Distribution\n(From ", sampleSize, "Sampled Alignments)")
      ) +
      plTheme + themeTilt + clFillScale
    report$ggsave("ipddist.png",
                  tp,
                  id = "ipd_violin",
                  title = "IPD Distribution - Violin Plot",
                  caption = "IPD Distribution - Violin Plot")
    
    tp = ggplot(cd2[cd2$ipd < maxIPD,], aes(x = Condition, y = ipd, fill = Condition)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white") +
      labs(
        y = paste("IPD (Truncated < ", maxIPD, ")", sep = ""),
        title = paste("IPD Distribution\n(From ", sampleSize, "Sampled Alignments)")
      ) +
      plTheme + themeTilt + clFillScale + facet_wrap( ~ ref)
    report$ggsave(
      "ipddistbybase_violin.png",
      tp,
      id = "ipd_violin_by_base",
      title = "IPD Distribution by Ref Base - Violin Plot",
      caption = "IPD Distribution by Ref Base - Violin Plot"
    )
    
    tp = ggplot(cd2[cd2$ipd < maxIPD,], aes(x = Condition, y = ipd, fill = Condition)) + geom_boxplot() +
      labs(
        y = paste("IPD (Truncated < ", maxIPD, ")", sep = ""),
        title = paste("IPD Distribution\n(From ", sampleSize, "Sampled Alignments)")
      ) +
      stat_summary(
        fun.y = median,
        colour = "black",
        geom = "text",
        show.legend = FALSE,
        vjust = -0.8,
        aes(label = round(..y.., digits = 3))
      ) +
      plTheme + themeTilt + clFillScale + facet_wrap( ~ ref)
    report$ggsave(
      "ipddistbybase_boxplot.png",
      tp,
      id = "ipd_boxplot_by_base",
      title = "IPD Distribution by Ref Base - Boxplot",
      caption = "IPD Distribution by Ref Base - Boxplot"
    )
    
    # PW Plots
    maxPW = 20
    cd2$Insertion = cd2$ref == "-"
    tp = ggplot(cd2[cd2$pw < maxPW,], aes(x = Condition, y = pw, fill = Insertion)) + geom_violin() +
      labs(
        y = paste("PW (Truncated < ", maxPW, ")", sep = ""),
        title = paste("PW Distribution\n(From ", sampleSize, "Sampled Alignments)")
      ) +
      stat_summary(
        fun.y = median,
        colour = "black",
        geom = "text",
        show.legend = FALSE,
        vjust = -0.8,
        aes(label = round(..y.., digits = 3))
      ) +
      plTheme + themeTilt + clFillScale
    report$ggsave(
      "pw_violin.png",
      tp,
      id = "pw_violin",
      title = "PW Distribution - Violin Plot",
      caption = "PW Distribution - Violin Plot"
    )
    
    tp = ggplot(cd2[cd2$pw < maxPW,], aes(x = Condition, y = pw, fill = Insertion)) + geom_boxplot() +
      labs(
        y = paste("PW (Truncated < ", maxPW, ")", sep = ""),
        title = paste("PW Distribution\n(From ", sampleSize, "Sampled Alignments)")
      ) +
      plTheme + themeTilt + clFillScale
    report$ggsave(
      "pw_boxplot.png",
      tp,
      id = "pw_boxplot",
      title = "PW Distribution - Boxplot",
      caption = "PW Distribution - Boxplot"
    )
    
    tp = tp + facet_wrap( ~ ref)
    report$ggsave(
      "pw_boxplot_by_base.png",
      tp,
      id = "pw_boxplot_by_base",
      title = "PW Distribution By Base",
      caption = "PW Distribution"
    )
    
    # Make a median PW plot
    summaries = cd2[cd2$ipd < maxIPD,] %>% group_by(Condition, ref) %>% summarise(PW.Median = median(pw), IPD.Median = median(ipd)) %>% ungroup()
    
    report$write.table("medianIPD.csv",
                       data.frame(summaries),
                       id = "medianIPD",
                       title = "Median IPD/PW Values by Reference")
    
    # Duty Cycle plot
    tp = ggplot(cd3, aes(x = Condition, y = DutyCycle, fill = Condition)) + geom_boxplot() +
      labs(
        y = "median(PW)/(median(PW) + median(IPD))",
        title = paste("Duty Cycle\n(From ", sampleSize, "Sampled Alignments)")
      ) +
      stat_summary(
        fun.y = median,
        colour = "black",
        geom = "text",
        show.legend = FALSE,
        vjust = -0.8,
        aes(label = round(..y.., digits = 3))
      ) + plTheme + themeTilt + clFillScale
    report$ggsave(
      "dutycycle_boxplot.png",
      tp,
      id = "dutycycle_boxplot",
      title = "Duty Cycle - Boxplot",
      caption = "Duty Cycle - Boxplot"
    )
    
    # Local PolRate plot
    tp = ggplot(cd3, aes(x = Condition, y = 1 / DC, fill = Condition)) + geom_boxplot() +
      labs(
        y = "1/(median(PW) + median(IPD))",
        title = paste(
          "Local Ploymerization Rate\n(From ",
          sampleSize,
          "Sampled Alignments)"
        )
      ) +
      stat_summary(
        fun.y = median,
        colour = "black",
        geom = "text",
        show.legend = FALSE,
        vjust = -0.8,
        aes(label = round(..y.., digits = 3))
      ) + plTheme + themeTilt + clFillScale
    report$ggsave(
      "localpolrate_boxplot.png",
      tp,
      id = "localpolrate_boxplot",
      title = "Local PolRate - Boxplot",
      caption = "Local PolRate - Boxplot"
    )
    
    # Global/Local PolRate plot
    tp = ggplot(cd3, aes(
      x = PolRate * (medianpw + medianipd) / log(2) / 6400,
      colour = Condition
    )) + geom_density(alpha = .5) + xlim(0, 2) +
      labs(
        y = "density",
        x = "PolRate*(median(PW) + median(IPD))/ln(2)",
        title = paste(
          "John Eid's Global/Local Ploymerization Rate\n(From ",
          sampleSize,
          "Sampled Alignments)"
        )
      ) + plTheme + themeTilt + clScale
    report$ggsave(
      "global_localpolrate.png",
      tp,
      id = "global_localpolrate",
      title = "Global/Local PolRate",
      caption = "Global/Local PolRate"
    )
    
    # Now mismatch insertions
    errorRates = cd2[cd2$read != "-",] %>%
      group_by(Condition, snrCfac, read) %>%
      summarise(correct = sum(read == ref),
                incorrect = sum(read != ref)) %>%
      mutate(erate = incorrect / (correct + incorrect)) %>%
      ungroup()
    tp = ggplot(errorRates,
                aes(
                  x = snrCfac,
                  y = erate,
                  color = Condition,
                  group = Condition
                )) + geom_point() +
      geom_line()  + clScale + plTheme + themeTilt + labs(y = "Error Rate (per called BP)\nFrom Sampled Alignments",
                                                          x = "SNR C Bin",
                                                          title = "Error Rates By Called Base") +  facet_wrap(~ read)
    report$ggsave(
      "bperr_rate_by_snr.png",
      tp,
      id = "bp_err_rate_by_snr",
      title = "BP Error Rates by SNR",
      caption = "BP Error Rates by SNR"
    )
    
    # Now for mismatch rates
    mmRates = cd2[cd2$read != "-" & cd2$ref != "-",] %>%
      group_by(Condition, snrCfac, ref) %>%
      summarise(correct = sum(read == ref),
                incorrect = sum(read != ref)) %>%
      mutate(erate = incorrect / (correct + incorrect)) %>%
      ungroup()
    tp = ggplot(mmRates,
                aes(
                  x = snrCfac,
                  y = erate,
                  color = Condition,
                  group = Condition
                )) + geom_point() +
      geom_line()  + clScale + plTheme + themeTilt + labs(y = "Mismatch Rate (per ref BP)\nFrom Sampled Alignments",
                                                          x = "SNR C Bin",
                                                          title = "Mismatch Rates By Template Base") +  facet_wrap(~ ref)
    report$ggsave(
      "bpmm_rate_by_snr.png",
      tp,
      id = "bp_mm_err_rate_by_snr",
      title = "Mismatch Rates by SNR",
      caption = "Mismatch Rates by SNR"
    )
    
    # Table of the polymerization rate
    
    pr <- aggregate(DC ~ ref + Condition, cd2, median)
    pr.rs <-
      reshape(pr,
              idvar = 'Condition',
              timevar = 'ref',
              direction = 'wide')
    report$write.table("medianPolymerizationRate.csv",
                       pr.rs,
                       id = "medianPolymerizationRate",
                       title = "Median Polymerization Rate")
  }


makeErrorsBySNRPlots <- function(report, cd, conLevel = 0.95) {
  # Rearrange the data into SNR bins and get summaries
  CI = conLevel / 2 + 0.5
  cd$alnLength = as.numeric(cd$tend - cd$tstart)
  cd2 = cd %>% dplyr::mutate(snrCfac = cut(snrC, breaks = c(0, seq(3, 20), 50))) %>%
    dplyr::mutate(
      mmrate = mismatches / alnLength,
      insrate = inserts / alnLength,
      delrate = dels / alnLength,
      acc = 1 - (mismatches + inserts + dels) / alnLength
    ) %>%
    dplyr::group_by(Condition, snrCfac) %>%
    dplyr::summarise(
      mmratemean = mean(mmrate),
      mmrateci = sd(mmrate) / sqrt(n()) * qt(CI, n() - 1),
      insratemean = mean(insrate),
      insrateci = sd(insrate) / sqrt(n()) * qt(CI, n() - 1),
      delratemean = mean(delrate),
      delrateci = sd(delrate) / sqrt(n()) * qt(CI, n() - 1),
      accmean = mean(acc),
      accci = sd(acc) / sqrt(n()) * qt(CI, n() - 1),
      insdelratmean = sum(inserts) / sum(dels),
      insdelratci = sqrt(1 / (mean(insrate) ^ 2) * ((mean(delrate) / mean(insrate)) ^
                                                      2 * sd(insrate) ^ 2
                                                    + sd(delrate) ^ 2 - 2 * mean(delrate) /
                                                      mean(insrate)
                                                    * cov(insrate, delrate)
      )) / n() * qt(CI, n() - 1)
    ) %>%
    dplyr::ungroup()
  
  # Add error bars
  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- position_dodge(0.2) # move them .05 to the left and right
  
  tp = ggplot(cd2,
              aes(
                x = snrCfac,
                y = accmean,
                color = Condition,
                group = Condition
              )) + geom_point() + geom_line() +
    geom_errorbar(
      aes(ymin = accmean - accci, ymax = accmean + accci),
      width = .1,
      position = pd
    ) +
    plTheme + themeTilt + clScale + labs(x = "SNR C Bin", y = "Accuracy (1 - errors per template pos)")
  report$ggsave(
    "snrvsacc.png",
    tp,
    id = "snr_vs_acc",
    title = "SNR vs Accuracy",
    caption = "SNR vs. Accuracy"
  )
  
  tp = ggplot(cd2,
              aes(
                x = snrCfac,
                y = insratemean,
                color = Condition,
                group = Condition
              )) + geom_point() + geom_line() +
    geom_errorbar(
      aes(ymin = insratemean - insrateci, ymax = insratemean + insrateci),
      width = .1,
      position = pd
    ) +
    plTheme + themeTilt + clScale + labs(x = "SNR C Bin", y = "Insertion Rate")
  report$ggsave(
    "snrvsinsertion.png",
    tp,
    id = "snr_vs_ins",
    title = "SNR vs Insertion Rate",
    caption = "SNR vs. Insertion Rate"
  )
  
  tp = ggplot(cd2,
              aes(
                x = snrCfac,
                y = delratemean,
                color = Condition,
                group = Condition
              )) + geom_point() + geom_line() +
    geom_errorbar(
      aes(ymin = delratemean - delrateci, ymax = delratemean + delrateci),
      width = .1,
      position = pd
    ) +
    plTheme + themeTilt + clScale + labs(x = "SNR C Bin", y = "Deletion Rate")
  report$ggsave(
    "snrvsdeletion.png",
    tp,
    id = "snr_vs_del",
    title = "SNR vs Deletion Rate",
    caption = "SNR vs. Deletion Rate"
  )
  
  tp = ggplot(cd2,
              aes(
                x = snrCfac,
                y = mmratemean,
                color = Condition,
                group = Condition
              )) + geom_point() + geom_line() +
    geom_errorbar(
      aes(ymin = mmratemean - mmrateci, ymax = mmratemean + mmrateci),
      width = .1,
      position = pd
    ) +
    plTheme + themeTilt + clScale + labs(x = "SNR C Bin", y = "Mismatch Rate")
  report$ggsave(
    "snrvsmismatch.png",
    tp,
    id = "snr_vs_mm",
    title = "SNR vs Mismatch Rate",
    caption = "SNR vs. Mismatch Rate"
  )
  
  tp = ggplot(cd2,
              aes(
                x = snrCfac,
                y = insdelratmean,
                color = Condition,
                group = Condition
              )) + geom_point() + geom_line() +
    geom_errorbar(
      aes(ymin = insdelratmean - insdelratci, ymax = insdelratmean + insdelratci),
      width = .1,
      position = pd
    ) +
    plTheme + themeTilt + clScale + labs(x = "SNR C Bin", y = "Insertion Rate / Deletion Rate")
  report$ggsave(
    "snrvsindelrat.png",
    tp,
    id = "snr_vs_indel_rat",
    title = "SNR vs Relative Indels",
    caption = "SNR vs. Indel Rate / Deletion Rate"
  )
}

# The core function, change the implementation in this to add new features.
makeReport <- function(report) {

  # Let's load all the conditions with SNR data
  conditions = report$condition.table
  # Load the pbi index for each data frame
  dfs = lapply(as.character(unique(conditions$MappedSubreads)), function(s) {
    loginfo(paste("Loading alignment set:", s))
    loadPBI2(s)
  })
  # Now combine into one large data frame
  ##browser()
  cd = combineConditions(dfs, as.character(conditions$Condition))
  
  # Subsample cd to get a smaller data frame, load SNR values for this dataframe
  cd = loadSNRforSubset(cd)
  cd = as.data.table(cd)
  
  ## Let's set the graphic defaults
  n = length(levels(conditions$Condition))
  clFillScale <<- getPBFillScale(n)
  clScale <<- getPBColorScale(n)
  
  
  # Let's look at SNR distributions
  logging::loginfo("Making SNR Distribution Plots")
  snrs = cd[, .(Condition, hole, snrA, snrC, snrG, snrT)]
  colnames(snrs) = sub("snr", "", colnames(snrs))
  snrs = snrs %>% gather(channel, SNR, A, C, G, T)
  tp = ggplot(snrs, aes(x = Condition, y = SNR, fill = Condition)) + geom_violin() +
    geom_boxplot(width = 0.1, fill = "white") + plTheme + themeTilt  + clFillScale +
    facet_wrap(~ channel)
  report$ggsave(
    "snrViolin.png",
    tp,
    id = "snr_violin",
    title = "SNR Violin Plot",
    caption = "Distribution of SNR in Aligned Files (Violin plot)"
  )
  
  tp = ggplot(snrs, aes(x = SNR, colour = Condition)) + geom_density(alpha = .5) +
    plTheme + themeTilt  + clScale + facet_wrap(~ channel) +
    labs(x = "SNR", title = "Distribution of SNR in Aligned Files (Density plot)")
  report$ggsave(
    "snrDensity.png",
    tp,
    id = "snr_density",
    title = "SNR Density Plot",
    caption = "Distribution of SNR in Aligned Files (Density plot)"
  )
  
  tp = ggplot(snrs, aes(x = Condition, y = SNR, fill = Condition)) +
    geom_boxplot() + plTheme + themeTilt  + clFillScale +
    facet_wrap( ~ channel)
  report$ggsave(
    "snrBoxNoViolin.png",
    tp,
    id = "snr_boxplot",
    title = "SNR Box Plot",
    caption = "Distribution of SNR in Aligned Files (Boxplot)"
  )
  
  snrs = NULL # make available for GC
  tp = NULL
  
  # Get Errors by SNR plot
  makeErrorsBySNRPlots(report, cd)
  
  # Now plots from sampling alignments
  makeSamplingPlots(report, cd, conditions, sampleSize = 1000)
  
  # Make a median SNR table
  summaries = cd[, .(
    A.Median = median(snrA),
    C.Median = median(snrC),
    G.Median = median(snrG),
    T.Median = median(snrT)
  ),  by = Condition]
  report$write.table("medianSNR.csv",
                     summaries,
                     id = "medianSNR",
                     title = "Median SNR values")
  
  # Save the report object for later debugging
  save(report, file = file.path(report$outputDir, "report.Rd"))
  
  # Output error rates by SNR
  loginfo("Examining error rates by SNR Bin")
  # At the end of this function we need to call this last, it outputs the report
  report$write.report()
}

main <- function()
{
    report <- bh2Reporter(
        "condition-table.csv",
        "reports/PbiSampledPlots/report.json",
        "Sampled ZMW metrics")
    makeReport(report)
    0
}

## Leave this as the last line in the file.
logging::basicConfig()
main()
