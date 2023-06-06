#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(pals)
library(pheatmap)
library(tidyr)
library(vegan)

source("utilities.R")

reportFiles = list.files("kraken2_results/", pattern="_report")

parseReportFile = function(reportFile) {
    report = read.table(paste0("kraken2_results/", reportFile),
                        sep="\t", header=FALSE, row.names=NULL, quote="")
    levels = (nchar(gsub("\\S.*", "", report[ , 6])) / 2) + 1
    tokens = gsub("^\\s+", "", report[ , 6])
    taxa = matrix("", nrow=nrow(report), ncol=max(levels))
    activeTaxa = as.character(rep(NA, max(levels)))
    activeLevel = 1
    for (i in 1:nrow(report)) {
        if (levels[[i]] < activeLevel) {
            activeTaxa[(levels[[i]]+1):length(activeTaxa)] = NA
        }
        activeLevel = levels[[i]]
        activeTaxa[activeLevel] = tokens[[i]]
        taxa[i, ] = activeTaxa
    }
    colnames(taxa) = gsub("_(\\d)$", "_0\\1",
                          paste0("taxonomic_level_", 1:ncol(taxa)))
    return(data.frame(sample = gsub("^.*/|_report\\.tsv.*", "", reportFile),
                      key = report[ , 6],
                      pct = report[ , 1],
                      count = report[ , 2],
                      res = report[ , 4],
                      taxa))
}

allData = do.call(rbind, lapply(reportFiles, parseReportFile))

annotation = read.table("SampleUploadTemplateWright.tsv",
                        sep="\t", header=TRUE, row.names=NULL, check.names=FALSE)
annotation = annotation[!duplicated(annotation$SampleName), ]
annotation = set_index(annotation, "SampleName")

if (exists("swap") && swap) {
    annotation["52O", "Sample Description"] = "HostZero kit vaginal swab"
    annotation["52V", "Sample Description"] = "HostZero kit oral swab"
} else if (exists("exclude") && exclude) {
    annotation = annotation[!(rownames(annotation) %in% c("52O", "52V")), ]
    allData = allData[!grepl("^52[OV]_S", allData$sample), ]
}

summarizeData = function(level) {
    out = allData[allData$res == level, c("sample", "key", "pct")]
    out$key = gsub("^ +", "", out$key)
    out = out %>% pivot_wider(names_from=key, values_from=pct) %>% as.data.frame
    out[is.na(out)] = 0
    return(out)
}

ppct = summarizeData("P") %>% set_index("sample")
gpct = summarizeData("G") %>% set_index("sample")
spct = summarizeData("S") %>% set_index("sample")

## =============================================================================
summarizeCounts = function(level) {
    out = allData[allData$res == level, c("sample", "key", "count")]
    out$key = gsub("^ +", "", out$key)
    out = out %>% pivot_wider(names_from=key, values_from=count) %>% as.data.frame
    out[is.na(out)] = 0
    return(out)
}

pc = summarizeCounts("P") %>% set_index("sample")
gc = summarizeCounts("G") %>% set_index("sample")
sc = summarizeCounts("S") %>% set_index("sample")
