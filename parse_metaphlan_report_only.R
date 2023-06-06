#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(matrixStats)
library(pals)
library(pheatmap)
library(tidyr)
library(vegan)
library(WriteXLS)

source("utilities.R")

reportFiles = list.files("metaphlan_out/", pattern="\\.txt\\.gz$")

parseReportFile = function(reportFile) {
    report = read.table(paste0("metaphlan_out/", reportFile),
                        sep="\t", header=FALSE, row.names=NULL, quote="")
    colnames(report) = c("clade_name", "NCBI_tax_id",
                         "relative_abundance", "additional_species")
    levels = nchar(gsub("[^|]+", "", report$clade_name)) + 1
    tokens = strsplit(report$clade_name, "|", fixed=TRUE)
    taxa = matrix("", nrow=nrow(report), ncol=max(levels))
    for (i in 1:nrow(report)) {
        for (j in 1:length(tokens[[i]])) {
            taxa[i, j] = gsub("^[kpcofgs]__", "", tokens[[i]][[j]])
        }
    }
    colnames(taxa) = c("k", "p", "c", "o", "f", "g", "s")[1:ncol(taxa)]
    return(data.frame(sample = gsub("^.*/|_profiled_metagenome\\.txt\\.gz", "",
                                    reportFile),
                      key = report$clade_name,
                      pct = report$relative_abundance,
                      ncbi = report$NCBI_tax_id,
                      level = levels,
                      taxa))
}

allData = do.call(rbind, lapply(reportFiles, parseReportFile))

summarizeData = function(level, bacteriaOnly=FALSE) {
    if (bacteriaOnly) {allData = allData[allData$k == "Bacteria", ]}
    out = allData[allData$level == level, c("sample", "key", "pct")]
    out$key = gsub(".*__", "", out$key)
    out = out %>% pivot_wider(names_from=key, values_from=pct) %>% as.data.frame
    out[is.na(out)] = 0
    return(out)
}

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

kingdom = summarizeData(1)
phyla = summarizeData(2)
class = summarizeData(3)
order = summarizeData(4)
family = summarizeData(5)
genus = summarizeData(6)
species = summarizeData(7)

heatAnnot = data.frame(
    row.names = phyla$sample,
    preparation = annotation[gsub("_.*", "", phyla$sample), "Sample Description"]
)
heatAnnot$kit = gsub("Power(Soil)?", "PowerSoil",
                     gsub(" .*", "", heatAnnot$preparation))
heatAnnot$site = gsub(".*(extraction|[kK]it) ", "", heatAnnot$preparation)
heatColors = list(
    kit = c(DNeasy="white", HostZero="black",
            PowerSoil="goldenrod", PureLink="gray"),
    site = c(blank="white", mock="lightgray", "oral swab"="lightseagreen",
             "rectal swab"="orangered", "vaginal swab"="steelblue")
)
