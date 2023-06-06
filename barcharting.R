#!/usr/bin/env Rscript

library(cowplot)
library(dplyr)
library(forcats)
library(ggplot2); theme_set(theme_bw())
library(ggtext)
library(pals)
library(tidyr)

source("utilities.R")

kraken = read.table("kraken2_species_wide.tsv.gz",
                    sep="\t", header=TRUE, row.names=1, check.names=FALSE,
                    comment.char="", quote="")
rownames(kraken) = gsub("_.*", "", rownames(kraken))
colnames(kraken) = gsub("Limosilactobacillus fermentum",
                        "Lactobacillus fermentum",
                        colnames(kraken))

metaph = read.table("metaphlan_species_wide.tsv.gz",
                    sep="\t", header=TRUE, row.names=1, check.names=FALSE,
                    comment.char="", quote="")
rownames(metaph) = gsub("_.*", "", rownames(metaph))
colnames(metaph) = gsub("_", " ", colnames(metaph))
metaph = metaph[rownames(kraken), ]

annot = read.table("SampleUploadTemplateWright.tsv",
                   sep="\t", header=TRUE, row.names=NULL, check.names=FALSE)
annot = annot[!duplicated(annot$SampleName), ]
annot = set_index(annot, "SampleName")
annot = annot[rownames(kraken), ]
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
annot["52O", "Sample Description"] = "HostZero kit vaginal swab"
annot["52V", "Sample Description"] = "HostZero kit oral swab"
## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
annot$kit = annot$"Sample Description" %>%
    gsub(" .*", "", .) %>%
    gsub("Power(Soil)?", "PowerSoil", .)
annot$site = annot$"Sample Description" %>%
    gsub(".*(extraction|[kK]it) ", "", .)

krakenSpecMeans = data.frame(sapply(
    idmap(unique(annot$"Sample Description")),
    function(sd) {colMeans(kraken[annot$"Sample Description" == sd, ])}
), check.names=FALSE)
metaphSpecMeans = data.frame(sapply(
    idmap(unique(annot$"Sample Description")),
    function(sd) {colMeans(metaph[annot$"Sample Description" == sd, ])}
), check.names=FALSE)
nPerSD = 2
topK = Reduce(union, lapply(krakenSpecMeans, function(.) {
    names(sort(structure(., names=rownames(krakenSpecMeans)),
               decreasing=TRUE)[1:(nPerSD+1)])
}))
topM = Reduce(union, lapply(metaphSpecMeans, function(.) {
    names(sort(structure(., names=rownames(metaphSpecMeans)),
               decreasing=TRUE)[1:nPerSD])
}))

allTop = union(topK, topM)
fillVals = structure(rev(nthize(stepped(length(allTop)-1), 2)),
                     names = sort(setdiff(allTop, "Homo sapiens")))
fillVals[["<other identifiable>"]] = "gray60"
fillVals[["Homo sapiens"]] = "gray80"
colVals = structure(ifelse(1:length(fillVals) < length(fillVals) / 2,
                           "black", "white"),
                    names = names(fillVals))
sizeVals = structure(ifelse(colVals == "black", 0.25, 0.0),
                     names = names(colVals))

long = rbind(
    kraken %>%
        reset_index("sample") %>%
        pivot_longer(-sample, names_to="species", values_to="percent") %>%
        mutate(analysis = "kraken"),
    metaph %>%
        reset_index("sample") %>%
        pivot_longer(-sample, names_to="species", values_to="percent") %>%
        mutate(analysis = "metaphlan")
) %>% data.frame
long[ , c("kit", "site")] = annot[long$sample, c("kit", "site")]
abbrev = long
abbrev[!(abbrev$species %in% allTop), "species"] = "<other identifiable>"
abbrev = abbrev %>%
    group_by(sample, species, analysis, kit, site) %>%
    summarize(percent = sum(percent))
abbrev$species = factor(
    as.character(abbrev$species),
    levels = rev(names(fillVals))
)

gg = ggplot(abbrev, aes(x=sample, y=percent,
                        color=species, fill=species, size=species))
gg = gg + facet_grid(analysis ~ site + kit, scales="free_x", space="free")
gg = gg + geom_bar(stat="identity", position="stack")
gg = gg + scale_color_manual(values=colVals)
gg = gg + scale_fill_manual(values=fillVals)
gg = gg + scale_size_manual(values=sizeVals)
gg = gg + ylab("relative abundance (%)")
gg = gg + theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                strip.text.x = element_text(size=5.5))
## png("kraken2_and_metaphlan_barcharts.png",
##     h=960, w=960*3.6, res=144*1.2)
print(gg)
## garbage = dev.off()

ggize = function(letter, title, data, ..., ylab=FALSE, leg=FALSE) {
    data$sample.label = paste0(
            "<span style = 'color: ",
            ifelse(data$sample %in% c("52O", "52V"), "red", "white"),
            ";'>", data$sample, "</span>"
    )
    data$sample.label = fct_reorder(data$sample.label,
                                    as.character(data$sample))
    gg = ggplot(data, ...,
                mapping=aes(x=sample.label, y=percent,
                            color=species, fill=species, size=species))
    gg = gg + facet_grid(analysis ~ kit, scales="free_x", space="free")
    gg = gg + geom_bar(stat="identity", position="stack")
    gg = gg + scale_color_manual(values=colVals)
    gg = gg + scale_fill_manual(values=fillVals)
    gg = gg + scale_size_manual(values=sizeVals)
    gg = gg + ylab("relative abundance (%)")
    gg = gg + theme(axis.text.x = element_markdown(angle = 90),
                    axis.ticks.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    strip.text.x = element_text(angle=90))
    gg = gg + ggtitle(
        latex2exp::TeX(paste0("\\textbf{", letter, ":} ", title))
    )
    gg = gg + ylim(0, 100)
    gg = gg + xlab("")
    if (!ylab) {
        gg = gg + ylab("") + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())
    }
    if (!leg) {gg = gg + theme(legend.position="non",
                               strip.text.y = element_blank(),
                               strip.background.y = element_blank())
    }
    return(gg)
}

ggg = plot_grid(
    ggize("A", "blank", abbrev %>% filter(site=="blank"), ylab=TRUE),
    ggize("B", "mock", abbrev %>% filter(site=="mock")),
    ggize("C", "oral swab", abbrev %>% filter(site=="oral swab")),
    ggize("D", "rectal swab", abbrev %>% filter(site=="rectal swab")),
    ggize("E", "vaginal swab", abbrev %>% filter(site=="vaginal swab"), leg=TRUE),
    nrow = 1,
    rel_widths = c(0.55, 0.65, 0.9, 0.9, 2.35)
)
## png("kraken2_and_metaphlan_barcharts_v2.png",
##     h=960, w=960*3, res=144*1.2)
print(ggg)
## garbage = dev.off()
