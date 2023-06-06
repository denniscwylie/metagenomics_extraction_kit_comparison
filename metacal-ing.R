#!/usr/bin/env Rscript

library(forcats)
library(metacal)

source("parse_metaphlan_report_only.R")
source("metaphlan_merge_lanes.R")

spec = species %>% set_index("sample") %>% t / 100
specMock = spec[ , heatAnnot[colnames(spec), "site"] == "mock"]

rownames(specMock) = rownames(specMock) %>%
    gsub("Pseudomonas_aeruginosa_group", "Pseudomonas_aeruginosa", .) %>%
    gsub("Bacillus_intestinalis", "Bacillus_subtilis", .)

actualMock = readExcel("mock_composition.xlsx")[[1]]
actMock = actualMock %>%
    mutate(taxa = paste0(Genus, "_", Species)) %>%
    (\(.) {matrix(.$`Genomic DNA`,
                  nrow = nrow(.),
                  ncol = ncol(specMock),
                  dimnames = list(.$taxa, colnames(specMock)))}) / 100

offset = 0.001

## =============================================================================
annot = heatAnnot[colnames(specMock), ]
ebs = lapply(sapply(unique(annot$kit), identity), function(.) {
    estimate_bias(specMock[rownames(actMock), annot$kit == .] + offset,
                  actMock[ , annot$kit == .],
                  margin = 2,
                  boot = TRUE)
})

ebsums = lapply(ebs, summary)
ebsums = do.call(rbind, lapply(names(ebsums), function(kit) {
    data.frame(kit=kit, coef(ebsums[[kit]]))
}))

## ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## write.table(ebsums, "metaphlan_metacal_bias_estimates.tsv",
##             sep="\t", row.names=FALSE, quote=FALSE)
## ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ebsums %>%
    group_by(kit) %>%
    summarize(geo_mean = exp(mean(log(estimate))),
              geo_sd = exp(sd(log(estimate))))
## kit       geo_mean geo_sd
## -------------------------
## DNeasy           1   2.60
## HostZero         1   7.62
## PowerSoil        1   1.81
## PureLink         1   2.60

taxLevels = c("Saccharomyces_cerevisiae", "Cryptococcus_neoformans",
              "Pseudomonas_aeruginosa", "Escherichia_coli", "Salmonella_enterica",
              "Staphylococcus_aureus", "Listeria_monocytogenes",
              "Lactobacillus_fermentum", "Enterococcus_faecalis",
              "Bacillus_subtilis") %>% gsub("_", " ", .)

gg = ebsums %>%
    (function(.) {.$taxon = gsub("_", " ", .$taxon); .}) %>%
    ## mutate(taxon = fct_reorder(taxon, estimate)) %>%
    mutate(taxon = factor(as.character(taxon), levels=taxLevels)) %>%
    ggplot(aes(x=taxon, y=estimate, color=kit,
               ymin=estimate/gm_se^2, ymax=estimate*gm_se^2)) +
    geom_hline(yintercept=1, color="gray") +
    geom_pointrange(alpha=0.5, size=0.25) +
    scale_y_log10() +
    scale_color_manual(values=c(DNeasy="steelblue", HostZero="black",
                                PowerSoil="goldenrod", PureLink="gray")) +
    ylab("bias estimate") + xlab("") +
    coord_flip()
## png("metacal_bias_estimates_by_kit_merged_lanes.png",
##     h=960, w=1728, res=288)
print(gg)
## garbage = dev.off()
