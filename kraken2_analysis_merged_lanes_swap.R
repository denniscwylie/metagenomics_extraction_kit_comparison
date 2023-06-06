#!/usr/bin/env Rscript

swap = TRUE

source("parse_kraken2_report_only.R")
source("kraken2_merge_lanes.R")


## =============================================================================
rightSpecs = c(
    "Staphylococcus_aureus",
    "Listeria_monocytogenes",
    "Pseudomonas_aeruginosa",
    "Saccharomyces_cerevisiae",
    "Cryptococcus_neoformans",
    ## "Lactobacillus_fermentum",
    "Limosilactobacillus_fermentum",
    "Enterococcus_faecalis",
    "Salmonella_enterica",
    "Escherichia_coli",
    "Bacillus_subtilis"
)

actualMock = readExcel("mock_composition.xlsx")[[1]]
actualMock$species = paste0(actualMock$Genus, "_", actualMock$Species)
actualMock$percent = actualMock$"Genomic DNA"

spct_ = spct; colnames(spct_) = gsub(" ", "_", colnames(spct_))
isMockRow = heatAnnot$site == "mock"
mockSpecs = names(which(colSums(spct_[isMockRow, 2:ncol(spct_)]) > 0))
mockSpecData = cbind(heatAnnot[isMockRow, ], spct_[isMockRow, mockSpecs]) %>%
    reset_index("sample") %>%
    pivot_longer(c(-sample, -preparation, -kit, -site),
                 names_to="species", values_to="percent") %>%
    as.data.frame
mockSpecData = mockSpecData[mockSpecData$percent > 0, ]
table(mockSpecData[!(mockSpecData$species %in% rightSpecs), c(
                       "preparation", "species"
                   )])
mockSpecData$species = mockSpecData$species %>%
    gsub("Limosilactobacillus_fermentum", "Lactobacillus_fermentum", .)
rightSpecs = rightSpecs %>%
    gsub("Limosilactobacillus_fermentum", "Lactobacillus_fermentum", .)

gg = ggplot(mockSpecData %>% filter(species %in% rightSpecs),
            aes(x=species, y=percent))
gg = gg + geom_boxplot(aes(color=kit), position="identity")
## gg = gg + geom_point(data=actualMock, shape="|", size=10, color="red")
gg = gg + geom_point(aes(color=kit))
gg = gg + scale_y_log10()
gg = gg + scale_color_manual(values=c("steelblue", "black", "goldenrod", "darkgray"))
gg = gg + coord_flip()
print(gg)


## =============================================================================
library(metacal)
offset = 0.001
specMock = mockSpecData[ , c("sample", "species", "percent")] %>%
    filter(species %in% rightSpecs) %>%
    pivot_wider(names_from=species, values_from=percent, values_fill=0) %>%
    set_index("sample") %>%
    t / 100
actMock = actualMock %>%
    mutate(taxa = paste0(Genus, "_", Species)) %>%
    (\(.) {matrix(.$`Genomic DNA`,
                  nrow = nrow(.),
                  ncol = ncol(specMock),
                  dimnames = list(.$taxa, colnames(specMock)))}) / 100
actMock = actMock[rownames(specMock), ]
actMock = sweep(actMock, 2, colSums(actMock), `/`)
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
ebsums %>%
    group_by(kit) %>%
    summarize(geo_mean = exp(mean(log(estimate))),
              geo_sd = exp(sd(log(estimate))))
## kit       geo_mean geo_sd
## -------------------------
## DNeasy           1   3.09
## HostZero         1  11.2 
## PowerSoil        1   3.22
## PureLink         1   3.37

taxLevels = c(#"Saccharomyces_cerevisiae", "Cryptococcus_neoformans",
              "Pseudomonas_aeruginosa", "Escherichia_coli", "Salmonella_enterica",
              "Staphylococcus_aureus", "Listeria_monocytogenes",
              "Lactobacillus_fermentum", "Enterococcus_faecalis",
              "Bacillus_subtilis") %>% gsub("_", " ", .)
gg = ebsums %>%
    (function(.) {.$taxon = gsub("_", " ", .$taxon); .}) %>%
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
print(gg)


## =============================================================================
gg = ggplot(domainData[ , c("sample", "Viruses")] %>%
                (\(.) {cbind(., heatAnnot[.$sample, ])}),
            aes(x=site, y=log10(Viruses+0.005), color=site))
gg = gg + facet_wrap(~ kit, nrow=1)
gg = gg + geom_point() + geom_boxplot(alpha=0)
gg = gg + scale_color_manual(values=heatColors$site %>%
                                    gsub("white", "goldenrod", .) %>%
                                    gsub("lightgray", "darkgray", .) %>%
                                 gsub("steelblue", "dodgerblue", .))
gg = gg + xlab("") + theme(axis.ticks.x = element_blank(),
                           axis.text.x = element_blank())
## png("total_viral_levels_vs_site_by_kit_merged_lanes.png",
##     h=960, w=960*4, res=144*3)
print(gg)
## garbage = dev.off()


## =============================================================================
data.frame(phyla[ , "Chordata", drop=FALSE], heatAnnot[ , c("kit", "site")]) %>%
    group_by(kit, site) %>%
    summarize(med=median(Chordata), mean=mean(Chordata)) %>%
    arrange(site, mean)
## kit       site           med  mean
## ----------------------------------
## PureLink  blank        52.6  52.6 
## HostZero  blank        55.3  55.3 
## PowerSoil blank        66.4  66.4 
## DNeasy    blank        80.7  80.7
## ----------------------------------
## PureLink  mock          8.06  8.07
## HostZero  mock         14.8  16.8 
## DNeasy    mock         21.6  22.9 
## PowerSoil mock         31.3  32.0
## ----------------------------------
## HostZero  oral swab    15.0  16.4 
## PureLink  oral swab    66.5  67.7 
## PowerSoil oral swab    75.9  73.6 
## DNeasy    oral swab    86.8  85.4
## ----------------------------------
## HostZero  rectal swab  48.0  42.4 
## PureLink  rectal swab  54.8  58.7 
## DNeasy    rectal swab  53.1  60.5 
## PowerSoil rectal swab  66.0  66.9
## ----------------------------------
## HostZero  vaginal swab 83.7  74.9 
## PureLink  vaginal swab 96.7  83.3 
## DNeasy    vaginal swab 97.8  87.5 
## PowerSoil vaginal swab 96.1  87.9 

species = summarizeData("S")
data.frame(species[ , "Homo sapiens", drop=FALSE], heatAnnot[ , c("kit", "site")],
           check.names=FALSE) %>%
    group_by(kit, site) %>%
    summarize(med=median(`Homo sapiens`), mean=mean(`Homo sapiens`)) %>%
    arrange(site, mean)
## kit       site           med  mean
## ----------------------------------
## PureLink  blank        52.6  52.6 
## HostZero  blank        55.3  55.3 
## PowerSoil blank        66.4  66.4 
## DNeasy    blank        80.7  80.7
## ----------------------------------
## PureLink  mock          8.06  8.07
## HostZero  mock         14.8  16.8 
## DNeasy    mock         21.6  22.9 
## PowerSoil mock         31.3  32.0
## ----------------------------------
## HostZero  oral swab    15.0  16.4 
## PureLink  oral swab    66.5  67.7 
## PowerSoil oral swab    75.9  73.6 
## DNeasy    oral swab    86.8  85.4
## ----------------------------------
## HostZero  rectal swab  48.0  42.4 
## PureLink  rectal swab  54.8  58.7 
## DNeasy    rectal swab  53.1  60.5 
## PowerSoil rectal swab  66.0  66.9
## ----------------------------------
## HostZero  vaginal swab 83.7  74.9 
## PureLink  vaginal swab 96.7  83.3 
## DNeasy    vaginal swab 97.8  87.5 
## PowerSoil vaginal swab 96.1  87.9 


## =============================================================================
relAbundBar = function(level, n=5, annot=annotation) {
    x = summarizeData(toupper(substr(level, 1, 1)))
    ggd = as.data.frame(pivot_longer(
        x, -sample, names_to=level, values_to="relative abundance"
    ))
    keepers = sort(head(unique(
        ggd[order(ggd$`relative abundance`, decreasing=TRUE), level]
    ), n))
    ggd[!(ggd[[level]] %in% keepers), level] = "other"
    ggd$taxon = factor(ggd[[level]], levels=c(keepers, "other"))
    ggd = ggd %>%
        group_by(sample, taxon) %>%
        summarize(relative_abundance = sum(`relative abundance`)) %>%
        as.data.frame
    colnames(ggd) = gsub("^taxon", level, colnames(ggd))
    ggd$preparation = annot[gsub("_.*", "", ggd$sample), "Sample Description"]
    ggd$kit = gsub("Power(Soil)?", "PowerSoil",
                   gsub(" .*", "", ggd$preparation))
    ggd$site = gsub(".*(extraction|[kK]it) ", "", ggd$preparation)
    gg = ggplot(ggd, aes_string(x="sample", y="relative_abundance", fill=level))
    gg = gg + facet_grid(. ~ site + kit, scales="free_x", space="free")
    gg = gg + geom_bar(stat="identity", position="stack")
    gg = gg + scale_fill_manual(values=c(rev(tol(length(keepers))), "gray"))
    gg = gg + ylab("relative abundance (%)")
    gg = gg + theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    strip.text.x = element_text(size=6))
    return(gg)
}

gg = relAbundBar("species", 12)
print(gg)

## spct[ , "Homo sapiens", drop=FALSE] %>%
##     reset_index("sample") %>%
##     write.table("kraken2_merged_homo_sapiens_percentages.tsv",
##                 sep="\t", row.names=FALSE, quote=FALSE)


## ## -----------------------------------------------------------------------------
## write.table(data.frame(species=rownames(spct), spct, check.names=FALSE),
##             "kraken2_species_wide.tsv",
##             sep="\t", row.names=FALSE, quote=FALSE)
## R.utils::gzip("kraken2_species_wide.tsv", overwrite=TRUE)


## =============================================================================
gg = doNMDS(wisconsin(sqrt(summarizeData("S") %>% set_index("sample"))),
            autotransform=FALSE, trace=TRUE, try=500)
print(gg)
center = function(x) {x - mean(x)}
gg$data %>%
    group_by(kit, site) %>%
    mutate(dist = sqrt(mean(center(MDS1)^2) + mean(center(MDS2)^2))) %>%
    summarize(RMS = round(mean(dist), 3)) %>%
    pivot_wider(names_from=site, values_from=RMS) %>%
    as.data.frame
##       kit blank  mock oral swab rectal swab vaginal swab
## --------------------------------------------------------
##    DNeasy 0.335 0.038     0.316       0.363        0.976
##  HostZero 0.650 0.382     0.113       0.729        1.429
## PowerSoil 0.044 0.035     0.174       0.532        0.875
##  PureLink 0.233 0.062     0.306       0.352        1.010

## ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## write.table(gg$data, "kraken2_nmds_swap.tsv",
##             sep="\t", row.names=FALSE, quote=FALSE)
## ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
