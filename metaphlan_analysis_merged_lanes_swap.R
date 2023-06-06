#!/usr/bin/env Rscript

library(lmerTest)

swap = TRUE

source("parse_metaphlan_report_only.R")
source("metaphlan_merge_lanes.R")


## =============================================================================
gg = doNMDS(wisconsin(sqrt(species %>% set_index("sample"))),
            autotransform=FALSE, trace=TRUE, try=500)
## png("bray_curtis_nmds_metaphlan_species_merged_lanes_swap.png",
##     h=960, w=960*1.3, res=144*1.2)
print(gg)
## garbage = dev.off()

## ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## write.table(gg$data, "metaphlan_nmds_swap.tsv",
##             sep="\t", row.names=FALSE, quote=FALSE)
## ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

center = function(x) {x - mean(x)}
gg$data %>%
    group_by(kit, site) %>%
    mutate(dist = sqrt(mean(center(MDS1)^2) + mean(center(MDS2)^2))) %>%
    summarize(RMS = round(mean(dist), 3)) %>%
    pivot_wider(names_from=site, values_from=RMS) %>%
    as.data.frame
## --------------------------------------------------------
##       kit blank  mock oral swab rectal swab vaginal swab
##    DNeasy 0.059 0.003     0.560       0.622        1.103
##  HostZero 0.000 0.613     0.492       1.203        1.907
## PowerSoil 0.000 0.031     0.401       0.526        0.913
##  PureLink 0.439 0.176     0.465       0.460        1.039


## =============================================================================
rightSpecs = c(
    "Staphylococcus_aureus",
    "Listeria_monocytogenes",
    ## "Pseudomonas_aeruginosa",
    "Pseudomonas_aeruginosa_group",
    "Saccharomyces_cerevisiae",
    "Cryptococcus_neoformans",
    "Lactobacillus_fermentum",
    "Enterococcus_faecalis",
    "Salmonella_enterica",
    "Escherichia_coli",
    ## "Bacillus_subtilis"
    "Bacillus_intestinalis"
)

actualMock = readExcel("mock_composition.xlsx")[[1]]
actualMock$species = paste0(actualMock$Genus, "_", actualMock$Species)
actualMock$percent = actualMock$"Genomic DNA"

isMockRow = heatAnnot$site == "mock"
mockSpecs = names(which(colSums(species[isMockRow, 2:ncol(species)]) > 0))
mockSpecData = cbind(heatAnnot[isMockRow, ], species[isMockRow, mockSpecs]) %>%
    reset_index("sample") %>%
    pivot_longer(c(-sample, -preparation, -kit, -site),
                 names_to="species", values_to="percent") %>%
    as.data.frame
mockSpecData = mockSpecData[mockSpecData$percent > 0, ]
table(mockSpecData[!(mockSpecData$species %in% rightSpecs), c(
                       "preparation", "species"
                   )])
mockSpecData$species = mockSpecData$species %>%
    gsub("Bacillus_intestinalis", "Bacillus_subtilis", .) %>%
    gsub("Pseudomonas_aeruginosa_group", "Pseudomonas_aeruginosa", .)
rightSpecs = rightSpecs %>%
    gsub("Bacillus_intestinalis", "Bacillus_subtilis", .) %>%
    gsub("Pseudomonas_aeruginosa_group", "Pseudomonas_aeruginosa", .)    

mockSpecData$species = factor(
    as.character(mockSpecData$species),
    levels = c("Saccharomyces_cerevisiae", "Cryptococcus_neoformans",
               "Pseudomonas_aeruginosa", "Escherichia_coli", "Salmonella_enterica",
               "Staphylococcus_aureus", "Listeria_monocytogenes",
               "Lactobacillus_fermentum", "Enterococcus_faecalis",
               "Bacillus_subtilis",
               "Sphingomonas_melonis", "Burkholderia_contaminans")
)

gg = ggplot(mockSpecData %>% filter(species %in% rightSpecs),
            aes(x=species, y=percent))
gg = gg + geom_boxplot(aes(color=kit), position="identity")
## gg = gg + geom_point(data=actualMock, shape="|", size=10, color="red")
gg = gg + geom_point(aes(color=kit))
gg = gg + scale_y_log10()
gg = gg + scale_color_manual(values=c("steelblue", "black", "goldenrod", "darkgray"))
gg = gg + coord_flip()
## png("main_mock_figure_swap.png", h=1920*0.5, w=1920, res=288*0.75)
print(gg)
## garbage = dev.off()


## =============================================================================
fungalPhyla = c(
    "Ascomycota", "Basidiomycota"## ,
    ## "Chytridiomycota", "Zygomycota", "Glomeromycota"
)
heatAnnot[rowSums(phyla[ , fungalPhyla]) > 0, ]
##                   preparation       kit site
## M1_S25 PowerSoil Pro Kit mock PowerSoil mock
## M2_S26 PowerSoil Pro Kit mock PowerSoil mock
## M3_S27 PowerSoil Pro Kit mock PowerSoil mock
## M4_S61      HostZero kit mock  HostZero mock
## M5_S62      HostZero kit mock  HostZero mock
## M6_S63      HostZero kit mock  HostZero mock
## M7_S41      PureLink kit mock  PureLink mock
## M8_S42      PureLink kit mock  PureLink mock
## M9_S43      PureLink kit mock  PureLink mock


## =============================================================================
krakenHS = read.table("kraken2_merged_homo_sapiens_percentages.tsv",
                      sep="\t", header=TRUE, row.names=1, check.names=FALSE)
allBact = structure(kingdom$Bacteria, names=kingdom$sample)
nonBact = 100 - allBact
specedBact = allData %>%
    filter(level == 7 & k == "Bacteria") %>%
    group_by(sample) %>%
    summarize(pct=sum(pct)) %>%
    (function(.) {structure(.$pct, names=.$sample)[names(allBact)]})
            

relAbundBar = function(level, n=5, annot=annotation,
                       omitControl=TRUE,
                       addHuman=FALSE, expansiveOther=FALSE,
                       bactOnly=FALSE,
                       otherName=ifelse(bactOnly, "other known bacteria", "other"),
                       keep=NULL) {
    levMap = c(k=1, p=2, c=3, o=4, f=5, g=6, s=7)
    x = summarizeData(levMap[tolower(substr(level, 1, 1))])
    if (omitControl) {
        x = x[!grepl("blank|mock",
                     annot[gsub("_.*", "", x$sample), "Sample Description"]), ]
    }
    if (bactOnly) {
        bactSpecs = unique(allData[allData$k == "Bacteria", "s"])
        x = x[ , colnames(x) %in% c("sample", bactSpecs)]
    }
    colnames(x) = gsub("_", " ", colnames(x))
    ggd = as.data.frame(pivot_longer(
        x, -sample, names_to=level, values_to="relative abundance"
    ))
    keepers = head(unique(
        ggd[order(ggd$`relative abundance`, decreasing=TRUE), level]
    ), n)
    if (length(keep) > 0) {
        keep = setdiff(keep, keepers)
        keepers = c(keepers[1:(n-length(keep))], keep)
    }
    keepers = sort(keepers)
    ggd[!(ggd[[level]] %in% keepers), level] = otherName
    ggd$taxon = ggd[[level]]
    ggd = ggd %>%
        group_by(sample, taxon) %>%
        summarize(relative_abundance = sum(`relative abundance`)) %>%
        as.data.frame
    if (expansiveOther) {
        addOther = 100 - structure(rowSums(x[ , -1]), names=x[ , 1])
        if (addHuman) {
            addOther = addOther - krakenHS[names(addOther), "Homo sapiens"]
        }
        idxOther = which(ggd$taxon == otherName)
        ggd[idxOther, "relative_abundance"] =
                ggd[idxOther, "relative_abundance"] +
                addOther[ggd[idxOther, "sample"]]
    }    
    colnames(ggd) = gsub("^taxon", level, colnames(ggd))
    ggd$preparation = annot[gsub("_.*", "", ggd$sample), "Sample Description"]
    ggd$kit = gsub("Power(Soil)?", "PowerSoil",
                   gsub(" .*", "", ggd$preparation))
    ggd$site = gsub(".*(extraction|[kK]it) ", "", ggd$preparation)
    if (addHuman) {
        samp = rownames(krakenHS)
        sampDesc = annotation[gsub("_.*$", "", samp), "Sample Description"]
        ggd = rbind(ggd, data.frame(
            sample = samp,
            species = "Homo sapiens",
            relative_abundance = krakenHS[ , "Homo sapiens"],
            preparation = sampDesc,
            kit = gsub("^Power.*", "PowerSoil", gsub(" [kK]it .*$", "", sampDesc)),
            site = gsub(".*blank", "blank", gsub("^.* [kK]it ", "", sampDesc))
        )[samp %in% ggd$sample, ])
        keepers = c("Homo sapiens", keepers)
        ggd[[level]] = factor(
            ggd[[level]],
            levels = c(otherName, setdiff(keepers, "Homo sapiens"), "Homo sapiens")
        )
    } else {
        ggd[[level]] = factor(ggd[[level]], levels=c(keepers, otherName))
    }
    gg = ggplot(ggd, aes_string(x="sample", y="relative_abundance", fill=level))
    gg = gg + facet_grid(. ~ site + kit, scales="free_x", space="free")
    gg = gg + geom_bar(stat="identity", position="stack")
    palf = if (length(keepers) <= 12) {tol} else {
        if (length(keepers) <= 20) {stepped2} else {cols25}
    }
    if (addHuman) {
        gg = gg + scale_fill_manual(
            values = c("lightgray", rev(palf(length(keepers)-1)), "skyblue")
        )
    } else {
        gg = gg + scale_fill_manual(values=c(rev(palf(length(keepers))), "gray"))
    }
    gg = gg + ylab("relative abundance (%)")
    gg = gg + theme(axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    strip.text.x = element_text(size=ifelse(omitControl, 10, 6)))
    return(gg)
}

gg = relAbundBar("species", 19,
                 omitControl=FALSE,
                 addHuman=TRUE, expansiveOther=FALSE, bactOnly=TRUE,
                 keep="Lactobacillus iners")
## png("relative_abundance_barcharts_including_blanks_mocks.png",
##     h=1920*0.5, w=1920*2, res=144*1.25)
print(gg)
## garbage = dev.off()

gg = relAbundBar("species", 19,
                 omitControl=TRUE,
                 addHuman=TRUE, expansiveOther=FALSE, bactOnly=TRUE,
                 keep="Lactobacillus iners")
## png("relative_abundance_barcharts_excluding_blanks_mocks.png",
##     h=1920*0.5, w=1920*2, res=144*1.25)
print(gg)
## garbage = dev.off()


## ## =============================================================================
## write.table(species, "metaphlan_species_wide.tsv",
##             sep="\t", row.names=FALSE, quote=FALSE)
## R.utils::gzip("metaphlan_species_wide.tsv", overwrite=TRUE)

