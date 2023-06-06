allData$sample = gsub("_L00[12]$", "", allData$sample)
colnames(allData) = gsub("count$", "count_", colnames(allData))
colnames(allData) = gsub("pct$", "pct_", colnames(allData))
allData = allData %>%
    group_by(sample, key, res,
             taxonomic_level_01, taxonomic_level_02, taxonomic_level_03,
             taxonomic_level_04, taxonomic_level_05, taxonomic_level_06,
             taxonomic_level_07, taxonomic_level_08, taxonomic_level_09,
             taxonomic_level_10, taxonomic_level_11, taxonomic_level_12,
             taxonomic_level_13, taxonomic_level_14, taxonomic_level_15,
             taxonomic_level_16, taxonomic_level_17, taxonomic_level_18,
             taxonomic_level_19, taxonomic_level_20, taxonomic_level_21,
             taxonomic_level_22, taxonomic_level_23, taxonomic_level_24,
             taxonomic_level_25, taxonomic_level_26, taxonomic_level_27,
             taxonomic_level_28, taxonomic_level_29, taxonomic_level_30,
             taxonomic_level_31, taxonomic_level_32
             ) %>%
    summarize(pct = sum(pct_)/2, count=sum(count_)/2, det_in=length(count_)) %>%
    as.data.frame

domainData = allData[allData$res == "D", c("sample", "key", "pct")]
domainData$key = gsub("^ +", "", domainData$key)
domainData = domainData %>%
    pivot_wider(names_from=key, values_from=pct) %>%
    as.data.frame

ppct = summarizeData("P") %>% set_index("sample")
gpct = summarizeData("G") %>% set_index("sample")
spct = summarizeData("S") %>% set_index("sample")

pc = summarizeCounts("P") %>% set_index("sample")
gc = summarizeCounts("G") %>% set_index("sample")
sc = summarizeCounts("S") %>% set_index("sample")

phyla = summarizeData("P")
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
