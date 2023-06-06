allData$sample = gsub("_L00[12]$", "", allData$sample)
colnames(allData) = gsub("pct$", "pct_", colnames(allData))
allData = allData %>%
    group_by(sample, key, ncbi, level, k, p, c, o, f, g, s) %>%
    summarize(pct = sum(pct_) / 2, det_in=length(pct_)) %>%
    as.data.frame

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
