set_index = function(x, index='index') {
    return(data.frame(x[ , setdiff(colnames(x), index), drop=FALSE],
                      row.names=x[[index]],
                      check.names=FALSE, stringsAsFactors=FALSE))
}

reset_index = function(x, index='index') {
    out = data.frame(index=rownames(x), x,
                     check.names=FALSE, stringsAsFactors=FALSE)
    colnames(out) = gsub('^index$', index, colnames(out))
    return(out)
}

readExcel = function(filename) {
    require(readxl)
    sheets = readxl::excel_sheets(filename)
    x = lapply(sheets, function(.) {
        readxl::read_excel(filename, sheet=., na=c("", "NA", "NaN"))
    })
    x = lapply(x, as.data.frame)
    names(x) = sheets
    return(x)
}

doNMDS = function(x, annot=annotation, seed=123,
                  autotransform=FALSE, trace=FALSE, try=100) {
    require(ggrepel)    
    set.seed(seed)
    mdsOut = metaMDS(x,
                     autotransform=autotransform, trace=trace,
                     try=try, trymax=try)
    ggd = reset_index(mdsOut$points, index="sample")
    ggd$label = as.character(NA)
    ggd$preparation = annot[gsub("_.*", "", ggd$sample), "Sample Description"]
    ggd$kit = gsub("Power(Soil)?", "PowerSoil",
                   gsub(" .*", "", ggd$preparation))
    ggd$site = gsub(".*(extraction|[kK]it) ", "", ggd$preparation)
    ggd = ggd %>% group_by(kit, site) %>% mutate(cent1=mean(MDS1),
                                                 cent2=mean(MDS2))
    ggd = as.data.frame(ggd)
    ggd$dist = sqrt((ggd$MDS1-ggd$cent1)^2 + (ggd$MDS2-ggd$cent2)^2)
    toLabel = (!is.na(ggd$dist)) & (ggd$dist >= 2)
    ggd[toLabel, "label"] = paste0(
        ggd[toLabel, "kit"], ": ", gsub(" swab$", "", ggd[toLabel, "site"])
    )
    gg = ggplot(ggd, aes(x=MDS1, y=MDS2, color=site, shape=kit, label=label))
    gg = gg + geom_point()
    gg = gg + geom_text_repel(max.overlaps=Inf, show.legend=FALSE)
    gg = gg + scale_color_manual(values=heatColors$site %>%
                                        gsub("white", "goldenrod", .) %>%
                                        gsub("lightgray", "darkgray", .) %>%
                                        gsub("steelblue", "dodgerblue", .))
    gg = gg + scale_shape_manual(values=c(24, 25, 22, 23))
    gg = gg + theme(panel.grid.minor=element_blank())
    return(gg)
}

idmap = function(v) {structure(v, names=v)}

`%[[%` = function(df, col) {structure(df[ , col], names=rownames(df))}

nthize = function(v, n) {
    out = rep(-1, length(v))
    offset = 0
    for (i in 1:n) {
        oldInds = seq(i, length(v), n)
        out[offset + (1:length(oldInds))] = v[oldInds]
        offset = offset + length(oldInds)
    }
    return(out)
}
