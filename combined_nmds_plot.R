#!/usr/bin/env Rscript

library(ggplot2); theme_set(theme_bw())
library(ggrepel)
library(magrittr)

heatColors = list(
    kit = c(DNeasy="white", HostZero="black",
            PowerSoil="goldenrod", PureLink="gray"),
    site = c(blank="white", mock="lightgray", "oral swab"="lightseagreen",
             "rectal swab"="orangered", "vaginal swab"="steelblue")
)

metaphlan = read.table("metaphlan_nmds_swap.tsv",
                       sep="\t", header=TRUE, row.names=NULL)
metaphlan$analysis = "metaphlan"
kraken = read.table("kraken2_nmds_swap.tsv",
                    sep="\t", header=TRUE, row.names=NULL)
kraken$analysis = "kraken"
ggd = rbind(metaphlan, kraken)

gg = ggplot(ggd, aes(x=MDS1, y=MDS2, color=site, shape=kit, label=label))
gg = gg + facet_wrap(~ analysis, scales="free")
gg = gg + geom_point()
gg = gg + geom_text_repel(max.overlaps=Inf, show.legend=FALSE)
gg = gg + scale_color_manual(values=heatColors$site %>%
                                    gsub("white", "goldenrod", .) %>%
                                    gsub("lightgray", "darkgray", .) %>%
                                    gsub("steelblue", "dodgerblue", .))
gg = gg + scale_shape_manual(values=c(24, 25, 22, 23))
gg = gg + theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank())
## png("bray_curtis_nmds_kraken2_and_metaphlan.png",
##     h=960, w=960*2.225, res=144*1.5)
print(gg)
## garbage = dev.off()
