# *Supplementary code for:* Comparison of commercial DNA extraction kits for whole metagenome sequencing of human oral, vaginal, and rectal microbiome samples

This repository includes bash and R scripts used for the analysis of the data described in the paper ``Comparison of commercial DNA extraction kits for whole metagenome sequencing of human oral, vaginal, and rectal microbiome samples``.

## Initial processing with `kraken2` and `metaphlan`

The scripts [run\_kraken2.sh](run_kraken2.sh) and [run\_metaphlan.sh](run_metaphlan.sh) contain the commands used to run [**kraken2**](https://ccb.jhu.edu/software/kraken2/) and [**MetaPhlAn**](https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-3.0), respectively.

## R scripts used for downstream analysis of processed results

The scripts

- [kraken2\_analysis\_merged\_lanes\_swap.R](kraken2_analysis_merged_lanes_swap.R),
- [metaphlan\_analysis\_merged\_lanes\_swap.R](metaphlan_analysis_merged_lanes_swap.R), and
- [metacal-ing.R](metacal-ing.R)

which in turn rely on the code in the supporting R source files [utilities.R](utilities.R), [parse\_kraken2\_report\_only.R](parse_kraken2_report_only.R), [kraken2\_merge\_lanes.R](kraken2_merge_lanes.R), [parse\_metaphlan\_report\_only.R](parse_metaphlan_report_only.R), and [metaphlan\_merge\_lanes.R](metaphlan_merge_lanes.R), contain the bulk of the analyses discussed.

The remaining scripts

- [barcharting.R](barcharting.R) and
- [combined\_nmds\_plot.R](combined_nmds_plot.R)

contain the code used to generate the barchart and NMDS plot figures.
