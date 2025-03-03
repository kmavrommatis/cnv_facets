def_rnd <- "The name of the input file"
def_nbhd <- list("auto", 250)
def_purity_cval <- NULL # 1000
xargs <- list(
    "out" = "/lab/ni/P-20180608-0001/testmethods/facets/EA835175test",
    "snp_tumour" = "/lab/ni/P-20180608-0001/data/BAM/EA835175.sorted.dedup.realigned.recal.hg38.bam",
    "snp_normal" = "/lab/ni/P-20180608-0001/data/BAM/EA1012091A6.sorted.dedup.realigned.recal.hg38.bam",
    "snp_mapq" = 15,
    "snp_baq" = 10,
    "snp_count_orphans" = FALSE,
    "snp_nprocs" = 26,
    "pileup" = "/lab/ni/P-20180608-0001/data/funnel/FACETS_tumor_matched_normal/15afddbd-5da4-4fc9-951f-f26799d9b013/cnv_facets/EA1012091A6.sorted.dedup.realigned.recal.hg38-EA835175.sorted.dedup.realigned.recal.hg38.csv.gz",
    "depth" = c(25, 4000),
    "targets" = "/lab/ni/P-20180608-0001/data//WholeGenome/ExonCapture/S04380110_Covered.bed",
    "cval" = c(25, 0),
    "purity_cval" = 0,
    "nbhd_snp" = "auto",
    "annotation" = NULL,
    "gbuild" = "hg38",
    "unmatched" = FALSE,
    "no_cov_plot" = FALSE,
    "dipLogR" = NULL,
    "rnd_seed" = "The name of the input file2 "
)
