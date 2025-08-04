avg_insert_size <- function(bam, default) {
    idx <- data.table(idxstatsBam(bam))[mapped > 0][order(-mapped)]
    if (nrow(idx) > 10) {
        idx <- idx[1:10, ]
    }
    ins_size <- c()
    n_reads <- c()
    hdr_len <- length(scanBamHeader(bam)[[1]][["text"]])
    for (chrom in idx$seqnames) {
        cmd <- sprintf(
            "samtools view -q 3 -f 3 -F 3840 -h %s %s | head -n %s | samtools stats | grep ^SN | cut -f 2-",
            bam, chrom, 200000 + hdr_len
        )
        stats <- system(cmd, intern = TRUE)
        size <- grep("insert size average", stats, value = TRUE)
        size <- as.numeric(unlist(strsplit(size, "\t"))[2])
        n <- grep("reads properly paired", stats, value = TRUE)
        n <- as.numeric(unlist(strsplit(n, "\t"))[2])
        if (n == 0) {
            next
        }
        ins_size <- c(ins_size, size)
        n_reads <- c(n_reads, n)
    }
    if (length(n_reads) == 0) {
        return(as.numeric(default))
    }
    return(weighted.mean(ins_size, n_reads))
}

exec_snp_pileup <- function(chrom, snp_vcf, output, normal_bam, tumour_bam, mapq, baq, pseudo_snp, keep_orphans) {
    # Execute snp-pileup on chromosome `chrom`
    # Send tmp output to the same directory of the final output so we are sure
    # we can write there.
    d <- dirname(output)
    chrom_vcf <- file.path(d, paste0(sub("\\.vcf\\.gz$|\\.vcf\\.bgz$", "", basename(snp_vcf)), ".", chrom, ".vcf"))
    chrom_nbam <- file.path(d, paste0(sub("\\.bam", "", basename(normal_bam)), ".", chrom, ".bam"))
    chrom_tbam <- file.path(d, paste0(sub("\\.bam", "", basename(tumour_bam)), ".", chrom, ".bam"))

    if (keep_orphans) {
        orphans <- "--count-orphans"
    } else {
        orphans <- ""
    }

    cmd <- c(
        "#!/bin/bash", "\n",
        "\n",
        "set -eo pipefail", "\n",
        "\n",
        "mkfifo", chrom_vcf, "\n",
        "mkfifo", chrom_nbam, "\n",
        "mkfifo", chrom_tbam, "\n",
        "\n",
        "bcftools view --output-type u", snp_vcf, chrom, ">", chrom_vcf, "&\n",
        "pid_bcf=$!", "\n",
        "\n",
        "samtools view -u", normal_bam, chrom, ">", chrom_nbam, "&\n",
        "pid_nbam=$!", "\n",
        "\n",
        "samtools view -u", tumour_bam, chrom, ">", chrom_tbam, "&\n",
        "pid_tbam=$!", "\n",
        "\n",
        "snp-pileup",
        "--gzip",
        "--pseudo-snps", pseudo_snp,
        "--min-map-quality", mapq,
        "--min-base-quality", baq,
        "--max-depth 10000000",
        "--min-read-counts", "0,0",
        orphans,
        chrom_vcf, output, chrom_nbam, chrom_tbam,
        "\n",
        "set +e", "\n",
        "\n",
        "wait $pid_bcf; exit_code=$?", "\n",
        "if [[ $exit_code != 0 && $exit_code != 141 ]]; then exit $exit_code; fi", "\n",
        "\n",
        "wait $pid_nbam; exit_code=$?", "\n",
        "if [[ $exit_code != 0 && $exit_code != 141 ]]; then exit $exit_code; fi", "\n",
        "\n",
        "wait $pid_tbam; exit_code=$?", "\n",
        "if [[ $exit_code != 0 && $exit_code != 141 ]]; then exit $exit_code; fi", "\n"
    )
    cmd <- paste(cmd, collapse = " ")
    script <- file.path(d, sprintf("snp_pileup.%s.sh", chrom))
    write(cmd, script)
    status <- system2("/bin/bash", args = script)
    if (status != 0) {
        stop(sprintf("\nError in computing snp pileup. Exit code %s from execution of %s\n", status, script))
    }
    unlink(c(chrom_vcf, chrom_nbam, chrom_tbam))
    return(cmd)
}

exec_snp_pileup_parallel <- function(snp_vcf, output, normal_bam, tumour_bam, mapq, baq, pseudo_snp, nprocs, keep_orphans) {
    dtm <- format(Sys.time(), "%y%m%d-%H%M%S")
    tmpdir <- tempfile(pattern = paste0(basename(sub("\\.cvs\\.gz$", "", output)), "_", dtm, "_"), tmpdir = dirname(output))
    dir.create(tmpdir)

    chroms <- headerTabix(snp_vcf)$seqnames

    cl <- makeCluster(nprocs, type = "FORK")
    chrom_csv <- parLapply(cl, chroms, function(chrom) {
        chrom_csv <- file.path(tmpdir, paste0(chrom, ".csv.gz"))
        cmd <- exec_snp_pileup(
            chrom = chrom,
            snp_vcf = snp_vcf,
            output = chrom_csv,
            normal_bam = normal_bam,
            tumour_bam = tumour_bam,
            mapq = mapq,
            baq = baq,
            pseudo_snp = pseudo_snp,
            keep_orphans = keep_orphans
        )
        return(chrom_csv)
    })
    stopCluster(cl)
    concat_csv(chrom_csv, output, tmpdir = tmpdir)
    unlink(tmpdir, recursive = TRUE)
}

concat_csv <- function(csv_list, xfile, tmpdir) {
    isGzip <- ifelse(grepl("\\.gz$", xfile), TRUE, FALSE)
    if (isGzip == TRUE) {
        conn <- file.path(tmpdir, sub("\\.gz", "", basename(xfile)))
    } else {
        conn <- xfile
    }
    col.names <- TRUE
    for (csv in csv_list) {
        options(datatable.fread.input.cmd.message = FALSE)
        fwrite(x = fread(sprintf("gzip -d -c %s", csv)), file = conn, sep = ",", col.names = col.names, row.name = FALSE, quote = FALSE, append = !isTRUE(col.names))
        options(datatable.fread.input.cmd.message = TRUE)
        col.names <- FALSE
    }
    if (isGzip == TRUE) {
        system2(c("gzip", "-c", conn), stdout = xfile)
        unlink(conn)
    }
}

readSnpMatrix2 <- function(pileup, gbuild) {
    xf <- file(pileup, open = "r")
    if (summary(xf)$class == "gzfile") {
        conn <- sprintf("gzip -d -c %s", pileup)
    } else {
        conn <- pileup
    }
    close(xf)

    options(datatable.fread.input.cmd.message = FALSE)
    rcmat <- fread(conn, select = c("Chromosome", "Position", "File1R", "File1A", "File2R", "File2A"))
    options(datatable.fread.input.cmd.message = TRUE)

    setnames(
        rcmat, c("File1R", "File1A", "File2R", "File2A"),
        c("NOR.RD", "NOR.DP", "TUM.RD", "TUM.DP")
    )
    rcmat[, NOR.DP := NOR.DP + NOR.RD]
    rcmat[, TUM.DP := TUM.DP + TUM.RD]

    chr_prefix <- any(rcmat$Chromosome %in% c(paste0("chr", 1:22), "chrX"))

    # Unfortunately, facets needs numeric chromsomes. X will be converted later by facets
    rcmat[, Chromosome := sub("^chr", "", Chromosome)]
    setcolorder(rcmat, c("Chromosome", "Position", "NOR.DP", "NOR.RD", "TUM.DP", "TUM.RD"))

    # We only keep the major chromosomes
    if (gbuild %in% c("hg19", "hg38")) {
        rcmat <- rcmat[Chromosome %in% c(1:22, "X")]
    } else if (gbuild %in% c("mm9", "mm10")) {
        rcmat <- rcmat[Chromosome %in% c(1:19, "X")]
    } else {
        write(sprintf("Invalid genome build: %s", gbuild), stderr())
        quit(status = 1)
    }
    return(list(pileup = rcmat, chr_prefix = chr_prefix))
}

facetsRecordToVcf <- function(x) {
    stopifnot(is.data.table(x))
    stopifnot(nrow(x) == 1)
    # Convert the annotated facets record to a VCF record.
    vcf <- vector(length = 8)
    vcf[1] <- x$chrom
    vcf[2] <- format(x$start + 1, scientific = FALSE) # TODO: Check +1 is correct
    vcf[3] <- x$seg
    vcf[4] <- "N"
    vcf[5] <- "<CNV>"
    vcf[6] <- "."
    vcf[7] <- ifelse(x$type == "NEUTR", "neutral", "PASS")

    if (abs(x$start - x$end) < 2) {
        return(NULL)
    }
    if (x$start > x$ end) {
        tt <- x$start
        x$start <- x$end
        x$end <- tt
    }

    # INFO field: Keep consistent with header
    vcf[8] <- paste0(
        "SVTYPE=", x$type,
        ";SVLEN=", x$end - x$start,
        ";END=", x$end,
        ";NUM_MARK=", x$num.mark,
        ";NHET=", x$nhet,
        ";CNLR_MEDIAN=", ifelse(is.na(x$cnlr.median), ".", round(x$cnlr.median, 3)),
        ";MAF_R=", ifelse(is.na(x$mafR), ".", round(x$mafR, 3)),
        ";SEGCLUST=", x$segclust,
        ";CNLR_MEDIAN_CLUST=", ifelse(is.na(x$cnlr.median.clust), ".", round(x$cnlr.median.clust, 3)),
        ";MAF_R_CLUST=", ifelse(is.na(x$mafR.clust), ".", round(x$mafR.clust, 3)),
        ";CF_EM=", ifelse(is.na(x$cf.em), ".", round(x$cf.em, 3)),
        ";TCN_EM=", ifelse(is.na(x$tcn.em), ".", x$tcn.em),
        ";LCN_EM=", ifelse(is.na(x$lcn.em), ".", x$lcn.em),
        ";CNV_ANN=", ifelse(is.na(x$annotation) || is.null(x$annotation) || x$annotation == "", ".", x$annotation)
    )
    return(vcf)
}

getScriptName <- function() {
    opt <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    name <- basename(sub("^--file=", "", opt))
    return(name)
}

annotate <- function(cnv, bed_file) {
    # Annotate data.table cnv with the features in bed_file

    ann <- data.table(read.table(bed_file, comment.char = "#", header = FALSE, sep = "\t", stringsAsFactors = FALSE, na.strings = ""))
    stopifnot(ncol(ann) >= 4)
    ann <- ann[, 1:4]
    setnames(ann, names(ann), c("chrom", "start", "end", "name"))
    ann <- ann[!is.na(name)]

    # Make feature names VCF compliant. From https://samtools.github.io/hts-specs/VCFv4.1.pdf:
    #
    #     String, no white-space, semi-colons, or equals-signs permitted; commas are
    #     permitted only as delimiters for lists of values
    #
    # Convert special characters using URL encoding
    ann[, name := gsub("%", "%25", name, fixed = TRUE)]
    ann[, name := gsub(",", "%2C", name, fixed = TRUE)]
    ann[, name := gsub("=", "%3D", name, fixed = TRUE)]
    ann[, name := gsub(";", "%3B", name, fixed = TRUE)]
    ann[, name := gsub("|", "%7C", name, fixed = TRUE)]
    ann[, name := gsub(" ", "_", name, fixed = TRUE)] # NB: Not URL encoding for spaces!

    ann <- makeGRangesFromDataFrame(ann,
        seqnames.field = "chrom",
        start.field = "start", end.field = "end", keep.extra.columns = TRUE,
        starts.in.df.are.0based = TRUE
    )

    gcnv <- makeGRangesFromDataFrame(cnv, seqnames.field = "chrom", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
    suppressWarnings({
        ovl <- findOverlaps(query = gcnv, subject = ann, ignore.strand = TRUE)
    })
    hits <- data.table(queryHits = ovl@from, subjectHits = ovl@to, feature = ann$name[ovl@to])
    hits <- hits[, list(feature = paste(feature, collapse = ",")), by = queryHits]
    gcnv <- as.data.table(gcnv)
    gcnv[, annotation := NA]
    gcnv$annotation[hits$queryHits] <- hits$feature
    gcnv$annotation[gcnv$type == "NEUTR"] <- NA
    setnames(gcnv, "seqnames", "chrom")
    gcnv[, chrom := as.character(chrom)]

    stopifnot(cnv$chrom == gcnv$chrom)
    stopifnot(cnv$start == gcnv$start)
    stopifnot(cnv$end == gcnv$end)
    cnv[, annotation := gcnv$annotation]
    return(cnv)
}

make_header <- function(gbuild, genomes, is_chrom_prefixed, cmd, extra) {
    # extra: Named vector of additional information. E.g. c(purity=0.5, ploidy= 2.1)
    header <- c(
        "##fileformat=VCFv4.2",
        sprintf("##reference=%s", gbuild)
    )
    # Contigs
    # --------------------------
    chrom_size <- NULL
    if (gbuild == "hg19") {
        chrom_size <- genomes$HG19
    }
    if (gbuild == "hg38") {
        chrom_size <- genomes$HG38
    }
    if (gbuild == "mm9") {
        chrom_size <- genomes$MM9
    }
    if (gbuild == "mm10") {
        chrom_size <- genomes$MM10
    }

    stopifnot(!is.null(chrom_size))

    for (i in 1:length(chrom_size)) {
        size <- chrom_size[i]
        name <- names(chrom_size)[i]
        if (is_chrom_prefixed != TRUE) {
            name <- sub("^chr", "", name)
        }
        header <- c(header, sprintf("##contig=<ID=%s,length=%s>", name, size))
    }
    # --------------------------

    header <- c(header, '##FILTER=<ID=PASS,Description="All filters passed">')
    header <- c(header, '##FILTER=<ID=neutral,Description="Copy number neutral">')
    header <- c(header, '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header <- c(header, '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
    header <- c(header, '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">')
    header <- c(header, '##INFO=<ID=NUM_MARK,Number=1,Type=Integer,Description="Number of SNPs in the segment">')
    header <- c(header, '##INFO=<ID=NHET,Number=1,Type=Integer,Description="Number of SNPs that are deemed heterozygous">')
    header <- c(header, '##INFO=<ID=CNLR_MEDIAN,Number=1,Type=Float,Description="Median log-ratio (logR) of the segment. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal">')
    header <- c(header, '##INFO=<ID=CNLR_MEDIAN_CLUST,Number=1,Type=Float,Description="Median log-ratio (logR) of the segment cluster. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal">')
    header <- c(header, '##INFO=<ID=MAF_R,Number=1,Type=Float,Description="Log-odds-ratio (logOR) summary for the segment. logOR is defined by the log-odds ratio of the variant allele count in the tumor versus in the normal">')
    header <- c(header, '##INFO=<ID=MAF_R_CLUST,Number=1,Type=Float,Description="Log-odds-ratio (logOR) summary for the segment cluster. logOR is defined by the log-odds ratio of the variant allele count in the tumor versus that in the normal">')
    header <- c(header, '##INFO=<ID=SEGCLUST,Number=1,Type=Integer,Description="Segment cluster to which the segment belongs">')
    header <- c(header, '##INFO=<ID=CF_EM,Number=1,Type=Float,Description="Cellular fraction, fraction of DNA associated with the aberrant genotype. Set to 1 for normal diploid">')
    header <- c(header, '##INFO=<ID=TCN_EM,Number=1,Type=Integer,Description="Total copy number. 2 for normal diploid">')
    header <- c(header, '##INFO=<ID=LCN_EM,Number=1,Type=Integer,Description="Lesser (minor) copy number. 1 for normal diploid">')
    header <- c(header, '##INFO=<ID=CNV_ANN,Number=.,Type=String,Description="Annotation features assigned to this CNV">')
    header <- c(header, '##ALT=<ID=CNV,Description="Copy number variable region">')


    header <- c(header, cmd)

    stopifnot(!is.null(names(extra)))
    for (i in 1:length(extra)) {
        key <- names(extra)[i]
        value <- extra[i]
        stopifnot(key != "")
        header <- c(header, sprintf("##%s=%s", key, value))
    }
    header <- c(header, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return(paste(header, collapse = "\n"))
}

classify_cnv <- function(dat) {
    # Classify CNV. See also https://github.com/mskcc/facets/issues/62
    # dat is a data.table modifed in-place
    dat[, type := NA]
    dat[, type := ifelse((tcn.em == 2 & (lcn.em == 1 | is.na(lcn.em))), "NEUTR", type)]
    dat[, type := ifelse(is.na(type) & tcn.em == 2 & lcn.em == 2, "DUP", type)]
    dat[, type := ifelse(is.na(type) & tcn.em == 0, "DEL", type)]
    dat[, type := ifelse(is.na(type) & tcn.em > 2 & (lcn.em > 0 | is.na(lcn.em)), "DUP", type)]
    dat[, type := ifelse(is.na(type) & tcn.em == 1, "HEMIZYG", type)]
    dat[, type := ifelse(is.na(type) & tcn.em == 2 & lcn.em == 0, "LOH", type)]
    dat[, type := ifelse(is.na(type) & tcn.em > 2 & lcn.em == 0, "DUP-LOH", type)]
    stopifnot(all(!is.na(dat$type))) # Everything has been classified

    return(dat)
}

prep_coverage_data <- function(rcmat) {
    # Reformat the rcmat data.table to make it suitable for plotting histogram
    dhist <- melt(
        data = rcmat, id = "Position", measure.vars = c("NOR.DP", "TUM.DP"),
        variable.name = "sample", value.name = "depth"
    )
    dhist[, Position := NULL]
    dhist[, sample := ifelse(sample == "NOR.DP", "Normal", "Tumour")]
    dhist <- dhist[depth > 0]
    dhist <- dhist[, list(depth = ifelse(.SD$depth > quantile(.SD$depth, 0.99), quantile(.SD$depth, 0.99), .SD$depth)), by = sample]
    nsites <- dhist[, .N, by = sample]
    dhist <- merge(dhist, nsites)
    dhist[, label := paste(sample, "N=", N)]
    dhist[, sample := NULL]
    return(dhist)
}

plot_coverage <- function(rcmat, rcmat_flt, fname, title) {
    if(nrow(rcmat)> 1e7) {
        idx=sample( seq(1, nrow(rcmat)), 1e7, replace = FALSE)
        rcmat= rcmat[idx, ]
    }
    if(nrow(rcmat_flt)> 1e7) {
        idx=sample( seq(1, nrow(rcmat_flt)), 1e7, replace = FALSE)
        rcmat_flt= rcmat_flt[idx, ]
    }

    xall <- ggplot(data = prep_coverage_data(rcmat), aes(x = depth)) +
        geom_histogram(bins = 20, colour = "white", fill = "#F8766D") +
        facet_wrap(~label) +
        xlab("All positions [capped]") +
        theme(axis.title.x = element_text(colour = "#F8766D"))

    xflt <- ggplot(data = prep_coverage_data(rcmat_flt), aes(x = depth)) +
        geom_histogram(bins = 20, colour = "white", fill = "#00BFC4") +
        facet_wrap(~label) +
        xlab("Filtered positions [capped]") +
        theme(axis.title.x = element_text(colour = "#00BFC4"))

    pdf(NULL) # Prevent Rplots.pdf to be generated
    gg <- arrangeGrob(xall, xflt, top = title)
    ggsave(fname, gg, width = 16, height = 18, units = "cm")
    rm(xall, xflt, gg)
    x_ <- gc(verbose = FALSE)
}

filter_rcmat <- function(rcmat, min_ndepth, max_ndepth, target_bed) {
    rcmat_flt <- rcmat[NOR.DP >= min_ndepth & NOR.DP < max_ndepth]
    if (is.null(target_bed) || nrow(rcmat_flt) == 0) {
        return(rcmat_flt)
    }
    targets <- makeGRangesFromDataFrame(target_bed,
        seqnames.field = "V1",
        start.field = "V2", end.field = "V3", starts.in.df.are.0based = TRUE
    )

    # Convert rcmat to GRanges
    rcmat_flt <- makeGRangesFromDataFrame(rcmat_flt,
        seqnames.field = "Chromosome",
        start.field = "Position", end.field = "Position", keep.extra.columns = TRUE
    )

    hits <- findOverlaps(query = rcmat_flt, subject = targets, ignore.strand = TRUE)

    rcmat_flt <- as.data.table(rcmat_flt[unique(hits@from)])
    rcmat_flt[, start := NULL]
    rcmat_flt[, width := NULL]
    rcmat_flt[, strand := NULL]
    setnames(rcmat_flt, c("seqnames", "end"), c("Chromosome", "Position"))
    # This is to print (i.e., for debugging) the data.table obj without calling
    # print twice. See 2.23 at
    # https://cran.r-project.org/web/packages/data.table/vignettes/datatable-faq.html#why-do-i-have-to-type-dt-sometimes-twice-after-using-to-print-the-result-to-console
    rcmat_flt[]
    return(rcmat_flt)
}

run_facets <- function(pre_rcmat,
                       pre_gbuild,
                       pre_snp.nbhd,
                       pre_het.thresh,
                       pre_cval,
                       pre_deltaCN,
                       pre_unmatched,
                       pre_ndepth,
                       pre_ndepthmax,
                       proc_cval,
                       proc_min.nhet,
                       proc_dipLogR,
                       emcncf_unif,
                       emcncf_min.nhet,
                       emcncf_maxiter,
                       emcncf_eps) {
    # Run the core functions of facets for segmentation, purity etc.
    # Here is where the actual CNV discovery happen.
    # Param prefix matches the facets function they go to.
    write(sprintf("[%s] Preprocessing sample...", Sys.time()), stderr())
    set.seed(1234)
    xx <- preProcSample(
        rcmat = pre_rcmat,
        gbuild = pre_gbuild,
        snp.nbhd = pre_snp.nbhd,
        het.thresh = pre_het.thresh,
        cval = pre_cval,
        deltaCN = pre_deltaCN,
        unmatched = pre_unmatched,
        ndepth = pre_ndepth,
        ndepthmax = pre_ndepthmax
    )
    rm(pre_rcmat)
    x_ <- gc(verbose = FALSE)



    write(sprintf("[%s] Processing sample...", Sys.time()), stderr())
    proc_out <- procSample(xx,
        cval = proc_cval,
        min.nhet = proc_min.nhet,
        dipLogR = proc_dipLogR
    )
    proc_out$jointseg <- data.table(proc_out$jointseg)
    proc_out$out <- data.table(proc_out$out)

    write(sprintf("[%s] Fitting model...", Sys.time()), stderr())
    emcncf_fit <- emcncf(
        x = proc_out,
        unif = emcncf_unif,
        min.nhet = emcncf_min.nhet,
        maxiter = emcncf_maxiter,
        eps = emcncf_eps
    )
    names(emcncf_fit$purity) <- NULL
    names(emcncf_fit$ploidy) <- NULL
    emcncf_fit[["cncf"]] <- data.table(emcncf_fit[["cncf"]])[order(chrom, start)]
    return(list(proc_out = proc_out, emcncf_fit = emcncf_fit))
}

reset_chroms <- function(cncf, gbuild, chr_prefix) {
    # Reset chomosome names *in-place*
    # Reset chrom X
    if (gbuild %in% c("hg19", "hg38")) {
        stopifnot(length(unique(cncf$chrom)) <= 23)
        cncf[, chrom := ifelse(chrom == 23, "X", chrom)]
    } else if (gbuild %in% c("mm9", "mm10")) {
        stopifnot(length(unique(cncf$chrom)) <= 20)
        cncf[, chrom := ifelse(chrom == 20, "X", chrom)]
    } else {
        stop(sprintf("Invalid genome build: %s", gbuild))
    }
    # Reset chrom names
    stopifnot(is.logical(chr_prefix))
    if (chr_prefix == TRUE) {
        cncf[, chrom := paste0("chr", chrom)]
    }
}
