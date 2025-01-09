options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("librarian", quietly = TRUE)) {
	install.packages("librarian")
}

librarian::shelf(argparse, purrr, rsnps)

parser <- ArgumentParser()

parser$add_argument("-mode", "--mode", default = "gwas", help = "Method used to curate variant list and/or info. Options: gwas (default), coord, info")
parser$add_argument("-intab", "--in_tab", nargs = 1,
    help = "If using mode 'gwas': Path to delimited file of signals. Min cols required: 'signal_index_rsid', 'ancestry'. If using mode 'coord': Path to delimited file of coordinates. Min cols required: chr, start, end. If using mode 'info': Path to delimited file of variants. Min cols required: rsid")
parser$add_argument("-ldlinktoken", "--ldlink_token", help = "LDLink personal access token - see https://ldlink.nih.gov/?tab=apiaccess; required for 'gwas' mode")
parser$add_argument("-ncbikey", "--ncbi_key", nargs = 1, help = "NCBI API key - see https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us")
parser$add_argument("-genomebuild", "--genome_build", default = "grch38_high_coverage", 
	help = "Preferred genome build for LDProxy results. Options: 'grch37', 'grch38', 'grch38_high_coverage' (default); required for 'gwas' mode")
parser$add_argument("-r2cutoff", "--r2_cutoff", default = 0.7, help = "R2 cutoff to filter LDLink results (default = 0.7); required for 'gwas' mode")
parser$add_argument("-outdir", "--out_dir", nargs = 1, help = "Full path to output directory for variant lists and info")

args <- parser$parse_args()

# set NCBI API token
Sys.setenv("ENTREZ_KEY" = args$ncbi_key)

# set up NCBI dbSNP function call
snp_function = possibly(ncbi_snp_query, otherwise = "failed")

# create output directory, if needed
if(!file.exists(args$out_dir)) {
    dir.create(args$out_dir)
    }

out_dir = normalizePath(args$out_dir)

if(args$mode == "gwas") {
    
    librarian::shelf(tidyverse, LDlinkR, argparse, httr)

    # set LDLink API token
    Sys.setenv("LDLINK_TOKEN" = args$ldlink_token)

    # read in signal list
    signals <- read.delim(args$in_tab)
    signals <- signals %>%
        select(signal_index_rsid, ancestry) %>%
        mutate(ancestry = tolower(ancestry))

    # filter signals with no rsid for manual curation
    signals_norsid <- signals %>%
        filter(!(str_detect(signal_index_rsid, regex("rs[0-9A-z]*", ignore_case = T))))
    write.table(signals_norsid, file.path(out_dir, "missing_rsid.txt"),
        sep = "\t", row.names = F, quote = F)

    signals <- signals %>%
        filter(str_detect(signal_index_rsid, regex("rs[0-9A-z]*", ignore_case = T)))


    # prep for AFR
    if(any(signals$ancestry %in% c("african", "afr", "multi", "te", "ta"))) {
        ldlink_AFR <- signals %>%
            filter(signals$ancestry %in% c("african", "afr", "multi", "te", "ta")) %>%
            unique()
    }

    # prep for AMR
    if(any(signals$ancestry %in% c("american", "amerindigenous", "amr", "multi", "te", "ta"))) {
        ldlink_AMR <- signals %>%
            filter(signals$ancestry %in% c("american", "amerindigenous", "amr", "multi", "te", "ta")) %>%
            unique()
    }

    # prep for EAS
    if(any(signals$ancestry %in% c("eastasian", "east_asian", "eas", "multi", "te", "ta"))) {
        ldlink_EAS <- signals %>%
            filter(signals$ancestry %in% c("eastasian", "eas", "multi", "te", "ta")) %>%
            unique()
    }

    # prep for EUR
    if(any(signals$ancestry %in% c("european", "eur", "multi", "te", "ta"))) {
        ldlink_EUR <- signals %>%
            filter(signals$ancestry %in% c("european", "eur", "multi", "te", "ta")) %>%
            unique()
    }

    # prep for SAS
    if(any(signals$ancestry %in% c("southasian", "south_asian", "sas", "multi", "te", "ta"))) {
        ldlink_SAS <- signals %>%
            filter(signals$ancestry %in% c("southasian", "south_asian", "sas", "multi", "te", "ta")) %>%
            unique()
    }

    # LDLink API can only be accessed sequentially!
    dir.create(file.path(out_dir, "ldlink_out"))
    setwd(file.path(out_dir, "ldlink_out"))

    if(exists("ldlink_AFR")) {
        LDproxy_batch(ldlink_AFR$signal_index_rsid,
        pop = "AFR",
        r2d = "r2",
        token = Sys.getenv("LDLINK_TOKEN"),
        append = TRUE,
        genome_build = args$genome_build)
        file.rename(paste0("combined_query_snp_list_", args$genome_build, ".txt"),
            paste0("combined_query_snp_list_AFR.txt"))
    }

    if(exists("ldlink_AMR")) {
        LDproxy_batch(ldlink_AMR$signal_index_rsid,
        pop = "AMR",
        r2d = "r2",
        token = Sys.getenv("LDLINK_TOKEN"),
        append = TRUE,
        genome_build = args$genome_build)
        file.rename(paste0("combined_query_snp_list_", args$genome_build, ".txt"),
            paste0("combined_query_snp_list_AMR.txt"))
    }

    if(exists("ldlink_EAS")) {
        LDproxy_batch(ldlink_EAS$signal_index_rsid,
        pop = "EAS",
        r2d = "r2",
        token = Sys.getenv("LDLINK_TOKEN"),
        append = TRUE,
        genome_build = args$genome_build)
        file.rename(paste0("combined_query_snp_list_", args$genome_build, ".txt"),
            paste0("combined_query_snp_list_EAS.txt"))
    }

    if(exists("ldlink_EUR")) {
        LDproxy_batch(ldlink_EUR$signal_index_rsid,
        pop = "EUR",
        r2d = "r2",
        token = Sys.getenv("LDLINK_TOKEN"),
        append = TRUE,
        genome_build = args$genome_build)
        file.rename(paste0("combined_query_snp_list_", args$genome_build, ".txt"),
            paste0("combined_query_snp_list_EUR.txt"))
    }

    if(exists("ldlink_SAS")) {
        LDproxy_batch(ldlink_SAS$signal_index_rsid,
        pop = "SAS",
        r2d = "r2",
        token = Sys.getenv("LDLINK_TOKEN"),
        append = TRUE,
        genome_build = args$genome_build)
        file.rename(paste0("combined_query_snp_list_", args$genome_build, ".txt"),
            paste0("combined_query_snp_list_SAS.txt"))
    }

    files <- list.files(file.path(out_dir, "ldlink_out"), full.names = TRUE)

	proxies <- list()
	for (i in 1:length(files)){
		proxies[[i]] <- files[[i]] %>%
		setNames(nm = .) %>%
		map_df(~read_delim(.x, col_types = cols(), col_names = TRUE, delim = "\t"), .id = "file_name")
		colnames(proxies[[i]]) <- c(colnames(proxies[[i]])[1], "row", colnames(proxies[[i]])[2:13])
	}

	for (i in 1:length(files)){
		proxies[[i]] <- proxies[[i]] %>%
		filter(., R2 >= args$r2_cutoff) %>%
		mutate(ancestry = str_extract(file_name, pattern = "[A-Z]*(?=.txt)"))
	}

	# bind into union table, filter for unique variants while keeping distinct ref panel info
	proxies_union <- bind_rows(proxies, .id = "column_label") %>%
		group_by(RS_Number) %>%
		reframe(signal = query_snp,
			coord = Coord,
			alleles = Alleles,
			MAF = toString(sort(unique(MAF)), .groups='drop'),
			R2 = toString(sort(unique(R2)), .groups='drop'),
			ancestry = toString(sort(unique(ancestry)), .groups='drop')) %>%
		distinct(RS_Number, .keep_all = TRUE) %>%
		rename(rsid = RS_Number)

	write.table(proxies_union, file.path(out_dir, "all_var.tsv"), sep = "\t", quote = F, row.names = F)

	# query results, return 'failed' if query throws an error
	results = map(proxies_union$rsid, snp_function)
	names(results) = proxies_union$rsid

	for (i in 1:length(names(results))){
		ifelse(dim(results[[i]][1]) == 0, results[[i]][1,] <- c(names(results)[i],rep(NA, 14)), results[[i]][1,])
	}
} else if(args$mode == "coord") {
    librarian::shelf(tidyverse, LDlinkR, SNPlocs.Hsapiens.dbSNP155.GRCh38, GenomicRanges, httr)

    # Check chr and start/end formatting
    df <- read.delim(args$in_tab)
    df <- df %>%
    dplyr::mutate(chr = case_when(
        chr == "23" | chr == "chr23" ~ "X",
        chr == "24" | chr == "chr24" ~ "Y",
        grepl("^chr", chr) ~ sub("^chr", "", chr),
        TRUE ~ chr
    ),
    start = as.numeric(start),
    end = as.numeric(end))

    gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start,
                                end = df$end))

    var_gr <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)

    var_df <- data.frame("chr" = var_gr@seqnames,
                        "bp_hg38" = var_gr@ranges,
                        "rsid" = var_gr@elementMetadata$RefSNP_id)

    write.table(var_df, file.path(out_dir, "var_by_coord.txt"),
        sep = "\t", row.names = F, quote = F)

	# query results, return 'failed' if query throws an error
	results = map(var_df$rsid, snp_function)
	names(results) = var_df$rsid

    for (i in 1:length(names(results))){
		ifelse(dim(results[[i]][1]) == 0, results[[i]][1,] <- c(names(results)[i],rep(NA, 14)), results[[i]][1,])
	}
} else if(args$mode == "info") {
    librarian::shelf(tidyverse)
    var_df <- read.delim(args$in_tab)

    # query results, return 'failed' if query throws an error
	results = map(var_df$rsid, snp_function)
	names(results) = var_df$rsid
} else {
    return("Please re-run and specify variant gathering mode! (gwas, coord)")
}

write.table(bind_rows(results)[,c(1:14)], file.path(out_dir, "var_info.tsv"), sep="\t", quote = F, row.names = F)