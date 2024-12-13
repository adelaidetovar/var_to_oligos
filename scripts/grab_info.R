options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("rsnps", quietly = TRUE)) {
	remotes::install_github("ropensci/rsnps")
}
librarian::shelf(tidyverse, argparse, rsnps, purrr, httr)

parser <- ArgumentParser()

parser$add_argument("-indir", "--in_dir", help = "Path to list of LDLink results")
parser$add_argument("-token", "--ncbi_key", help = "NCBI API key - see https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us")
parser$add_argument("-outdir", "--out_dir", help = "Output directory for variant info")
parser$add_argument("-r2cutoff", "--r2_cutoff", default = 0.7, help = "R2 cutoff to filter LDLink results")

args <- parser$parse_args()

# set NCBI API token
Sys.setenv("ENTREZ_KEY" = args$ncbi_key)

# read in all variants from LDLink query
files <- list.files(args$in_dir, pattern = "combined*", full.names = TRUE)

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

write.table(proxies_union, paste0(args$in_dir, "/proxies_union.tsv"), sep = "\t", quote = F, row.names = F)

# query results, return 'failed' if query throws an error
snp_function = possibly(ncbi_snp_query, otherwise = "failed")
results = map(proxies_union$rsid, snp_function)
names(results) = proxies_union$rsid

for (i in 1:length(names(results))){
  ifelse(dim(results[[i]][1]) == 0, results[[i]][1,] <- c(names(results)[i],rep(NA, 14)), results[[i]][1,])
}

dir.create(args$out_dir)
write.table(bind_rows(results)[,c(1:14)], paste0(args$out_dir, "/proxies_info.tsv"), sep="\t", quote = F, row.names = F)

