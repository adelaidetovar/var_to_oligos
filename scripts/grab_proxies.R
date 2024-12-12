options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("librarian")

librarian::shelf(tidyverse, LDlinkR, argparse)

# parse cmd line arguments
parser <- ArgumentParser()

parser$add_argument("-intab", "--signal_list", default = "signal_list.txt", help = "Path to list of signals. Min cols required: 'signal_index_rsid', 'ancestry'")
parser$add_argument("-token", "--ldlink_token", help = "LDLink personal access token - see https://ldlink.nih.gov/?tab=apiaccess")
parser$add_argument("-outdir", "--out_dir", help = "Output directory for LDProxy results")
parser$add_argument("-genomebuild", "--genome_build", default = "grch38_high_coverage", 
	help = "Preferred genome build for LDProxy results. Options: 'grch37', 'grch38', 'grch38_high_coverage' (default)")

args <- parser$parse_args()


# set LDLink API token
Sys.setenv("LDLINK_TOKEN" = args$ldlink_token)

# read in signal list
signals <- read.delim(args$signal_list)
signals <- signals %>%
	select(signal_index_rsid, ancestry) %>%
	mutate(ancestry = tolower(ancestry))


# filter signals with no rsid for manual curation
signals_norsid <- signals %>%
	filter(!(str_detect(signal_index_rsid, regex("rs[0-9A-z]*", ignore_case = T))))
write.table(signals_norsid, paste0(args$out_dir, "/missing_rsid.txt"),
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

dir.create(args$out_dir)
setwd(args$out_dir)

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

