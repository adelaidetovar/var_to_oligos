# Variants to MPRA oligos

This repo contains R scripts to make MPRA oligos from GWAS variants. It requires one of the three following items: a list of GWAS signals and corresponding ancestries to perform LD expansion (mode: 'gwas'); a set of hg38 coordinates to scan for common variants (mode: 'coord'); OR a complete list of variants for which you want to gather info about position, alleles, etc (mode: 'info').

Example inputs, job submission scripts, and outputs for each of the three usage modes are available in the `examples` directory.

Prior to usage, you will need to obtain an NCBI API key following instructions [here](https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us) for use with the `rsnps` package. Additionally, if you are using 'gwas' mode, you will need to obtain an LDLink API token following instructions [here](https://ldlink.nih.gov/?tab=apiaccess).


## To start and note on use

Clone repo to wherever you like to store software.

`git clone https://github.com/adelaidetovar/var_to_oligos.git`

If you don't already have a working tidyverse install in R, check to make sure that you've removed any reference to anaconda, miniconda, etc. in your `PATH` variable prior to running this code. Otherwise, it'll think you don't have the right C/C++ libraries installed.

## Grab variants (optional) and variant info

The `grab_proxies.R` script is meant to be run from the command line wrapped with `sbatch` or in a job submission script. Your input signal list should have two columns: `signal_index_rsid` and `ancestry`, where `ancestry` corresponds to one of the 1000 Genomes superpopulation IDs (AFR, AMR, EAS, EUR, SAS) or from a trans-ethnic/-ancestry analysis (TE/TA). This script will need to be modified if you want to use one of the subpopulations.

For example:

```
sbatch --account=[acct] --ntasks=1 --mem-per-cpu=4GB --time=00:10:00 --wrap "Rscript /path/to/grab_proxies.R -intab /path/to/signal_list.txt -token [ldlink_token] -genomebuild [select one: grch37, grch38, grch38_high_coverage (default)] -outdir /path/to/snp_outdir"
```

Adjust paths and parameters (in square brackets) as needed. 

## Design oligos

Unlike the previous script, the oligo design process is outlined in an R Notebook (`make_oligos.Rmd`). You should run this process interactively because there are several steps that require proofreading and (potentially) manual correction. This notebook includes an array of scripts to correct any restriction enzyme sites based on a protocol developed by Adelaide Tovar and Jacob Kitzman -- if you use different enzymes and want to use these cut site correction functions, you'll need to adjust the referenced enzyme recognition sequences. Alternatively you can omit the `cut_correct` function and make some adjustments to variable names to design sequences without cut site correction.
