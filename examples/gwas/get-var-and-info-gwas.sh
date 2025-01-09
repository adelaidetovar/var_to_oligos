#!/bin/bash
#SBATCH --time="24:00:00"
#SBATCH --mem=10G
#SBATCH --output=logs/slurm-%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type END,FAIL
#SBATCH --signal=B:TERM@60
#SBATCH --job-name=var-to-oligos
#SBATCH --account=[acct]
#SBATCH --mail-user [username]@umich.edu

Rscript /path/to/var_to_oligos/scripts/get_var_and_info.R -intab signal_list.txt -ldlinktoken [ldlinktoken] -ncbikey [ncbikey] -outdir var_out_gwas
