#!/usr/bin/env Rscript

# /bioinfo/local/build/R/R-3.3.2_centos/bin/R
library(optparse)
# load tracksig
library(reshape)
library(ggplot2)
src_files <- setdiff(grep(".*R$", list.files(paste0( "/data/users/jabecass/dl_tools_centos/TrackSig/src"),full.names = T), value = T), 
                     c(paste0( "/data/users/jabecass/dl_tools_centos/TrackSig/src/compute_mutational_signatures.R"),
                       paste0( "/data/users/jabecass/dl_tools_centos/TrackSig/src/header.R")))
for (file in src_files)
{
  source(file)
}

option_list <- list(
    make_option(c("-f", "--folder"), type="character", default=NULL, 
              help="path to folder with data input", metavar="character")
); 
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
if (is.null(opt$folder)){
  print_help(opt_parser)
  stop("The folder name is a mandatory argument.\n", call.=FALSE)
}

folder_path <- opt$folder
#folder_path = '20190623_simulations_clonesig_cn_cancer_type/type5-perc_diploid2000-nb_clones2-nb_mut300'


# get MU
MU_table = read.table(paste(folder_path, 'MU_matrix.csv', sep='/'), sep='\t', header=T)
MU = t(MU_table[, 2:dim(MU_table)[2]])


# get cancer signature match
subMU_table = read.table(paste(folder_path, 'sub_MU_matrix.csv', sep='/'), sep='\t', header=T)
subMU = t(subMU_table[, 2:dim(subMU_table)[2]])

# load data (inspired from the tracksig package)
vcfFile <- paste(folder_path, 'tracksig', 'batch_100_pattern96.csv', sep='/')
vcfData <- tryCatch(read.table(vcfFile), error=function(e) NULL) # 96 trinucleotide counts are read as input
vcfData[,1] <- NULL # Filenames on the first column are deleted
if (sum(vcfData[,1] %% 1  != 0) > 0)
{
# second column represents phi values
phis <- vcfData[,1]
vcfData[,1] <- NULL
}
# Hack because of the previous bug in the code that assigned 101 mutations to each time points
rows_keep <- (apply(vcfData, 1, sum) == 101) | (apply(vcfData, 1, sum) == 100)
vcfData <- vcfData[rows_keep, ]
if (!is.null(phis))
{
phis <- phis[rows_keep]
}
vcf <- t(vcfData)
phis_for_plot <- phis_sliding_window <- phis
colnames(vcf) <- round(phis_sliding_window, 3)
assigns_phylo_nodes_sw <- assigns_phylo_nodes <- NULL

suffix = '_cancertype'
start_time <- Sys.time()
list[changepoints, mixtures] <- find_changepoints_pelt(vcf, t(subMU))
end_time <- Sys.time()
runtime = end_time - start_time
write.csv(mixtures, file=paste(folder_path, 'tracksig', paste("tracksig_mixtures", suffix, ".csv", sep=''), sep='/'))
n_col <- ifelse(length(changepoints) > 0, length(changepoints), 1)
write(changepoints, file=paste(folder_path, 'tracksig', paste("tracksig_changepoints", suffix, ".txt", sep=''), sep='/'), ncolumns=n_col)
write.csv(runtime, file=paste(folder_path, 'tracksig', paste("tracksig_runtime", suffix, ".csv", sep=''), sep='/'))

suffix = '_all'
start_time <- Sys.time()
list[changepoints, mixtures] <- find_changepoints_pelt(vcf, t(MU))
end_time <- Sys.time()
runtime = end_time - start_time
write.csv(mixtures, file=paste(folder_path, 'tracksig', paste("tracksig_mixtures", suffix, ".csv", sep=''), sep='/'))
n_col <- ifelse(length(changepoints) > 0, length(changepoints), 1)
write(changepoints, file=paste(folder_path, 'tracksig', paste("tracksig_changepoints", suffix, ".txt", sep=''), sep='/'), ncolumns=n_col)
write.csv(runtime, file=paste(folder_path, 'tracksig', paste("tracksig_runtime", suffix, ".csv", sep=''), sep='/'))

suffix = '_prefit'
start_time <- Sys.time()
premixtures <- fit_mixture_of_multinomials_EM(rowSums(vcf), t(MU))
prefit_mu = t(MU)[,premixtures>0.05]
list[changepoints, mixtures] <- find_changepoints_pelt(vcf, prefit_mu)
end_time <- Sys.time()
runtime = end_time - start_time
write.csv(mixtures, file=paste(folder_path, 'tracksig', paste("tracksig_mixtures", suffix, ".csv", sep=''), sep='/'))
write.csv(premixtures, file=paste(folder_path, 'tracksig', paste("tracksig_premixtures", suffix, ".csv", sep=''), sep='/'))
n_col <- ifelse(length(changepoints) > 0, length(changepoints), 1)
write(changepoints, file=paste(folder_path, 'tracksig', paste("tracksig_changepoints", suffix, ".txt", sep=''), sep='/'), ncolumns=n_col)
write.csv(runtime, file=paste(folder_path, 'tracksig', paste("tracksig_runtime", suffix, ".csv", sep=''), sep='/'))



