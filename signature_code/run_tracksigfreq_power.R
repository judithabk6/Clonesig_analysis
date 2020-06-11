#!/usr/bin/env Rscript

# /bioinfo/local/build/R/R-3.3.2_centos/bin/R
library(optparse)
library(TrackSig)
library(ggplot2)

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
#folder_path = '20200210_simulations_clonesig_cn_cancer_type/type5-perc_diploid20-nb_clones2-nb_mut1000'
#folder_path = '20200210_simulations_clonesig_cn_cancer_type/type12-perc_diploid20-nb_clones1-nb_mut1000'
vcaf_file  <- paste(folder_path, 'tracksigfreq', 'vcaf.csv', sep='/')
trinut_count_file <- paste(folder_path, 'tracksigfreq', 'batch_100_pattern96.csv', sep='/')

trinut_count_data <- tryCatch(read.table(trinut_count_file, header=TRUE, row.names=1), error=function(e) NULL)
vcaf_full <- read.table(vcaf_file, header=TRUE)


scoreMethod = "SigFreq"
binSize = 100
desiredMinSegLen = NULL

MU_table = read.table("external_data/sigProfiler_SBS_signatures_2018_03_28.csv", sep=',', header=TRUE)
MU = MU_table[, 3:dim(MU_table)[2]]

# get cancer signature match
cancer_type_sig = read.table("external_data/curated_match_signature_cancertype_tcgawes_literature.csv", sep='\t', header=TRUE, row.names=1)
bigselect = rowSums(cancer_type_sig)
MU = MU[,bigselect>=1]

row.names(MU) <- row.names(trinut_count_data)

vcaf = vcaf_full[,c('phi', 'qi', 'bin')]

#results <- TrackSig:::getChangepointsPELT(trinut_count_data, subMU, vcaf, scoreMethod, binSize, desiredMinSegLen)



suffix = '_all'
start_time <- Sys.time()
results <- TrackSig:::getChangepointsPELT(trinut_count_data, MU, vcaf, scoreMethod, binSize, desiredMinSegLen)
end_time <- Sys.time()
runtime = end_time - start_time
write.csv(results$mixtures, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_mixtures", suffix, ".csv", sep=''), sep='/'))
n_col <- ifelse(length(results$changepoints) > 0, length(results$changepoints), 1)
write(results$changepoints, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_changepoints", suffix, ".txt", sep=''), sep='/'), ncolumns=n_col)
write.csv(runtime, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_runtime", suffix, ".csv", sep=''), sep='/'))

suffix = '_prefit'
start_time <- Sys.time()
premixtures <- TrackSig:::fitMixturesEM(rowSums(trinut_count_data), MU)
prefit_mu = MU[,premixtures>0.05]
results <- TrackSig:::getChangepointsPELT(trinut_count_data, prefit_mu, vcaf, scoreMethod, binSize, desiredMinSegLen)
end_time <- Sys.time()
runtime = end_time - start_time
write.csv(results$mixtures, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_mixtures", suffix, ".csv", sep=''), sep='/'))
write.csv(premixtures, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_premixtures", suffix, ".csv", sep=''), sep='/'))
n_col <- ifelse(length(results$changepoints) > 0, length(results$changepoints), 1)
write(results$changepoints, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_changepoints", suffix, ".txt", sep=''), sep='/'), ncolumns=n_col)
write.csv(runtime, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_runtime", suffix, ".csv", sep=''), sep='/'))


