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
#folder_path = 'PhylogicNDT500/Sim_500_19_var'
vcaf_file  <- paste(folder_path, 'tracksigfreq', 'vcaf.csv', sep='/')
trinut_count_file <- paste(folder_path, 'tracksigfreq', 'batch_100_pattern96.csv', sep='/')

trinut_count_data <- tryCatch(read.table(trinut_count_file, header=TRUE, row.names=1), error=function(e) NULL)
vcaf_full <- read.table(vcaf_file, header=TRUE, comment.char = "%", sep='\t')


scoreMethod = "SigFreq"
binSize = 100
desiredMinSegLen = NULL


# get cancer signature match
subMU_table = read.table(paste(folder_path, 'subMU.csv', sep='/'), sep='\t', header=T)
subMU = subMU_table[, 2:dim(subMU_table)[2]]

row.names(subMU) <- row.names(trinut_count_data)

vcaf = vcaf_full[,c('phi', 'qi', 'bin')]

#results <- TrackSig:::getChangepointsPELT(trinut_count_data, subMU, vcaf, scoreMethod, binSize, desiredMinSegLen)


suffix = '_cancertype'
start_time <- Sys.time()
results <- TrackSig:::getChangepointsPELT(trinut_count_data, subMU, vcaf, scoreMethod, binSize, desiredMinSegLen)
end_time <- Sys.time()
runtime = end_time - start_time
write.csv(results$mixtures, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_mixtures", suffix, ".csv", sep=''), sep='/'))
n_col <- ifelse(length(results$changepoints) > 0, length(results$changepoints), 1)
write(results$changepoints, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_changepoints", suffix, ".txt", sep=''), sep='/'), ncolumns=n_col)
write.csv(runtime, file=paste(folder_path, 'tracksigfreq', paste("tracksigfreq_runtime", suffix, ".csv", sep=''), sep='/'))

