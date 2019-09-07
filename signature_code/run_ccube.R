#!/usr/bin/env Rscript

# /bioinfo/local/build/R/R-3.3.2_centos/bin/R
library(ccube)
library(optparse)

# get arguments
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

mydata = read.table(paste(folder_path, 'ccube', 'input.tsv', sep='/'), sep='\t', header=TRUE)

numOfClusterPool = 1:8
numOfRepeat = 1
results <- RunCcubePipeline(ssm = mydata, 
                            numOfClusterPool = numOfClusterPool, 
                            numOfRepeat = numOfRepeat,
                            runAnalysis = T, 
                            runQC = F)
write.table(results$ssm, paste(folder_path, 'ccube', 'ssm_clusters.csv', sep='/'),
    quote = FALSE, sep = "\t", row.names = FALSE)