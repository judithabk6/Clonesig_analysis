#!/usr/bin/env Rscript

# /bioinfo/local/build/R/R-3.3.2_centos/bin/R
library(optparse)
library(deconstructSigs)

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

# get cancer signature match
subMU_table = read.table(paste(folder_path, 'subMU.csv', sep='/'), sep='\t', header=T)
subMU = data.frame(t(subMU_table[, 2:dim(subMU_table)[2]]))
colnames(subMU) = colnames(signatures.cosmic)

tt = read.table(paste(folder_path, 'deconstructsigs', 'pattern96.csv', sep='/'), sep='\t', header=TRUE, check.names=FALSE)
rownames(tt)=c(1)
colnames(tt) = colnames(signatures.cosmic)
start_time <- Sys.time()
test = whichSignatures(tumor.ref = tt, 
                       signatures.ref = subMU, 
                       tri.counts.method = 'default',
                       contexts.needed=TRUE)
end_time <- Sys.time()
runtime = end_time - start_time
write.table(test$weights, paste(folder_path, 'deconstructsigs', 'signatures_cancertype.csv', sep='/'), quote=FALSE, row.names=FALSE)
write.csv(runtime, file=paste(folder_path, 'deconstructsigs', "deconstructsig_runtime_cancertype.csv", sep='/'))

