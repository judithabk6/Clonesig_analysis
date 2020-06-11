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
#folder_path = '20190623_simulations_clonesig_cn_cancer_type/type5-perc_diploid2000-nb_clones2-nb_mut300'

cancer_type = as.integer(strsplit(strsplit(strsplit(folder_path, '/')[[1]][2], 'type')[[1]][2], '-')[[1]][1])

# get MU
MU_table = read.table(paste(folder_path, 'MU_matrix.csv', sep='/'), sep='\t', header=T)
MU = data.frame(t(MU_table[, 2:dim(MU_table)[2]]))
colnames(MU) = colnames(signatures.cosmic)


# get cancer signature match
subMU_table = read.table(paste(folder_path, 'sub_MU_matrix.csv', sep='/'), sep='\t', header=T)
subMU = data.frame(t(subMU_table[, 2:dim(subMU_table)[2]]))
colnames(subMU) = colnames(signatures.cosmic)

tt = read.table(paste(folder_path, 'deconstructsigs', 'pattern96.csv', sep='/'), sep='\t', header=TRUE, check.names=FALSE)
rownames(tt)=c(1)
colnames(tt) = colnames(signatures.cosmic)
start_time <- Sys.time()
test = whichSignatures(tumor.ref = tt, 
                       signatures.ref = MU, 
                       tri.counts.method = 'default',
                       contexts.needed=TRUE)
end_time <- Sys.time()
runtime = end_time - start_time
write.table(test$weights, paste(folder_path, 'deconstructsigs', 'signatures_all.csv', sep='/'), quote=FALSE, row.names=FALSE)
write.csv(runtime, file=paste(folder_path, 'deconstructsigs', "deconstructsig_runtime_all.csv", sep='/'))


start_time <- Sys.time()
test = whichSignatures(tumor.ref = tt, 
                       signatures.ref = subMU, 
                       tri.counts.method = 'default',
                       contexts.needed=TRUE)
end_time <- Sys.time()
runtime = end_time - start_time
write.table(test$weights, paste(folder_path, 'deconstructsigs', 'signatures_cancertype.csv', sep='/'), quote=FALSE, row.names=FALSE)
write.csv(runtime, file=paste(folder_path, 'deconstructsigs', "deconstructsig_runtime_cancertype.csv", sep='/'))

