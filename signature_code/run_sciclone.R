#!/usr/bin/env Rscript

# /bioinfo/local/build/R/R-3.3.2_centos/bin/R
library(sciClone)
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


cnv_file = read.table(paste(folder_path, 'cnv_table.csv', sep='/'), sep='\t', header=TRUE)
cnv_file$total_cn = cnv_file$major + cnv_file$minor
if ("chromosome" %in% colnames(cnv_file)) {
    cnv_to_use = cnv_file[,c('chromosome', 'start', 'end', 'total_cn')]
} else {
    cnv_file$chromosome = cnv_file$chr
    cnv_file$start = cnv_file$startpos
    cnv_file$end = cnv_file$endpos
    cnv_file$total_cn = cnv_file$major_cn + cnv_file$minor_cn
    cnv_to_use = cnv_file[,c('chromosome', 'start', 'end', 'total_cn')]
}
tt = read.table(paste(folder_path, 'input_t.tsv', sep='/'), sep='\t', header=TRUE, comment.char='@')
tt$vaf = tt$var_counts/(tt$var_counts + tt$ref_counts) * 100
gg = tt[,c('chromosome', 'position', 'ref_counts', 'var_counts', 'vaf')]
gg$position = as.numeric(as.character(gg$position))
gg$chromosome = gsub('chr', '', gg$chromosome)

sc = sciClone(vafs=list(gg), copyNumberCalls=list(cnv_to_use), sampleNames=c('tumor'), minimumDepth=3, clusterMethod="bmm", clusterParams="no.apply.overlapping.std.dev.condition", cnCallsAreLog2=FALSE, useSexChrs=TRUE, doClustering=TRUE, verbose=TRUE,  copyNumberMargins=0.25, maximumClusters=10, annotation=NULL, doClusteringAlongMargins=TRUE, plotIntermediateResults=0)
out_path = paste(folder_path, 'sciclone', sep='/')
dir.create(out_path, showWarnings = FALSE,recursive=TRUE)
writeClusterTable(sc, paste(out_path, 'clusters1', sep='/'))
sc.plot1d(sc,paste(out_path, "clusters1_1d.pdf", sep='/'))
