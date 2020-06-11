#!/usr/bin/env Rscript

# /bioinfo/local/build/R/R-3.3.2_centos/bin/R
library(optparse)
# load palimpsest
library(Palimpsest)
options(stringsAsFactors = FALSE)


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
#folder_path = '20190623_simulations_clonesig_cn_cancer_type/type5-perc_diploid20-nb_clones2-nb_mut300'


# get MU
MU_table = read.table("external_data/sigProfiler_SBS_signatures_2018_03_28.csv", sep=',', header=TRUE)
MU = MU_table[, 3:dim(MU_table)[2]]

# get cancer signature match
cancer_type_sig = read.table("external_data/curated_match_signature_cancertype_tcgawes_literature.csv", sep='\t', header=TRUE, row.names=1)
bigselect = rowSums(cancer_type_sig)
MU = MU[,bigselect>=1]

vcf_table = read.table(paste(folder_path, 'palimpsest', 'vcf_table.csv', sep='/'), sep='\t', header=TRUE)
cna_table = read.table(paste(folder_path, 'palimpsest', 'cna_table.csv', sep='/'), sep='\t', header=TRUE)
annot_table = read.table(paste(folder_path, 'palimpsest', 'annot_table.csv', sep='/'), sep='\t', header=TRUE)

propMutsByCat <- palimpsestInput(vcf = vcf_table,type = "SNV",sample.col = "Sample", mutcat.col = "mutcat3", proportion = TRUE)
propMutsByCat = as.matrix(propMutsByCat)
colnames(propMutsByCat) = c(unique(vcf_table$Sample)[1])
vcf_ccf <- cnaCCF_annot(vcf=vcf_table,annot_data = annot_table,cna_data = cna_table,CCF_boundary = 0.95)
vcf.clonal <- vcf_ccf[which(vcf_ccf$Clonality=="clonal"),]
vcf.subclonal <- vcf_ccf[which(vcf_ccf$Clonality=="subclonal"),]
propMutsByCat.clonal <- palimpsestInput(vcf = vcf.clonal,type="SNV",
                                        sample.col = "Sample",
                                        mutcat.col = "mutcat3",
                                        proportion = TRUE)
propMutsByCat.clonal = as.matrix(propMutsByCat.clonal)
colnames(propMutsByCat.clonal) = c(unique(vcf.clonal$Sample)[1])
propMutsByCat.subclonal <- palimpsestInput(vcf = vcf.subclonal,type="SNV",
                                        sample.col = "Sample",
                                        mutcat.col = "mutcat3",
                                        proportion = TRUE)
propMutsByCat.subclonal = as.matrix(propMutsByCat.subclonal)
colnames(propMutsByCat.subclonal) = c(unique(vcf.subclonal$Sample)[1])

for (suffix in c('prefit', 'all')) {
    start_time <- Sys.time()
    if (suffix == 'prefit') {
        signatures_exp <- deconvolution_fit(vcf = vcf_ccf, type = "SNV",
                                           input_data = propMutsByCat,
                                           threshold = 6,
                                           input_signatures = t(MU),
                                           sig_cols = mycol,plot = F,
                                           resdir = resdir)
        mu_to_use = t(MU)[signatures_exp$sig_props>0,]
    } else if (suffix == 'all') {
        mu_to_use = t(MU)

    } else {
        mu_to_use = t(subMU)

    }
    signatures_exp_clonal <- deconvolution_fit(vcf = vcf.clonal, type = "SNV",
                                               input_data = propMutsByCat.clonal,
                                               threshold = 6,
                                               input_signatures = mu_to_use,
                                               sig_cols = mycol,plot = F,
                                               resdir = resdir)

    signatures_exp_subclonal <- deconvolution_fit(vcf = vcf.subclonal, type = "SNV",
                                               input_data = propMutsByCat.subclonal,
                                               threshold = 6,
                                               input_signatures = mu_to_use,
                                               sig_cols = mycol,plot = F,
                                               resdir = resdir)
    vcf.subclonal.origin <- palimpsestOrigin(vcf.subclonal,type = "SNV",
                                             sample.col = "Sample",
                                             mutcat.col = "mutcat3",
                                             signature_contribution = signatures_exp_subclonal$sig_nums,
                                             input_signatures = mu_to_use)
    vcf.clonal.origin <- palimpsestOrigin(vcf.clonal,type = "SNV",
                                          sample.col = "Sample",
                                          mutcat.col = "mutcat3",
                                          signature_contribution = signatures_exp_clonal$sig_nums,
                                          input_signatures = mu_to_use)

    vcf_ccf$origin = ''
    vcf_ccf[which(vcf_ccf$Clonality=="subclonal"), 'origin'] = vcf.subclonal.origin$Sig.max
    vcf_ccf[which(vcf_ccf$Clonality=="clonal"), 'origin'] = vcf.clonal.origin$Sig.max

    end_time <- Sys.time()
    runtime = end_time - start_time

    mixtures = rbind(signatures_exp_clonal$sig_props, signatures_exp_subclonal$sig_props)
    rownames(mixtures) = c('clonal', 'subclonal')
    write.table(mixtures, file=paste(folder_path, 'palimpsest', paste("palimpsest_mixtures_", suffix, ".csv", sep=''), sep='/'), quote = FALSE, sep = "\t")
    write.table(vcf_ccf, file=paste(folder_path, 'palimpsest', paste("palimpsest_mut_data_", suffix, ".csv", sep=''), sep='/'), quote = FALSE, sep = "\t")
    write.csv(runtime, file=paste(folder_path, 'palimpsest', paste("palimpsest_runtime_", suffix, ".csv", sep=''), sep='/'))

    if (suffix == 'prefit') {
        n_col <- dim(mu_to_use)[1]
        write(rownames(mu_to_use), file=paste(folder_path, 'palimpsest', paste("palimpsest_premixtures_", suffix, ".txt", sep=''), sep='/'), ncolumns=n_col)
    }
}