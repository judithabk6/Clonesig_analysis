library(survival)
library(xtable)


clonesig_res = read.table('20200519_tcga_surv.csv', sep='\t', header=TRUE)

bold <- function(x){paste0('{\\textbf{\\bfseries ', x, '}}')}

clonesig_res$nbClones = clonesig_res$nb_clones
clonesig_res$sigChange = clonesig_res$sig_change
for (loc in unique(clonesig_res$cancer_loc_x)) {
    sub = clonesig_res[clonesig_res$cancer_loc_x==loc,]
    mod1 = coxph(Surv(survival_months_15y, binary_vital_status_15y) ~ AGE + SEX + stage + group, data=sub)
    m = summary(mod1)
    cap = paste('Cox model for ', m$n, ' of the TCGA ', loc, ' cohort, with ', m$nevent, ' events, censured at 15 years.', sep='')
    df = m$coefficients
    print(xtable(df, caption=cap, type='latex', display=c('f', 'f', 'f', 'f', 'f', 'e')), sanitize.rownames.function = bold, sanitize.colnames.function = bold)
}

mod1 = coxph(Surv(survival_months_15y, binary_vital_status_15y) ~ AGE + SEX + stage + group, data=clonesig_res)
m = summary(mod1)
cap = paste('Cox model for ', m$n, ' of the protected TCGA cohort, with ', m$nevent, ' events, censured at 15 years.', sep='')
df = m$coefficients
print(xtable(df, caption=cap, type='latex', display=c('f', 'f', 'f', 'f', 'f', 'e')), sanitize.rownames.function = bold, sanitize.colnames.function = bold)
<<<<<<< HEAD



mod1 = coxph(Surv(survival_months_15y, binary_vital_status_15y) ~ group + strata(cancer_loc_x), data=clonesig_res)

=======
>>>>>>> update code clonesig review
