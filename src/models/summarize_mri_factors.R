library(dplyr)
library(ncdf4)
library(readr)
library(readxl)
library(tidyr)

mri1 <- read_xlsx('data/raw/mri-features.xlsx') %>%
    rename(margins_patient=MARGINSstudyNr) %>%
    select(-Comment, -MultiFocal, -PCE_top10percent, -mean_vox_val, -variance_vox_val)

mri2 <- read_xlsx('data/raw/Features_not_sequenced.xlsx') %>%
    select(-NewPatID, -in_study_tycho, -var_sharpness) %>%
    rename(top_init_enhancement=top_init_enhancment, top_late_enhancement=top_late_enhancment,
           vol_init_enhancement_GT100=vol_init_enhancment_GT100,
           ld_init_enhancement_GT100=ld_init_enhancment_GT100, mean_sharpness_uptake=mean_sharpness,
           mean_smoothness_all_timeframes=rad_grad_ind_all_timeframes,
           variation_smoothness_all_timeframes=std_rgh_val_all_timeframes,
           mean_smoothness_uptake=mean_smoothness, variation_smoothness_uptake=variation_smoothness,
           var_sharpness_uptake=variation_sharpness)

mri <- bind_rows(mri1, mri2) %>% filter_all(~ !is.na(.))

mri_scaled <- mri %>%
    mutate_at(vars(starts_with('vol')), ~ .^(1/3)) %>%
    mutate_at(vars(starts_with('var')), sqrt)

nc_read_matrix <- function(ds, var_name) {
    mat = ncvar_get(ds, var_name)

    row_dim <- ds[['var']][[var_name]][['dim']][[1]]
    col_dim <- ds[['var']][[var_name]][['dim']][[2]]
    mat_dim <- list(c(row_dim[['vals']]), c(col_dim[['vals']]))
    names(mat_dim) <- c(row_dim[['name']], col_dim[['name']])
    dimnames(mat) <- mat_dim
    mat
}

f <- nc_open('data/processed/mri-features-all-fa.nc')
mri_loadings <- nc_read_matrix(f, 'loadings')
nc_close(f)


f = nc_open('data/processed/gene-expression.nc')
gexp_cases <- ncvar_get(f, 'case')
nc_close(f)

summarize_mri <- function(mri) {
     mri %>%
         select(-margins_patient) %>%
         summarize_all(list(mean=mean,
                            sd=sd,
                            q1=~quantile(., .25),
                            median=median,
                            q3=~quantile(., .75))) %>%
         pivot_longer(everything(),
                      names_to=c('variable', 'statistic'),
                      names_pattern='(.*)_(.*)') %>%
         pivot_wider(names_from='statistic', values_from='value')
}

mri_summary_all <- summarize_mri(mrif)
mri_summary_gexp <- summarize_mri(filter(mrif, margins_patient %in% gexp_cases))

write_tsv(mri_summary_all, 'models/mrif_summary_all.tsv')
write_tsv(mri_summary_gexp, 'models/mrif_summary_gexp.tsv')
