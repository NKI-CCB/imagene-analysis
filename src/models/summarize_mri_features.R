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

mri_summary_all <- summarize_mri(mri)
mri_summary_gexp <- summarize_mri(filter(mri, margins_patient %in% gexp_cases))

write_tsv(mri_summary_all, 'models/mri_summary/mri_summary_all.tsv')
write_tsv(mri_summary_gexp, 'models/mri_summary/mri_summary_gexp.tsv')
