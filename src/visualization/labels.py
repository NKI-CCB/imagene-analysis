feature_order = [
    'volume', 'largest_diameter',
    'vol_init_enhancement_GT100', 'ld_init_enhancement_GT100',
    'vol_late_LT0', 'ld_late_LT0',
    'mean_smoothness_uptake', 'variation_smoothness_uptake',
    'mean_smoothness_all_timeframes', 'variation_smoothness_all_timeframes',
    'mean_sharpness_uptake', 'var_sharpness_uptake',
    'mean_sharpness_all_timeframes', 'var_sharpness_all_timeframes',
    'uptake_speed',
    'top_init_enhancement', 'top_late_enhancement',
    'ser', 'washout',
    'circularity',
    'irregularity',
]

feature_display_names = {
    'volume': "Volume",
    'largest_diameter': "Diameter",
    'vol_init_enhancement_GT100': "Volume Initial Enhancement > 100",
    'ld_init_enhancement_GT100': "Diameter Initial Enhancement > 100",
    'vol_late_LT0': "Volume Late Enhancement < 0",
    'ld_late_LT0': "Diameter Late Enhancement < 0",
    'mean_smoothness_uptake': "Smoothness Uptake (mean)",
    'variation_smoothness_uptake': "Smoothness Uptake (variation)",
    'mean_smoothness_all_timeframes': "Smoothness Maximum (mean)",
    'variation_smoothness_all_timeframes': "Smoothness Maximum (variation)",
    'mean_sharpness_uptake': "Sharpness Uptake (mean)",
    'var_sharpness_uptake': "Sharpness Uptake (variation)",
    'mean_sharpness_all_timeframes': "Sharpness Maximum (mean)",
    'var_sharpness_all_timeframes': "Sharpness Maximum (variation)",
    'uptake_speed': "Uptake Speed",
    'top_init_enhancement': "Top Initial Enhancement",
    'top_late_enhancement': "Top Late Enhancement",
    'ser': "Signal Enhancement Ratio",
    'washout': "Washout",
    'mean_vox_val': "Average Pre-Contrast Voxel Value",
    'variance_vox_val': "Variance Pre-Contrast Voxel Value",
    'circularity': "Circularity",
    'irregularity': "Irregularity",
}

factor_display_names = {
    'size': "Tumor Size",
    'shape': "Tumor Shape",
    'smoothness': "Smoothness",
    'sharpness_var': "Sharpness (variation)",
    'late_enhancement': "Late Enhancement",
    'init_enhancement': "Initial Enhancement",
    'circularity': "Circularity",
    'sharpness_mean': "Sharpness (mean)",
    'irregularity': "Irregularity",
}

clin_display_names = {
    'grade': "Grade",
    'ihc_subtype': "IHC Subtype",
    'age_at_diagnosis': "Age",
}
