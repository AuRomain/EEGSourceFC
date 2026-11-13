function names = spectral_profile_measures()
%SPECTRAL_PROFILE_MEASURES Names of the numeric measures of the spectral profile (spectral_profile.csv).
%   names = spectral_profile_measures(): the columns that
%   spectral_profile_summary computes and read_subject_table joins to the
%   subject table with S.spectral_profile (IAF, IAF_cov, alpha_pw_post,
%   exponent_all, ...; IAF_source and IAF_flag, text, are not listed). Used
%   to recognise a covariate that comes from the spectral profile.
names = {'IAF', 'IAF_PAF', 'IAF_CoG', 'IAF_cov', 'IAF_IRASA', 'IAF_specparam', 'alpha_pw_post', 'alpha_peak_share', ...
    'alpha_cf_post', 'alpha_front_post', 'exponent_all', 'offset_all', 'exponent_post', 'fit_r2', ...
    'log_alpha_post', 'rel_alpha_post', 'theta_alpha_post', 'purity_fixed', 'purity_iaf', 'capture_fixed', ...
    'capture_iaf', 'hf_temporal', 'eye_frontopolar'};
end
