function prepare_spectral_covariates(S)
%PREPARE_SPECTRAL_COVARIATES Compute the spectral profile first when covariates need it.
%   prepare_spectral_covariates(S), called by the PSD statistics (stage 1b),
%   which run before stage 3: when a covariate of S (S.covariates, a
%   correlation set or a covariate set) is a measure of the spectral profile
%   (spectral_profile_measures: alpha_pw_post, exponent_all, IAF, ...) and
%   not a column of the subject table (nor of iaf.csv), the profile is
%   computed now (spectral_profile_features: its cache, which stage 3 then
%   reuses) and written to results/Connectivity/spectral_profile.csv, which
%   read_subject_table joins. On a fresh run the covariate was otherwise
%   missing, and every PSD test failed. With S.spectral_profile = false the
%   measures are never joined: an error says so.
c = cellstr(S.covariates);
for k = 1:numel(S.correlations), c = [c, cellstr(S.correlations(k).covariates)]; end %#ok<AGROW>
for k = 1:numel(S.covariate_sets), c = [c, cellstr(S.covariate_sets(k).covariates)]; end %#ok<AGROW>
need = intersect(unique(c), spectral_profile_measures());
if isempty(need) || isempty(S.subs_info) || ~isfile(S.subs_info), return; end
Sx = S;
Sx.spectral_profile = false;
T0 = read_subject_table(Sx);                       % the table (and iaf.csv with individual bands)
need = setdiff(need, T0.Properties.VariableNames);
if isempty(need), return; end
if ~S.spectral_profile
    error('prepare_spectral_covariates:off', ['Covariate(s) %s are measures of the spectral profile ' ...
        '(spectral_profile_features), not columns of the subject table, and S.spectral_profile is false: ' ...
        'set it to true, or leave them out'], strjoin(need, ', '));
end
fprintf('Covariate(s) %s come from the spectral profile: computed first (stage 3 reuses it)\n', strjoin(need, ', '));
spectral_profile_features(S);
end
