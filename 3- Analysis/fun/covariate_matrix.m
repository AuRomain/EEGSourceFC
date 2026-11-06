function [Z, names] = covariate_matrix(T, cols, categorical)
%COVARIATE_MATRIX Covariates of the subject table as a numeric matrix.
%   [Z, names] = covariate_matrix(T, cols): numeric columns as they are;
%   text or categorical columns (e.g. sex) as 0/1 indicator columns, one
%   per level but the first (rows with a missing level get NaN). Rows with
%   any NaN are left out of each test by the test functions, and the number
%   of subjects of every test is reported. The original handled only numeric
%   covariates, through within-group regressions (S2).
%   [Z, names] = covariate_matrix(T, cols, categorical): the columns named
%   in CATEGORICAL (S.categorical_covariates) are indicator columns even
%   when numeric (e.g. a site coded 1, 2, 3: one indicator per value but the
%   smallest, named 'site=2', 'site=3'), instead of one linear term.
if nargin < 3, categorical = {}; end
Z = zeros(height(T), 0);
names = {};
for k = 1:numel(cols)
    c = cols{k};
    if ~ismember(c, T.Properties.VariableNames)
        hint = '';
        if ismember(c, spectral_profile_measures())
            hint = [': it is a measure of the spectral profile, joined from results/Connectivity/spectral_profile.csv ' ...
                'with S.spectral_profile = true once stage 3 or the PSD statistics have written it'];
        end
        error('covariate_matrix:column', 'Covariate %s is not in the subject table%s', c, hint);
    end
    v = T.(c);
    if (isnumeric(v) || islogical(v)) && ismember(c, categorical)
        v = double(v);
        lev = unique(v(~isnan(v)));
        for l = 2:numel(lev)
            z = double(v == lev(l));
            z(isnan(v)) = NaN;
            Z(:, end + 1) = z; %#ok<AGROW>
            names{end + 1} = sprintf('%s=%g', c, lev(l)); %#ok<AGROW>
        end
    elseif isnumeric(v) || islogical(v)
        Z(:, end + 1) = double(v); %#ok<AGROW>
        names{end + 1} = c; %#ok<AGROW>
    else
        v = lower(strtrim(string(v)));
        miss = ismissing(v) | v == "";
        lev = unique(v(~miss));
        for l = 2:numel(lev)
            z = double(v == lev(l));
            z(miss) = NaN;
            Z(:, end + 1) = z; %#ok<AGROW>
            names{end + 1} = sprintf('%s=%s', c, lev(l)); %#ok<AGROW>
        end
    end
end
end
