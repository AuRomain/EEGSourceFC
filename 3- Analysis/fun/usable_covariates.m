function [Z, names, note] = usable_covariates(Z, names, what, group)
%USABLE_COVARIATES The covariates that can enter the model for these subjects.
%   [Z, names] = usable_covariates(Z, names, what, group) leaves out, one
%   by one, covariates
%     - missing for all, or constant among the subjects that have them;
%     - (when group is given) present in fewer than 2 subjects of either
%       group (e.g. a medication dose that controls do not have, which would
%       leave every test with N = 0);
%     - (when group is given) with a single value in one group that fewer
%       than 2 subjects of the other group share (e.g. a dose of 0 in every
%       control and above 0 in every patient): the group effect would rest
%       on a straight-line extrapolation between groups that never overlap,
%       and power is lost (in a simulation 0.84 without the dose, 0.20 with
%       it). A covariate with one value in a group that the other group also
%       has (all controls female, both sexes among the patients) is kept;
%     - collinear with the intercept, the group and the covariates kept
%       before them,
%   with a warning naming WHAT (the comparison or correlation set) and the
%   reason. Used by the group comparisons, the correlations, NBS and the
%   sample-size check; the names kept are written with every result.
%   [Z, names, note] = usable_covariates(...) returns the warning text in
%   NOTE ('' when nothing is left out) instead of issuing it.
drop = false(1, size(Z, 2));
why = strings(1, size(Z, 2));
for c = 1:size(Z, 2)
    v = Z(~isnan(Z(:, c)), c);
    drop(c) = numel(v) < 2 || all(v == v(1));
    if drop(c), why(c) = "missing or constant"; continue; end
    if isempty(group), continue; end
    has = ~isnan(Z(:, c));
    inB = has & group(:) ~= 0;  inA = has & group(:) == 0;
    if nnz(inB) < 2 || nnz(inA) < 2                % both groups need it
        drop(c) = true;  why(c) = "missing in a group";
    elseif one_value_apart(Z(inA, c), Z(inB, c)) || one_value_apart(Z(inB, c), Z(inA, c))
        drop(c) = true;  why(c) = "one value in a group, which the other group (almost) never has";
    end
end
keep = find(~drop);
base = ones(size(Z, 1), 1);
if ~isempty(group), base = [base, group(:)]; end
okr = all(~isnan(Z(:, keep)), 2);
X = base(okr, :);
used = [];
for c = keep
    okc = okr & ~isnan(Z(:, c));
    if rank([X(okc(okr), :), Z(okc, c)]) > rank(X(okc(okr), :))
        X = [X, Z(okr, c)]; %#ok<AGROW>
        used(end + 1) = c; %#ok<AGROW>
    else
        drop(c) = true;  why(c) = "collinear";
    end
end
note = '';
if any(drop)
    r = unique(why(drop), 'stable');
    parts = strings(1, numel(r));
    for i = 1:numel(r)
        parts(i) = strjoin(string(names(drop & why == r(i))), ', ') + " (" + r(i) + ")";
    end
    note = sprintf('%s: covariate(s) left out in these subjects: %s', what, strjoin(parts, '; '));
    if nargout < 3
        warning('permu:covariates', '%s', note);
    end
end
Z = Z(:, used);
names = names(used);
end

function tf = one_value_apart(a, b)
% A has a single value and fewer than 2 values of B equal it
tf = all(a == a(1)) && nnz(b == a(1)) < 2;
end
