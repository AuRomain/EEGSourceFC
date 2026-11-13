function q = fdr_bh(p)
%FDR_BH Benjamini-Hochberg adjusted p-values (q-values); NaN entries kept.
%   q = fdr_bh(p): q_(i) = min over j >= i of min(1, m p_(j) / j), the p
%   sorted, m the number of p that are not NaN (Benjamini & Hochberg 1995).
%   The same as MATLAB's mafdr(p, 'BHFDR', true) (Bioinformatics Toolbox,
%   not needed here) and R's p.adjust(p, 'BH') on the p that are not NaN;
%   The tests check it against mafdr when that toolbox is present.
q = nan(size(p));
ok = ~isnan(p);
pv = p(ok);
m = numel(pv);
if m == 0, return; end
[ps, order] = sort(pv(:));
qs = ps .* m ./ (1:m)';
qs = min(1, flipud(cummin(flipud(qs))));
tmp = zeros(m, 1);
tmp(order) = qs;
q(ok) = tmp;
end
