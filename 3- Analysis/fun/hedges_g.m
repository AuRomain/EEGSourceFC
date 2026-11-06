function g = hedges_g(x, y)
%HEDGES_G Hedges' g of two independent samples: (mean(x) - mean(y)) over the pooled SD, bias-corrected.
%   g = hedges_g(x, y)
%   Cohen's d with the pooled SD (df = nx + ny - 2) times Hedges' exact
%   correction J(df) = Gamma(df/2) / (sqrt(df/2) Gamma((df - 1)/2))
%   (Hedges 1981, J Educ Stat 6:107), as PERMUTOOLS' booteffectsize
%   computes it with 'effect' 'cohen', 'vartype' 'equal', 'correct' true
%   (Crosse et al. 2024; booteffectsize.m lines 257, 274 and 380); before
%   5 November 2026 the approximation J = 1 - 3 / (4 df - 1), 0.0003
%   larger at df = 10. NaN values are left out; NaN when a sample has
%   fewer than 2 values or the pooled SD is 0. tests/test_reference_statistics
%   checks it against booteffectsize.
x = x(~isnan(x));  y = y(~isnan(y));
nx = numel(x);  ny = numel(y);
g = NaN;
if nx < 2 || ny < 2, return; end
df = nx + ny - 2;
sp = sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / df);
J = exp(gammaln(df / 2) - log(sqrt(df / 2)) - gammaln((df - 1) / 2));
g = J * (mean(x) - mean(y)) / sp;
if ~isfinite(g), g = NaN; end
end
