function z = t_to_z(t, df)
%T_TO_Z A t (or Welch-type G) statistic on df degrees of freedom as the z of the same tail probability.
%   z = t_to_z(t, df)
%   t, df   arrays of the same size (or df scalar); df Inf gives z = t
%   z       sign(t) * the normal quantile of t's one-tailed probability;
%           NaN where t or df is NaN
%   PALM's palm_gtoz (Winkler et al. 2014, NeuroImage 92:381;
%   github.com/andersonwinkler/PALM, palm_gtoz.m, the branch for one
%   numerator degree of freedom), each tail on its own for precision:
%   z = erfcinv(2 tcdf(-|t|, df)) sqrt(2), with the sign of t (palm_gcdf is
%   Student's t cdf there). Statistics of different degrees of freedom are
%   on one scale once converted, so their maximum over tests (or over
%   permutations of a Welch-type statistic, whose df changes with the data)
%   is valid (Alberton et al. 2020, NeuroImage 216:116760).
%   tests/test_reference_statistics checks it against palm_gtoz.
if isscalar(df), df = repmat(df, size(t)); end
z = sign(t) .* erfcinv(2 * tcdf(-abs(t), df)) * sqrt(2);
end
