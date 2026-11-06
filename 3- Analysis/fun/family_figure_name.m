function f = family_figure_name(T, i)
%FAMILY_FIGURE_NAME File name of the figure of row i of a family's table.
%   f = family_figure_name(T, i): T a group-comparison table (column
%   Comparison) or a correlation table (column CorrelationSet) of
%   permu_ttest_permu_corr_BNS, with its column Family:
%     <family>_<comparison>_<feature>_s<state>.png
%     <family>_<correlation set>_<feature>_<scale>_s<state>.png
%   (the correlation set in the name: a covariate set re-running a set of
%   correlations no longer overwrites its figures). Used to draw the figures
%   and, after the adjustment over families, to remove those of the tests
%   that are not significant (run_statistics).
v = @(s) matlab.lang.makeValidName(char(s));
if ismember('Comparison', T.Properties.VariableNames)
    f = sprintf('%s_%s_%s_s%d.png', T.Family(i), v(T.Comparison(i)), T.Parameter(i), T.state(i));
else
    f = sprintf('%s_%s_%s_%s_s%d.png', T.Family(i), v(T.CorrelationSet(i)), T.Parameter(i), v(T.Scale(i)), T.state(i));
end
end
