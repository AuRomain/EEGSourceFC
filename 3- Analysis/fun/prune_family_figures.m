function [kept, removed] = prune_family_figures(T, out_dir, alpha)
%PRUNE_FAMILY_FIGURES Remove the figures of the tests not significant after the adjustment over families.
%   [kept, removed] = prune_family_figures(T, out_dir, alpha), T a family's
%   table after family_adjust (columns p_maxT_family, significant) and
%   out_dir the folder it was written to: permu_ttest_permu_corr_BNS drew a
%   figure in out_dir/figures for every test with p_maxT_family <= alpha
%   (the tests that can be significant); those whose test is not
%   significant after the selection of families are deleted, so that the
%   figures left are the significant tests (family_figure_name gives the
%   files). Called by run_statistics.
kept = 0;  removed = 0;
if isempty(T) || height(T) == 0 || ~all(ismember({'p_maxT_family', 'significant'}, T.Properties.VariableNames))
    return
end
for r = find(T.p_maxT_family <= alpha)'
    f = fullfile(out_dir, 'figures', family_figure_name(T, r));
    if ~isfile(f), continue; end
    if T.significant(r)
        kept = kept + 1;
    else
        delete(f);
        removed = removed + 1;
    end
end
end
