function [head, lines, n_fig] = failure_summary(failures, minutes)
%FAILURE_SUMMARY  The last line of a run in which steps failed, and the list of those steps.
%   [head, lines, n_fig] = failure_summary(failures, minutes)
%   failures  string array (what, message), one row per step that failed,
%             as run_statistics and run_psd_analysis return them (not empty)
%   minutes   the time of the run
%   head      "Done with errors in ... min: N step(s) failed", with
%             "(M of them figures only)" when some of them only draw
%             figures; when every one of them only draws figures, "Done in
%             ... min: every result written, but M figure step(s) failed":
%             a figure that fails leaves out no table and no test
%   lines     "what: message" of each step, the steps that compute results
%             first, then the figure steps, marked "(figure only)"
%   n_fig     the number of figure steps
%   A step only draws figures when its name ends in "figure" or "figures"
%   (the group-matrix figure, the PSD figures, the state network figures).
fig = endsWith(failures(:, 1), ["figure", "figures"]);
n_fig = nnz(fig);
if all(fig)
    head = sprintf('Done in %.1f min: every result written, but %d figure step(s) failed', minutes, n_fig);
else
    head = sprintf('Done with errors in %.1f min: %d step(s) failed', minutes, size(failures, 1));
    if n_fig > 0, head = sprintf('%s (%d of them figures only)', head, n_fig); end
end
what = failures(:, 1);
what(fig) = what(fig) + " (figure only)";
lines = [what(~fig) + ": " + failures(~fig, 2); what(fig) + ": " + failures(fig, 2)];
end
