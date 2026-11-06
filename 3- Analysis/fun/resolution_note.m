function s = resolution_note(n_perm, alpha, F, label)
%RESOLUTION_NOTE What the number of permutations allows once families are selected (max-statistic, Benjamini-Bogomolov).
%   s = resolution_note(n_perm, alpha) states the rule, before the families
%   are known (start of stage 3).
%   s = resolution_note(n_perm, alpha, F) gives the numbers for F families
%   (end of stage 3); resolution_note(n_perm, alpha, F, label) adds a word
%   on F (start of the PSD statistics: "at most").
%   A permutation p is at least p_min = 1 / (n_perm + 1), and so are a
%   family's p (family_p, its largest statistic against the permutations'
%   largest) and a test's family-wise p (p_maxT_family), whatever the size
%   of the family. When one family has signal and the others none,
%   Benjamini-Hochberg over the F families selects it when its family_p <=
%   alpha / F, and its strongest test is then significant (p_bb =
%   p_maxT_family F / R <= alpha, with R = 1): a family is selected alone
%   when at most e = alpha (n_perm + 1) / F - 1 permutations reach its
%   largest statistic. With e < 0 no family can be selected (n_perm + 1 <
%   F / alpha); about 10 F / alpha permutations leave room (e = 9).
%   (family_adjust; docs/validation.md, "Permutation resolution".)
pmin = 1 / (n_perm + 1);
if nargin < 3
    s = sprintf(['Permutation resolution: the smallest p is 1/(n_perm + 1) = %.2g, whatever the size of a ' ...
        'family. Among F families, a family is selected, and its strongest test significant, only when its p is ' ...
        'at most alpha / F: n_perm >= F / alpha - 1 at the least, about 10 F / alpha to leave room. The numbers ' ...
        'of this run are given at the end.'], pmin);
    return
end
if nargin < 4, label = ""; end
if strlength(string(label)) > 0, label = " (" + string(label) + ")"; end
need = ceil(F / alpha - 1);
room = ceil(10 * F / alpha - 1);
e = floor(alpha * (n_perm + 1) / F - 1 + 1e-9);
if e < 0
    s = sprintf(['Permutation resolution: %d families%s; the smallest p (%.2g) is above alpha / F = %.2g: no ' ...
        'family can be selected. n_perm >= %d is needed, about %d to leave room.'], F, label, pmin, alpha / F, need, room);
    return
end
s = sprintf(['Permutation resolution: %d families%s; a family is selected alone, and its strongest test ' ...
    'significant, when at most %d of the %d permutations reach its largest statistic (p <= alpha / F = %.2g)'], ...
    F, label, e, n_perm, alpha / F);
if e < 9
    s = [s sprintf('; n_perm = %d would allow 9', room)];
end
s = [s '.'];
end
