function Gr = psd_groups(S, T, Lk)
%PSD_GROUPS The two groups of each comparison among a level's recordings (names, members, colours).
%   Gr = psd_groups(S, T, Lk) with T the subject table (read_subject_table)
%   and Lk a level of the PSD step (its recordings' ids, Lk.ids): per
%   comparison of S.comparisons ("A - B"), its name, the two group names,
%   their members (recordings x 2, select_subjects of cmp.A and cmp.B) and
%   colours (Okabe-Ito hues, one per group name across comparisons).
%   Without a table or comparisons: all recordings as one group ("all").
%   Used by psd_figures, psd_region_check and psd_vertex_maps.
pal = [213 94 0; 0 114 178; 170 51 119; 0 158 115; 230 159 0; 86 180 233; 102 102 102] / 255;
n = numel(Lk.ids);
Gr = struct('name', {}, 'names', {}, 'in', {}, 'colours', {});
if ~isempty(T) && ~isempty(S.comparisons)
    [ok, row] = ismember(string(Lk.ids), string(T.(S.id_column)));
    seen = strings(0);
    for cmp = S.comparisons(:)'
        names = split(string(cmp.name), " - ");
        if numel(names) ~= 2, names = string(cmp.name) + [" (A)"; " (B)"]; end
        in = false(n, 2);
        in(ok, 1) = select_subjects(T(row(ok), :), cmp.A);
        in(ok, 2) = select_subjects(T(row(ok), :), cmp.B);
        col = zeros(2, 3);
        for j = 1:2
            c = find(seen == names(j), 1);
            if isempty(c), seen(end + 1) = names(j); c = numel(seen); end %#ok<AGROW>
            col(j, :) = pal(min(c, end), :);
        end
        Gr(end + 1) = struct('name', string(cmp.name), 'names', names(:)', 'in', in, 'colours', col); %#ok<AGROW>
    end
end
if isempty(Gr)
    Gr = struct('name', "all recordings", 'names', "all", 'in', true(n, 1), 'colours', [0.3 0.3 0.3]);
end
end
