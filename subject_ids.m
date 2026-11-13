function ids = subject_ids(S, names)
%SUBJECT_IDS Subject ID of each recording: the first token of S.id_regexp on its name.
%   ids = subject_ids(S, names) with names the file names without extension
%   (e.g. SUB_001_eo_ready -> SUB_001); the whole name when the expression
%   does not match or S.id_regexp is empty. Used by get_subject_groups and
%   by stage 1 (iaf.csv).
names = string(names(:));
ids = names;
if ~isempty(S.id_regexp)
    tok = regexp(cellstr(names), S.id_regexp, 'tokens', 'once');   % one cell per name, also for a single name
    for i = 1:numel(names)
        if ~isempty(tok{i}), ids(i) = string(tok{i}{1}); end
    end
end
end
