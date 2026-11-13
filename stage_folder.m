function [folder, n] = stage_folder(root, name)
%STAGE_FOLDER The folder of a stage of the pipeline and the stage's number.
%   [folder, n] = stage_folder(root, name): the folder "<n>- <name>" in
%   root (the pipeline folder, S.root) and n, the number of the stage; ''
%   and [] when this pipeline has no such stage. name:
%     'EEG source connectivity'  stage 1
%     'PSD analysis'             stage 2
%     'Dynamic connectivity'     the brain network states (stage 3; absent
%                                from a static-only pipeline)
%     'Analysis'                 the statistics (stage 4; stage 3 without
%                                the states)
%   The stages are numbered in their order, without gaps, so the code finds
%   them by name, with or without the states.
d = dir(fullfile(root, ['*- ' name]));
d = d([d.isdir]);
d = d(~cellfun(@isempty, regexp({d.name}, ['^\d+- ' regexptranslate('escape', name) '$'], 'once')));
if isempty(d)
    folder = '';  n = [];
    return
end
if numel(d) > 1
    error('stage_folder:many', 'More than one folder of the stage %s in %s: %s', name, root, ...
        strjoin({d.name}, ', '));
end
folder = fullfile(root, d.name);
n = sscanf(d.name, '%d', 1);
end
