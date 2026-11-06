function [Scouts, C] = get_atlas(atlas, A, fields)
%GET_ATLAS The regions (Brainstorm scouts) of an atlas on the cortex of a head.
%   [Scouts, C] = get_atlas(atlas) on the shipped ICBM152 cortex
%   (inputs/tess_cortex_pial_low.mat, 15,002 vertices);
%   [Scouts, C] = get_atlas(atlas, A) on the cortex of the head A (get_anatomy).
%   atlas:
%     'Desikan_Killiany_68'  shipped head: inputs/scout_Desikan-Killiany_68.mat;
%                            another head: the 'Desikan-Killiany' atlas of
%                            its cortex (Brainstorm names FreeSurfer's and
%                            CAT12's aparc / aparc_DK40 so)
%     a scout file           Brainstorm scouts (Scouts, or Atlas(1).Scouts of a
%                            surface-atlas export) defined on that cortex
%     an atlas name          one of the cortex's atlases (its Atlas field), e.g.
%                            'Destrieux', 'DKT', 'Schaefer_100_17net',
%                            'Brainnetome', 'HCP_MMP1' (the shipped cortex has 22)
%   On another head the regions are put in the order of the shipped
%   cortex's atlas of the same name, so every head gives the same rows:
%   a region of it missing on this head is an error, regions it does not
%   have (e.g. 'unknown') are left out; an atlas the shipped cortex does
%   not have is sorted by label.
%   C: the cortex's VertNormals, and the fields asked for (e.g.
%   {'Vertices', 'Faces', 'Reg'}).
if nargin < 2 || isempty(A), A = get_anatomy(struct()); end
if nargin < 3, fields = {}; end
w = whos('-file', A.cortex);
fields = intersect([{'VertNormals'}, cellstr(fields)], {w.name});   % those the file has
C = load(A.cortex, fields{:});
in = fullfile(fileparts(mfilename('fullpath')), 'inputs');
name = char(atlas);
if strcmp(name, 'Desikan_Killiany_68'), name = 'Desikan-Killiany'; end
if isfile(atlas)
    Scouts = scout_file(atlas);
elseif A.shipped && strcmp(atlas, 'Desikan_Killiany_68')
    Scouts = scout_file(fullfile(in, 'scout_Desikan-Killiany_68.mat'));
else
    Scouts = by_name(A.cortex, name, atlas);
    if ~A.shipped                                  % the regions in the shipped head's order
        try
            ref = by_name(fullfile(in, 'tess_cortex_pial_low.mat'), name, atlas);
            [ok, row] = ismember({ref.Label}, {Scouts.Label});
            if ~all(ok)
                error('get_atlas:regions', '%s: atlas %s lacks %d region(s) of the shipped ICBM152''s: %s', ...
                    A.cortex, name, nnz(~ok), strjoin({ref(~ok).Label}, ', '));
            end
            Scouts = Scouts(row);
        catch err
            if ~strcmp(err.identifier, 'get_atlas:atlas'), rethrow(err); end
            [~, order] = sort({Scouts.Label});
            Scouts = Scouts(order);
        end
    end
end
Scouts = Scouts(:)';
v = [Scouts.Vertices];
if max(v) > size(C.VertNormals, 1)
    error('get_atlas:vertices', '%s is not defined on %s (%d vertices)', atlas, A.cortex, size(C.VertNormals, 1));
end
if numel(unique(v)) < numel(v)
    warning('get_scouts_timeseries:overlap', '%s: some vertices belong to several scouts', atlas);
end
end

function Scouts = scout_file(file)
F = load(file);
if isfield(F, 'Scouts'), Scouts = F.Scouts;
elseif isfield(F, 'Atlas'), Scouts = F.Atlas(1).Scouts;      % Brainstorm surface-atlas export
else, error('get_atlas:file', '%s has no Scouts', file);
end
end

function Scouts = by_name(cortex, name, atlas)
F = load(cortex, 'Atlas');
names = {};
if isfield(F, 'Atlas'), names = {F.Atlas.Name}; end
i = find(strcmpi(names, name), 1);
if isempty(i) || isempty(F.Atlas(i).Scouts)
    error('get_atlas:atlas', ['Unknown atlas %s: give Desikan_Killiany_68, a scout file, or an atlas of ' ...
        '%s (%s)'], atlas, cortex, strjoin(names(~cellfun(@isempty, {F.Atlas.Scouts})), ', '));
end
Scouts = F.Atlas(i).Scouts;
end
