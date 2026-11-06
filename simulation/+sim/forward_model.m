function fm = forward_model(S, leadfield_file)
%FORWARD_MODEL The forward model of stage 1, for simulating recordings.
%   fm = sim.forward_model(S): the cap S.montages{1} and atlas
%   S.atlases{1} of the settings, as stage 1 uses them:
%   fm.G        nChan x nVertices leadfield on the cortical normals
%               (tess_cortex_pial_low, ICBM152), average-referenced
%   fm.labels   channel labels of the cap
%   fm.src_pos, fm.src_nrm   vertex positions (m) and normals
%   fm.scouts   the atlas regions (Vertices, Seed, Label)
%   fm.network  Yeo-7 network of each region (S.network_table, rsn_yeo7)
%   fm = sim.forward_model(S, leadfield_file): the recordings are generated
%   with another leadfield (G_normal and labels, on the same cortex, e.g.
%   other electrode positions and skull conductivity: a forward model that
%   differs from the one stage 1 inverts), its channels in the cap's order.
in = fullfile(S.root, '1- EEG source connectivity', 'inputs');
C = load(fullfile(in, 'tess_cortex_pial_low.mat'), 'Vertices', 'VertNormals');
[~, ftLeadfield] = get_ftHeadmodel_and_ftLeadfield(S.montages{1}, S.atlases{1}, S);
ftElec = get_ftElec(S.montages{1});
nsrc = size(C.VertNormals, 1);
G = zeros(numel(ftElec.label), nsrc);
for i = 1:nsrc
    l = ftLeadfield.leadfield{i};
    if size(l, 2) == 3, l = l * C.VertNormals(i, :)'; end   % projected on the normal
    G(:, i) = l;
end
if nargin > 1 && ~isempty(leadfield_file)            % a mismatched forward model
    L = load(leadfield_file, 'G_normal', 'labels');
    [ok, row] = ismember(lower(ftElec.label(:)), lower(cellstr(L.labels(:))));
    if ~all(ok) || size(L.G_normal, 2) ~= nsrc
        error('forward_model:leadfield', '%s must have the channels of %s and %d sources', leadfield_file, S.montages{1}, nsrc);
    end
    G = double(L.G_normal(row, :));
end
fm.G = G - mean(G, 1);
fm.labels = ftElec.label(:)';
fm.src_pos = C.Vertices;
fm.src_nrm = C.VertNormals;
A = load(fullfile(in, 'scout_Desikan-Killiany_68.mat'));
if ~strcmp(S.atlases{1}, 'Desikan_Killiany_68'), A = load(S.atlases{1}); end
fm.scouts = A.Scouts;
T = readtable(S.network_table, 'TextType', 'string');
[ok, row] = ismember(lower(string({fm.scouts.Label})), lower(T.label));
if ~all(ok) || ~ismember('rsn_yeo7', T.Properties.VariableNames)
    error('forward_model:networks', 'S.network_table must give rsn_yeo7 for every region of the atlas');
end
fm.network = T.rsn_yeo7(row);
end
