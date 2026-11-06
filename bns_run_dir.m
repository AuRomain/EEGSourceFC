function d = bns_run_dir(S)
%BNS_RUN_DIR Folder of one brain-network-state analysis (stages 2 and 3).
%   d = bns_run_dir(S): results/Connectivity/BNS analysis/<run>, where <run>
%   is S.run_name or, by default,
%   <method>_<min k>_<max k>_run_1_<restarts>_restarts_k-<rule>_<fingerprint>
%   (<rule> = fixed<k> with S.k_fixed)
%   (followed by _<cap>_<atlas> for another cap or atlas than the study's).
%   The fingerprint (6 characters) covers every stage 2 setting
%   (stage2_settings: smoothing, seed, surrogates, ...), so runs with other
%   settings get their own folder instead of overwriting each other (V3-6).
%   Stage 2 writes Results/ and "Results - smoothing/" there, stage 3 Analysis/.
%   Runs with other stage 1 settings (filter, bands, inverse) need their own
%   S.results.
if isfield(S, 'run_name') && ~isempty(S.run_name)
    run_name = S.run_name;
else
    method = 'modkmeans';
    if isfield(S, 'clustering') && ~isempty(S.clustering), method = char(S.clustering); end
    rule = 'median_vote';
    if isfield(S, 'k_selection') && ~isempty(S.k_selection), rule = char(S.k_selection); end
    if isfield(S, 'k_fixed') && ~isempty(S.k_fixed), rule = sprintf('fixed%d', S.k_fixed); end
    fp = hash_text(jsonencode(stage2_settings(S)));
    run_name = sprintf('%s_%d_%d_run_1_%d_restarts_k-%s_%s', method, min(S.NCs), max(S.NCs), S.reps, rule, fp(1:6));
    % another cap or atlas than the study's: in the name, so runs do not overwrite each other
    if ~strcmp(S.montages{1}, 'EEGLAB_61') || ~strcmp(S.atlases{1}, 'Desikan_Killiany_68')
        run_name = [run_name '_' file_tag(S.montages{1}) '_' file_tag(S.atlases{1})];
    end
end
d = fullfile(S.results, 'Connectivity', 'BNS analysis', run_name);
end
