function S = stage_settings(S, stage)
%STAGE_SETTINGS The settings a stage runs with.
%   S = stage_settings(S, stage) returns S when it is given, otherwise
%   pipeline_settings(). A variable S left in the base workspace is never
%   used (it used to be, silently, so edits to pipeline_settings.m could be
%   ignored for a whole session); a note says how to pass it instead.
%   stage: the name of the calling stage function, for the note.
if nargin < 1 || isempty(S)
    try
        if evalin('base', 'exist(''S'', ''var'')')
            fprintf('Note: the variable S of the workspace is not used; to run with it: %s(S)\n', stage);
        end
    catch
    end
    S = pipeline_settings();
end
if ~isstruct(S) || ~isfield(S, 'root')
    error('stage_settings:S', '%s: S must be a settings struct (pipeline_settings or a copy of it)', stage);
end
end
