function g = default_group_files(groups)
%DEFAULT_GROUP_FILES The suffixes of the per-group files of stage 2 when S.group_files is {}.
%   g = default_group_files(groups): {'ctrl', 'psych'} for two groups, in the
%   order of S.groups (the original file names, kept so that the run folders
%   of stage 2, whose names follow the settings, stay the same), whatever
%   the groups are called; otherwise the group names, lower case, made valid
%   names (letters, digits, _: they are part of variable names). They name
%   files only: the tables and figures of stage 3 name the groups by
%   S.groups and S.comparisons. Used by check_settings and stage2_settings.
groups = cellstr(groups);
if numel(groups) == 2
    g = {'ctrl', 'psych'};
else
    g = matlab.lang.makeValidName(lower(groups));
end
end
