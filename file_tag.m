function t = file_tag(name)
%FILE_TAG Short name of a montage or atlas for file names.
%   t = file_tag(name): the name itself ('EEGLAB_61', 'Desikan_Killiany_68')
%   or, for a file given in the settings, its name without folder and
%   extension, made safe for file names.
if isfile(name)
    [~, name] = fileparts(name);
end
t = regexprep(char(name), '[^\w-]', '_');
end
