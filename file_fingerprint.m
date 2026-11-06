function h = file_fingerprint(files, max_bytes)
%FILE_FINGERPRINT Fingerprint of files: their names and contents.
%   h = file_fingerprint(files) with files a cell of paths: a hash of each
%   file's name and content (MD5 through Java), so a copy of the same file on
%   another computer or with another date gives the same fingerprint. Files
%   larger than max_bytes (default 500 MB) are represented by name and size
%   only. Missing files count as missing. Contents are hashed once per MATLAB
%   session (cached by path, size and modification time; files modified in
%   the last 2 s are always re-read).
if nargin < 2, max_bytes = 500e6; end
persistent cache
if isempty(cache), cache = containers.Map(); end
files = cellstr(files);
parts = cell(1, numel(files));
for i = 1:numel(files)
    [~, name, ext] = fileparts(files{i});
    d = dir(files{i});
    if isempty(d) || d(1).isdir
        parts{i} = [name ext ':missing'];
        continue
    end
    modified = d(1).datenum * 86400e3;                    % ms; Java's clock is finer than dir's
    try, modified = double(java.io.File(files{i}).lastModified()); catch, end
    key = sprintf('%s|%d|%.0f', files{i}, d(1).bytes, modified);
    recent = (now - d(1).datenum) * 86400 < 2;            % just written: may change again within the clock's step
    if isKey(cache, key) && ~recent
        parts{i} = cache(key);
        continue
    end
    if d(1).bytes > max_bytes
        v = sprintf('%s%s:%d', name, ext, d(1).bytes);
    else
        v = sprintf('%s%s:%s', name, ext, content_md5(files{i}));
    end
    if ~recent, cache(key) = v; end
    parts{i} = v;
end
h = hash_text(strjoin(parts, '|'));
end

function h = content_md5(file)
fid = fopen(file, 'r');
if fid < 0, h = 'unreadable'; return; end
c = onCleanup(@() fclose(fid));
try
    md = java.security.MessageDigest.getInstance('MD5');
    while true
        b = fread(fid, 2^24, '*uint8');
        if isempty(b), break; end
        md.update(b);
    end
    h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, []));
catch                                     % no Java: the text hash of the bytes
    frewind(fid);
    h = hash_text(char(fread(fid, inf, '*uint8')'));
end
end
