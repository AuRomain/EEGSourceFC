function file = write_eeg(file, X, fs, labels, fmt)
%WRITE_EEG Write continuous EEG (nChan x nSamples, microvolts) to a file.
%   sim.write_eeg(file, X, fs, labels, fmt) with fmt 'set' (EEGLAB, data in
%   the .set), 'edf' (16 bit), 'bdf' (24 bit), 'vhdr' (BrainVision, float32,
%   multiplexed) or 'mat' (variables data, srate, labels). The extension
%   of FILE is replaced by the format's own. Used for simulated datasets.
%   Stage 1 reads 'set', 'edf', 'bdf' and 'vhdr' (FieldTrip).
[folder, base] = fileparts(file);
if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
labels = cellstr(labels);
switch fmt
    case 'set'
        file = fullfile(folder, [base '.set']);
        EEG = struct('setname', base, 'filename', [base '.set'], 'filepath', folder, ...
            'subject', '', 'group', '', 'condition', '', 'session', [], 'comments', 'simulated', ...
            'nbchan', size(X, 1), 'trials', 1, 'pnts', size(X, 2), 'srate', fs, 'xmin', 0, ...
            'xmax', (size(X, 2) - 1) / fs, 'times', [], 'data', single(X), 'icaact', [], ...
            'icawinv', [], 'icasphere', [], 'icaweights', [], 'icachansind', [], ...
            'chanlocs', struct('labels', labels), 'urchanlocs', [], 'chaninfo', struct(), ...
            'ref', 'average', 'event', struct('type', {}, 'latency', {}), 'urevent', [], ...
            'epoch', [], 'etc', struct());
        save(file, 'EEG', '-v7');
    case {'edf', 'bdf'}
        file = fullfile(folder, [base '.' fmt]);
        write_edf(file, X, fs, labels, strcmp(fmt, 'bdf'));
    case 'vhdr'
        file = fullfile(folder, [base '.vhdr']);
        fid = fopen(fullfile(folder, [base '.eeg']), 'w', 'ieee-le');
        fwrite(fid, single(X), 'float32');
        fclose(fid);
        fid = fopen(file, 'w');
        fprintf(fid, 'Brain Vision Data Exchange Header File Version 1.0\n\n[Common Infos]\n');
        fprintf(fid, 'Codepage=UTF-8\nDataFile=%s.eeg\nMarkerFile=%s.vmrk\nDataFormat=BINARY\n', base, base);
        fprintf(fid, 'DataOrientation=MULTIPLEXED\nNumberOfChannels=%d\nSamplingInterval=%.10g\n\n', ...
            size(X, 1), 1e6 / fs);
        fprintf(fid, '[Binary Infos]\nBinaryFormat=IEEE_FLOAT_32\n\n[Channel Infos]\n');
        for k = 1:numel(labels)
            fprintf(fid, 'Ch%d=%s,,1,%s\n', k, labels{k}, char([181 86]));
        end
        fclose(fid);
        fid = fopen(fullfile(folder, [base '.vmrk']), 'w');      % marker file (FieldTrip needs one)
        fprintf(fid, 'Brain Vision Data Exchange Marker File, Version 1.0\n\n[Common Infos]\n');
        fprintf(fid, 'Codepage=UTF-8\nDataFile=%s.eeg\n\n[Marker Infos]\nMk1=New Segment,,1,1,0\n', base);
        fclose(fid);
    case 'mat'
        file = fullfile(folder, [base '.mat']);
        data = X; srate = fs; %#ok<NASGU>
        save(file, 'data', 'srate', 'labels');
    otherwise
        error('write_eeg:format', 'Unknown format %s', fmt);
end
end

function write_edf(file, X, fs, labels, is_bdf)
% One-second records; physical range = data range per channel.
[nch, n] = size(X);
rec = round(fs);
if abs(rec - fs) > 1e-9, error('write_edf:rate', 'EDF writer needs an integer sampling rate'); end
nrec = floor(n / rec);
X = X(:, 1:nrec * rec);
if is_bdf, dmax = 2^23 - 1; else, dmax = 32767; end
dmin = -dmax - 1;
[pmin_s, pmin] = arrayfun(@(v) fit8(v, -1), min(X, [], 2), 'UniformOutput', false);
[pmax_s, pmax] = arrayfun(@(v) fit8(v, 1), max(X, [], 2), 'UniformOutput', false);
pmin = [pmin{:}]';  pmax = [pmax{:}]';
D = round((X - pmin) ./ (pmax - pmin) * (dmax - dmin) + dmin);
fid = fopen(file, 'w', 'ieee-le');
pad = @(s, w) sprintf(['%-' num2str(w) 's'], s(1:min(end, w)));
num = @(v, w) pad(sprintf('%.8g', v), w);
if is_bdf
    fwrite(fid, 255, 'uint8');  fprintf(fid, '%s', pad('BIOSEMI', 7));
else
    fprintf(fid, '%s', pad('0', 8));
end
fprintf(fid, '%s%s%s%s', pad('X X X X', 80), pad('Startdate X X X X', 80), '01.01.26', '00.00.00');
fprintf(fid, '%s%s', num(256 * (nch + 1), 8), pad(ifelse(is_bdf, '24BIT', ''), 44));
fprintf(fid, '%s%s%s', num(nrec, 8), num(1, 8), num(nch, 4));
fprintf(fid, '%s', strjoin(cellfun(@(l) pad(l, 16), labels, 'UniformOutput', false), ''));
fprintf(fid, '%s', repmat(pad('', 80), 1, nch));
fprintf(fid, '%s', repmat(pad('uV', 8), 1, nch));
fprintf(fid, '%s', strjoin(cellfun(@(v) pad(v, 8), pmin_s, 'UniformOutput', false), ''));
fprintf(fid, '%s', strjoin(cellfun(@(v) pad(v, 8), pmax_s, 'UniformOutput', false), ''));
fprintf(fid, '%s', repmat(num(dmin, 8), 1, nch));
fprintf(fid, '%s', repmat(num(dmax, 8), 1, nch));
fprintf(fid, '%s', repmat(pad('', 80), 1, nch));
fprintf(fid, '%s', repmat(num(rec, 8), 1, nch));
fprintf(fid, '%s', repmat(pad('', 32), 1, nch));
R = reshape(permute(reshape(D, nch, rec, nrec), [2, 1, 3]), [], 1);
if is_bdf
    fwrite(fid, R, 'bit24');
else
    fwrite(fid, R, 'int16');
end
fclose(fid);
end

function [s, v] = fit8(v, direction)
% Physical limit written in 8 characters, rounded outwards so the data fit.
for d = 6:-1:0
    q = 10 ^ d;
    if direction < 0, w = floor(v * q) / q; else, w = ceil(v * q) / q; end
    s = sprintf('%.*f', d, w);
    if numel(s) <= 8, v = str2double(s); return; end
end
error('write_edf:range', 'Value %g does not fit an EDF header field', v);
end

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
