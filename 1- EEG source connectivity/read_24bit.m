function dat = read_24bit(filename, offset, numwords)
%READ_24BIT Signed 24-bit samples of a BioSemi BDF file, where FieldTrip has no mex file for the platform.
%   dat = read_24bit(filename, offset, numwords): numwords samples from
%   byte offset, as a 1 x numwords double row, as FieldTrip's
%   fileio/private/read_24bit mex file returns them. FieldTrip's BDF reader
%   (read_biosemi_bdf) calls that mex file, which the validated FieldTrip
%   (commit 0beabe6, April 2023) has in fileio/private for every platform
%   but Apple silicon (mexmaca64): there, every BDF file would fail with
%   "Undefined function 'read_24bit'". MATLAB looks
%   for a private function first, so this file is used only where the mex
%   file is missing; it reads the samples as FieldTrip's own code does for
%   files over 2 GB (fread 'bit24').
fid = fopen(filename, 'r', 'ieee-le');
if fid < 0, error('read_24bit:open', 'Cannot open %s', filename); end
c = onCleanup(@() fclose(fid));
if fseek(fid, offset, 'bof') ~= 0, error('read_24bit:seek', 'Cannot seek to byte %d of %s', offset, filename); end
[dat, n] = fread(fid, numwords, 'bit24=>double');
if n < numwords, error('read_24bit:short', '%s: %d of %d samples read', filename, n, numwords); end
dat = dat(:)';
end
