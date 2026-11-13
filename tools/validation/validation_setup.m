function [root, n_workers] = validation_setup()
%VALIDATION_SETUP The repository's statistics functions on the path, for the scripts of tools/validation.
%   [root, n_workers] = validation_setup()
%   root       the repository (two folders above this one); it and the
%              fun folder of its analysis stage ('4- Analysis', or '3- Analysis'
%              without the states) are added at the top of the MATLAB path,
%              so the scripts of this folder always test the code they
%              sit in
%   n_workers  the workers of a parallel pool that is already open; 0 when
%              none is open or without the Parallel Computing Toolbox. The
%              scripts that use parfor give it as parfor's largest number
%              of workers: they never open a pool themselves, and run
%              serially (in the same order, with the same results) without
%              one. To use several cores, open a pool first (parpool).
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root);
addpath(fullfile(stage_folder(root, 'Analysis'), 'fun'));
n_workers = 0;
if exist('gcp', 'file') == 2 && license('test', 'Distrib_Computing_Toolbox')
    try
        pool = gcp('nocreate');
        if ~isempty(pool), n_workers = pool.NumWorkers; end
    catch
        n_workers = 0;
    end
end
end
