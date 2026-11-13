function [ Alpha,gainWQ ] = sourceModelAssumption( gain, gridLoc, gridOrient,weightExp,weightLimit )

% needed method for wmne computing
% gain: nb_channels x (3*nb_sources) (x, y, z per source) with gridOrient
% nb_sources x 3, or nb_channels x nb_sources (orientation already fixed,
% e.g. normal to the cortex) with gridOrient nb_sources x 1 (+1 or -1).
% gridOrient empty: free orientation, the 3 columns of each source kept
% (S.orientation = 'free'); Alpha is then given per column (3 per source).
% Changes from the var-EEG-FC version (Allouch et al. 2023): the number of
% columns per source is read from the gain (var-EEG-FC assumed 3, so a
% fixed-orientation gain failed); vectorised (var-EEG-FC looped over every
% channel and source); unused Wq removed.

% This code was originally developped by Sahar Yassine and Ahmad Mheich based on the brainstorm toolbox codes
% (Tadel, F., Baillet, S., Mosher, J.C., Pantazis, D., Leahy, R.M., 2011.
% Brainstorm: A user-friendly application for MEG/EEG analysis. Computational Intelligence and Neuroscience 2011. https://doi.org/10.1155/2011/879716)

% contact: saharyassine94@gmail.com
%          mheich.ahmad@gmail.com

numDipoles = size(gridLoc,1);
blklength = size(gain,2) / numDipoles;          % 3 (free) or 1 (fixed)
free = isempty(gridOrient);
if blklength ~= round(blklength) || (~free && size(gridOrient,2) ~= blklength)
    error('sourceModelAssumption:size', 'gain has %d columns for %d sources and %d orientation columns', ...
        size(gain,2), numDipoles, size(gridOrient,2));
end
ColNorms2 = sum(gain .* gain);
SourceNorms2 = sum(reshape(ColNorms2,blklength,[]),1);
Alpha2 = SourceNorms2 .^ weightExp;
AlphaMax2 = max(Alpha2);
Alpha2 = max(Alpha2, AlphaMax2 ./ (weightLimit^2));
Alpha2 = AlphaMax2 ./ Alpha2;
Alpha = sqrt(Alpha2);
if free                                         % every column kept, depth-weighted by its source
    Alpha = repelem(Alpha, blklength);
    gainWQ = gain .* Alpha;
    return
end
orient = gridOrient ./ vecnorm(gridOrient, 2, 2);  % ensure unity
G3 = reshape(gain, size(gain,1), blklength, numDipoles);
gainWQ = reshape(sum(G3 .* reshape(orient', 1, blklength, numDipoles), 2), size(gain,1), numDipoles) .* Alpha;

end
