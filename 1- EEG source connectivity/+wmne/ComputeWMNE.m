function [ kernel ] = ComputeWMNE( noiseCov,gain,gridLoc,gridOrient,weightExp,weightLimit,SNR )
% Compute the kernel matrix by solving the inverse problem
% inputs: - noiseCov: the noise covariance matrix, nb_channels*nb_channels
%         - gain: the gain Matrix (leadfield matrix ) obtained after
%         computing the head model, nb_channels*(3*nb_regions), nb_regions denotes the number
%         of cortical sources, nb_channels denotes the number of EEG channels;
%         or nb_channels*nb_regions for a fixed orientation (gridOrient = ones(nb_regions,1))
%         - gridLoc: the position (x,y,x) of the sources, nb_regions*3
%         - gridOrient: the orientation of the sources, nb_regions*3 (or nb_regions*1)
%         - weighExp: parameter from brainstorm to compute the wMNE
%         - weightLimit: parameter from brainstorm to compute the wMNE
%         - SNR: signal to noise ratio
%
% Output: - kernel: nb_regions*nb_channels, nb_regions denotes the number
%           of cortical sources, nb_channels denotes the number of EEG channels.
%
% Changes from the var-EEG-FC version (Allouch et al. 2023): a gain with
% more channels than the noise covariance stops the run (var-EEG-FC
% silently kept its first channels); the depth weighting is applied in
% one product (vectorised).

% This code was originally developped by Sahar Yassine and Ahmad Mheich based on the brainstorm toolbox codes
% (Tadel, F., Baillet, S., Mosher, J.C., Pantazis, D., Leahy, R.M., 2011.
% Brainstorm: A user-friendly application for MEG/EEG analysis. Computational Intelligence and Neuroscience 2011. https://doi.org/10.1155/2011/879716)

% contact: saharyassine94@gmail.com
%          mheich.ahmad@gmail.com

%% Compute the inverse problem

if size(noiseCov,2) ~= size(gain,1)
    error('ComputeWMNE:size', 'unmatched dimensions between the gain matrix (%d channels) and the noise covariance (%d)', ...
        size(gain,1), size(noiseCov,2));
end
numberOfChannels = size(noiseCov,1);
montageMatrix = eye(numberOfChannels)-(ones(numberOfChannels,numberOfChannels)./numberOfChannels);
CNoise = (montageMatrix*noiseCov)*montageMatrix';
CNoise = (CNoise+CNoise')./2;
gain = montageMatrix*gain;

[alpha,gainWQ] = wmne.sourceModelAssumption(gain,gridLoc,gridOrient,weightExp,weightLimit);
[~,iW] = wmne.TruncateAndRegularizeCovariance(CNoise);
L = iW*gainWQ;
[UL,SL2] = svd((L*L'));
SL2 = diag(SL2);
SL = sqrt(SL2);
SNR = SNR*SNR;
lambda = SNR/mean(SL.^2);

kernel = lambda * L' * (UL * diag(1./(lambda * SL2 + 1)) * UL');
kernel = kernel * iW;
kernel = alpha(:) .* kernel;

end
