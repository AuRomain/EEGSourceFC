function ftElec = get_ftElec(montage)
%GET_FTELEC Electrode structure (FieldTrip format) of a montage.
%   ftElec = get_ftElec(montage) with montage
%     'EEGLAB_61'      the 61 channels of the FES/HVS recordings, corrected
%                      positions (E1: the original channel_EEGLab_61_bs.mat
%                      had been rotated by about 90 degrees)
%     'ASA_10-05_343'  Brainstorm's ICBM152 ASA 10-05 cap, 343 positions
%                      (including T3/T4/T5/T6, M1/M2, A1/A2): any recording
%                      with standard 10-20 / 10-10 / 10-05 labels
%     a file name      Brainstorm channel_*.mat, or CSV with columns
%                      label,x,y,z (metres, ICBM152 coordinates)
%   The original also listed the Biosemi, Elekta and 256-channel montages of
%   other datasets, with files not part of this study; a cap of another
%   dataset is now given as a file (or matched in the ASA cap by label), and
%   a subset of channels (the original Biosemi 32 and 19) by
%   S.channel_subset (get_eeg_montage).
in = fullfile(fileparts(mfilename('fullpath')), 'inputs');
switch montage
    case 'EEGLAB_61'
        ftElec = BS_to_ft_channels(fullfile(in, 'channel_EEGLab_61_corrected_bs.mat'), 1:61);
    case 'ASA_10-05_343'
        ftElec = BS_to_ft_channels(fullfile(in, 'channel_ASA_10-05_343_bs.mat'), []);
    otherwise
        if isfile(montage)
            ftElec = BS_to_ft_channels(montage, []);
        else
            error('get_ftElec:montage', ['Unknown montage %s: use ''EEGLAB_61'', ''ASA_10-05_343'' ' ...
                'or a channel file'], montage);
        end
end
end
