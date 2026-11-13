function wpli = calc_wpli(eeg_data)
%WPLI calculate weighted PLI (reference implementation, kept for the tests;
%   plv_wpli_sliding_window computes the same quantity much faster)
%   Input:
%       eeg_data: time-series to calculate wpli
%      
%   Output:
%       wpli: wPLI matrix
%% Setting up variables
    data = eeg_data';
    number_channel = size(data,2); 
    a_sig = hilbert(data);
    wpli = zeros(number_channel,number_channel);   % zero diagonal

    for channel_i=1:number_channel-1
        for channel_j=channel_i+1:number_channel
            c_sig=a_sig(:,channel_i).*conj(a_sig(:,channel_j));

            numerator=abs(mean(imag(c_sig))); % average of imaginary
            denominator=mean(abs(imag(c_sig))); % average of abs of imaginary

            wpli(channel_i,channel_j) = numerator/denominator;
            wpli(channel_j,channel_i) = wpli(channel_i,channel_j);
        end
    end 
end