function plv = compute_plv(data)

nchannels = size(data,1);
plv = zeros(nchannels,nchannels);
    for channelCount = 1:nchannels-1
        channelData = squeeze(data(channelCount,:));
        for compareChannelCount = channelCount+1:nchannels
            compareChannelData = squeeze(data(compareChannelCount,:));
                diff=channelData(:, :) - compareChannelData(:, :);
                diff=diff';
               plv(channelCount,compareChannelCount) = abs(sum(exp(1i*diff)))/length(diff);  
        end
    end
plv = plv + plv';    
end

