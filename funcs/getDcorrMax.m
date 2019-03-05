% function [ind,snr] = getDcorrMax(d)

function [ind,snr] = getDcorrMax(d)

[snr,ind] = max(d);
t = d;
% arbitrary peak significance parameter imposed by numerical noise
% this becomes relevant especially when working with post-processed data
dt = 0.001;

while ind == length(t)
    t(end) = [];
    if isempty(t)
        snr = 0;
        ind = 1;
    else
        [snr,ind] = max(t);
        % check if the peak is significantly larger than the former minimum
        if t(ind) - min(d(ind:end)) > dt 
            break
        else
            t(ind) = min(d(ind:end));
            ind = length(t);
        end
    end
     
end