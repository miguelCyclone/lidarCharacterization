% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the maximum range that the Lidar can get

function out=maxRange(SNR, R, gamma)
out = zeros(length(SNR),0);
for i = 1:length(SNR)
    if SNR(i) <= gamma
        out(i) = R(i);
    else
        out(i) = R(length(R));
    end    
end
end