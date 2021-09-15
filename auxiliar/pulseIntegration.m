% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script enables a puls eintegration to the laser to achieve a higher range and accuracy (SNR)

function [t_obs, maxRangeObsTime]=pulseIntegration(SNR, SNRN, PRF, maxObserbationTime, log)

aux = 0;
maxRangeObsTime = 1;
for i = 1:length(SNR)
    N=(10^(SNRN/log))/(10^(SNR(i)/log));%this is in db, log = 20 if dBV, or 10 if it is dBW
    N=N^2;
    N=ceil(N);
    t_obs(i)=N*(1/PRF);
    if(t_obs(i) > maxObserbationTime)
        t_obs(i) = maxObserbationTime;
        if(aux == 0)
            maxRangeObsTime = i;   
            aux = 1;
        end
    end
        
end
end