% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script runs the Q9 that aims to obtain a pulse integration vs observation time range

function [SNR]=question9(gamma, PRF, tObs_simul, rvNet, P)
SNR = zeros(length(tObs_simul),0);
gamma = 10^(gamma/20); % from dbV to linear
for i = 1:length(tObs_simul)
    N= tObs_simul(i) / (1/PRF);
    N=ceil(N);
    SNRN(i) = sqrt(N); %SNR with "pulse integration"
    
end

end