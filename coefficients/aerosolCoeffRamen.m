% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the aerosocl coefficient from the Ramen Lidar

% No backscqatter from the aerosol for Rammen
% Aerosol species follow the lambda to the -1 to -1.8 from experimental results

function alpha_aer_lambdaR = aerosolCoeffRamen(alpha_aer,em_wL,lambda_r)
k = 1.8; %
alpha_aer_lambdaR = alpha_aer * (lambda_r/em_wL)^-k; % []
end