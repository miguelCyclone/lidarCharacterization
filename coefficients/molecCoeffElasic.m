% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the molecular coefficient from the Elastic lidar

% the US model combines the theoriteical result from Lord Rayleigh for
% moleculas species follow the lambda to the -4, thats the -4 in the US model

function [alpha_mol, beta_mol] = molecCoeffElasic(S_r,R)
alpha_mol = (1.2593*10^-2)-(7.7599*10^-4)*R; % Rayleigh (molec) extinction [km^-1]
beta_mol = alpha_mol / S_r;
end