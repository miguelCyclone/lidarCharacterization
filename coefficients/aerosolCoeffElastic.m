% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This function calculates the aerosol coefficient for the Elastic Lidar

function [alpha_aer,beta_aer] = aerosolCoeffElastic(Vm,Sm)
alpha_aer = 3.912/Vm;
beta_aer = alpha_aer/Sm;
end