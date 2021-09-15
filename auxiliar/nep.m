% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the NEP from the Lidar

%{
Laser is off! Only thermal and dark current noises are present

NEP= photodiode NEP, only the photodiode, increase lidar power until it matches the dark shot signal
NEP_s= All systems NEP, all system as described on slide 3 on link ppt

%}

function [NEP,NEP_s,quantum_eficiency] = nep(sigma_photodiode_sh_d,sigma_sh_d,sigma_th, Ri,Rv_net, M,q,h,c,em_wL)

NEP = sigma_photodiode_sh_d/Ri; % [A/sqrt(HZ)] / [A/W] = [W/sqrt(Hz)]
NEP_s = sqrt((sigma_sh_d^2) + (sigma_th^2))/(Rv_net);

quantum_eficiency = (Ri*h*c) / (q*M*em_wL);
quantum_eficiency=quantum_eficiency*100;%pass the efficiency to percentage

% We pass the values from Watts to Femto Watts
NEP=NEP*10^15;
NEP_s=NEP_s*10^15;

end