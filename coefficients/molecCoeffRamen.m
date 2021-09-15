% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the molecular coefficient from the Ramen Lidar

function [beta_mol_lambdaR, alpha_mol_lamdaR] = molecCoeffRamen(backscatterCrosssection,R)
% backscattering cross-section [km^2*sr^-1]
% molecNumbDensity [km^-3]
% beta_mol_lambdaR [km^-1*sr^-1]

molecNumbDensity = 2.1145*10^34 - 2.0022*10^33 * R + 5.4585*10^31 * R^2;
beta_mol_lambdaR = molecNumbDensity * backscatterCrosssection; 

alpha_mol_lamdaR = (7.3219*10^-3)-(4.5204*10^-4)*R;  %[km-1]

end