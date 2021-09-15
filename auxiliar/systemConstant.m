% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the system constant

%{
E= Laser energy
c=speed of light
d=telescope lens diamater, pass from inch to m
dsh=dsecondary=dShade

K= W*Km^3

%}

function [K,area_ef] = systemConstant(E,d,dsh,c)

area_ef= pi*( (d/2)^2 - (dsh/2)^2 ); %m^2
K=(E*c/2)*area_ef; % this is in W*Km
K=K*10^-9; %now is in W*Km^3

end

