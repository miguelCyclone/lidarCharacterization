% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the background power

%{
Lradiation = background radiance component, it can be moon or sun
Lradiation units = W / (cm^2 * nm * str)
%}
function pBack = Pback(Lradiation,area_ef,fl,Dapd,if_bw)

area_ef =area_ef*10^4; %from m^2 to cm^2

rd=Dapd/2; %check! what is rd?
FOV=rd/fl;
omega=pi*sin(FOV)^2; % in rad! and actually it is: area/radius^2, but the area is pi(R*sin(FOV))^2

pBack=Lradiation*area_ef*if_bw*omega;

end

