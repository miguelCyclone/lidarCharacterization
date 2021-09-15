% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the receviving chain voltage from the Elastic lidar

%{
P[W]+L[.]*Ri[A/W]*Gt[V/A]=V

Voltage responsivity = Rv[V/W] = Ri[A/W]*Gt[V/A] ; Ri = Rio[A/W]*M[.]
GT=Gt*Gacn

%}
function [Rv,Rv_net, Ri] = receiverChainVoltage(R_io,M,G_T,L)
Ri=R_io*M; %A/W photodiode current responsivity
Rv=Ri*G_T;%V/W voltage responsivity
Rv_net = Rv*L;%V/W net voltage responsivity
end

