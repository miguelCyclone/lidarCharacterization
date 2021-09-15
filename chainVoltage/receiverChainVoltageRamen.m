% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This script calculates the chain voltage from the Ramen lidar

function [Rv_ramen,Rv_net_ramen, Rio_ramen] = receiverChainVoltageRamen(R_i_ramen, M_ramen, G_T_ramen, L_ramen)

Rio_ramen= R_i_ramen / M_ramen; % A/W Rio = photodiode intrinsic current responsivity, Ri = photodiode current responsivity 
Rv_ramen= R_i_ramen * G_T_ramen;%V/W voltage responsivity
Rv_net_ramen = Rv_ramen * L_ramen;%V/W net voltage responsivity

end

