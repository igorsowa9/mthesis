function [Z] = imp_cable_pi(Zrl,Zrc)
%IMP_CABLE_PI impedance of cable in PI model
%   Arguments are impedance of series RL and impedance of 
%   shunt RC, provided that both shut RC are the same
Za = Zrl*Zrc/(Zrl+Zrc);
Z = Za*Zrc/(Za+Zrc);
end

