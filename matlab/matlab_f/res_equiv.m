function [ Rout ] = res_equiv( Rin, Vin, Vout )
%RES_EQUIV Convert resistance to another voltage lvl
%   
Rout = Rin*(Vin/Vout)^2;
end

