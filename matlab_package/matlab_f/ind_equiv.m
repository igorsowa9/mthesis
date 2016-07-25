function [Lout] = ind_equiv(Lin,f,Vin,Vout)
%IND_EQUIV calculates inductance to another voltage level
%   inductance Lfrom is converted from voltage Vfrom to equivalent
%   value Lto, at voltage Vto. For frequency f.
w = 2*pi*f;
Xl_in = w*Lin;
n = (Vin/Vout);
Xl_out = Xl_in * n^2;
Lout = Xl_out/w;
end

