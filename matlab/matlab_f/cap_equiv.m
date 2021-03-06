function [ Cout ] = cap_equiv(Cin,f,Vin,Vout)
%CAP_EQUIV Scalculates capacitance to another voltage level
%   capacitance Cfrom is converted from voltage Vfrom to equivalent
%   value Cto, at voltage Vto. For frequency f.
w = 2*pi*f;
Xc_in = 1/(w*Cin);
n = (Vin/Vout);
Xc_out = Xc_in * n^2;
Cout = 1/(Xc_out*w);
end

