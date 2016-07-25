function [Zwt_p, Zwt_n, TFwt_p, TFwt_n, Zhvdc_p, Zhvdc_n, TFhvdc_p, TFhvdc_n] = convertersImpedanceModel(H,Sr_wt)
%ConverterImpedanceModel calculates impedance models of WT's converters and
%HVDC converter with parameters defined inside this function:
%   in WT: voltage and power at 400MVA and 150kV, L output changed to 1.2
%   like in Aalto data.
%   in HVDC: Y(s) consists of only Cf so tuned filter from Aalto data, L
%   input of converter is named phase reactor in Aalto data, value changed.
%   Control values original from Liu&Sun,2014.
%
%   INPUT: resolution of haronics analysis and rated power of WT converter,
%   depends on case considered.

f1 = 50;
w1 = 2*pi*f1;

%% WT converter
% Vlineline_peak = m/2 * Vdc * sqrt(3)

Vll_to = 150e3; % Vlineline_rms
V1_to = Vll_to*sqrt(2)/sqrt(3); % Vphase_peak
V0_to = 2*V1_to*sqrt(3)/(sqrt(3)*0.75); % Vdc
I1_to = Sr_wt/(sqrt(3)*(sqrt(3)*V1_to));

V0 = 1500; % Vdc
V1 = 563; % Vphase_peak
I1 = 236e3; % Iphase_peak

% L =  1.2 ;% 0.526e-6;
L = 1.2*(V0/V0_to)^2; % convert our 150kV value to the 1500Vdc level

% Hi(s) current control compensator
Kp_i = 0.44e-6;
Ki_i = 0.55e-3;
% Hp(s) PLL compensator
Kp_p = 0.239;
Ki_p = 45;

Hi = @(s) Kp_i+Ki_i/s;
Hp = @(s) (Kp_p+Ki_p/s)/s;
Tpll = @(s) V1*Hp(s)/(2*(1+V1*Hp(s)));

ZpM = zeros(1,length(H));
ZnM = zeros(1,length(H));
for hh=1:length(H)
    h=H(hh);
    s = h*1i*w1;
    
    if (s-1i*w1)==0
        bignum = 10e4;
        ZpM(hh) = (bignum*V0 + (s-1i*w1)*L)/...
            (1-(V1*bignum/(2*(1+V1*bignum)))*(1+bignum*I1*V0/V1));
    else
        ZpM(hh) = (Hi(s-1i*w1)*V0 + (s-1i*w1)*L)/...
            (1-Tpll(s-1i*w1)*(1+Hi(s-1i*w1)*I1*V0/V1));
    end
    ZnM(hh) = (Hi(s+1i*w1)*V0 + (s+1i*w1)*L)/...
        (1-Tpll(s+1i*w1)*(1+Hi(s+1i*w1)*I1*V0/V1));
end
% conversion from 1500dc to our 150kVac
Zwt_p = ZpM * (V1_to/V1)^2;
Zwt_n = ZnM * (V1_to/V1)^2;
% ----------------------
clear ZpM ZnM;

s = tf('s');
TFwt_p = (Hi(s-1i*w1)*V0 + (s-1i*w1)*L)/...
    (1-Tpll(s-1i*w1)*(1+Hi(s-1i*w1)*I1*V0/V1));
TFwt_n = (Hi(s+1i*w1)*V0 + (s+1i*w1)*L)/...
    (1-Tpll(s+1i*w1)*(1+Hi(s+1i*w1)*I1*V0/V1));
clear s
%% HVDC converter

Vdc = 300e3; % HVDC DC link voltage - this V/I levels correspond to 150 kV - no conversion
L = 19.3e-3; % our ph. inductance (PHR at 150kV)
Cf = 5.658e-6; % our Tuned filter capacitance

% values from liu sun 2014
% Lf = 0.43e-6;
% Cf = 17.7e-6;
% Rf = 10;

% Hi(s) current control compensator
Kp_i = 0.075e-3;
Ki_i = 0.094;
Kid = 23e-6; % current decoupling coefficient
% Hp(s) PLL compensator
Kp_v = 11.1e-3;
Ki_v = 8.388;

Hi = @(s) Kp_i+Ki_i/s;
Hv = @(s) Kp_v+Ki_v/s;
Tp = @(s) (Hi(s-1i*w1) + 1i*Kid)*Hv(s-1i*w1)*Vdc;
Tn = @(s) (Hi(s+1i*w1) - 1i*Kid)*Hv(s+1i*w1)*Vdc;
Yf = @(s) s*Cf; % aalto tuned C filter
% Yf = @(s) imp_parallel(s*Cf,(1/(s*Lf)+(1/Rf))); % other tuned RLC filter

ZpM = zeros(1,length(H));
ZnM = zeros(1,length(H));
for hh=1:length(H)
    h=H(hh);
    s = h*1i*w1;
    
    if (s-1i*w1)==0
        bignum = 10e4;
        ZpM(hh) = (bignum*Vdc + s*L)/...
            (1+Yf(s)*(bignum*Vdc+s*L)+((bignum+1i*Kid)*bignum*Vdc));
    else
        ZpM(hh) = (Hi(s-1i*w1)*Vdc + s*L)/...
            (1+Yf(s)*(Hi(s-1i*w1)*Vdc+s*L)+Tp(s));
    end
    ZnM(hh) = (Hi(s+1i*w1)*Vdc + s*L)/...
        (1+Yf(s)*(Hi(s+1i*w1)*Vdc+s*L)+Tn(s));
end
Zhvdc_p = ZpM;
Zhvdc_n = ZnM;
clear ZpM ZnM;

s = tf('s');
TFhvdc_p = (Hi(s-1i*w1)*Vdc + s*L)/...
    (1+Yf(s)*(Hi(s-1i*w1)*Vdc+s*L)+Tp(s));
TFhvdc_n = (Hi(s+1i*w1)*Vdc + s*L)/...
    (1+Yf(s)*(Hi(s+1i*w1)*Vdc+s*L)+Tn(s));
clear s

end

