clear all; close all; clc

%% System data:
% Subst. transformer (T2)
Tr2_HV = 115e3;
Tr2_LV = 34.5e3;
Tr2_S = 250e6;
Tr2_Zp = 0.18;
Tr2_XRratio = 12;
Tr2_n = 2;
% Turbine transformer (T1)
Tr1_HV = 34.5e3;
Tr1_LV = 690;
Tr1_S = 2.5e6;
Tr1_Zp = 0.05;
Tr1_XRratio = 5;
Tr1_n = 100;
% Induction machine
IM_V = 690;
IM_S = 2e6;
IM_Zp = 0.1712;
IM_XRratio = 8.5;
IM_n = 100;
% Cables
Cab_V = 34.5e3;
Cab_R0 = 0.13; %ohm/km
Cab_C0 = 0.25e-6; %F/km
Cab_XRratio = 18;
% Capacitor bank
Cap_V = 34.5e3;
Cap_S = 12e6; % VAr
Cap_n = 6;
% OHL
OHL_V = 115e3;
OHL_R0 = 0.36; % ohm/km
OHL_length = 0;
% Grid
Grid_V = 115e3;
Grid_Sf = 3500e6;
Grid_XRratio = 18;

% Their impedances (already to 34.5kV):
Tr2_Z = 0.07117 + 1i*0.8540;
Tr1_Z = 4.669 + 1i*23.34; % for all WTs (100)
IM_Z = 10.92 + 1i*92.83; % for all WTs (100)
Cab_Z = -1i*277.2;
OHL_Z = 0.05936 + 1i*0.5639; % both lines (one line /2 - parallel conn!)
Grid_Z = 0.01886 + 1i*0.3396;
Cap6x_Z = -1i*16.53;

% Their method:
H = 0:0.001:10;

% Series resonance - when minimal abs impedance seen from HV substation
Z_HVabsM = zeros(length(H),1);
U_MVdistM = zeros(length(H),1); % amplification of V distortion at MV subst.

for hh=1:length(H)
    h=H(hh);
    
    Z_OHLTr2 = h*(OHL_Z+imp_parallel(Tr2_Z,Tr2_Z)); % Connection like in paper
    Z_WF = h*(Tr1_Z+IM_Z)/100; % divided by 100 to obtain 1 branch
    Z_C = imp_parallel(Cab_Z/h,Cap6x_Z/h);
    Z_HV = imp_parallel(h*Grid_Z,imp_parallel(Z_C,Z_WF)+Z_OHLTr2);
    
    Z_HVabs = abs(Z_HV);
    Z_HVabsM(hh) = Z_HVabs;
    
    Umv = imp_parallel(Z_C,Z_WF);
    Uhv = imp_parallel(Z_C,Z_WF)+Z_OHLTr2;
    U_MVdistM(hh)=abs(Umv/Uhv);
end
figure(1)
plot(H,Z_HVabsM, 'LineWidth', 1)
axis([0 10 0 3]);
xlabel('harmonics');
ylabel('Impedance [ohms]');
title('Impedance seen from HV substation');

figure(2)
plot(H,U_MVdistM, 'LineWidth', 1)
axis([0 10 0 5]);
xlabel('harmonics');
ylabel('U_MV / U_HV');
title('Voltage distortion amplification');









