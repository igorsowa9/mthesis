clear all; clc; close all

% Example 1.1
%% Base case:
Utr_hi = 70e3; % V trafo high voltage
Utr_lo = 20e3; % V trafo low voltage (substation)
Str = 30e6; % VA trafo rated power 
Ztr_k = 0.17; % per-unit impedance
Str_fault_hi = 1850e6; % MVA trafo falut level at 70kV
Qcap = 5.8e6; % VAr capacitor bank
cable_l = 3; % km length of cables
cable_n = 4; % number of cables
cable_C0 = 0.21e-6; % uF/km cable capacitance
% load does not contribute to the falut level
f = 50; % Hz

Q_cable = 2*pi*f*cable_C0*cable_l*cable_n*Utr_lo^2;
Str_fault = Str/Ztr_k; % analogicaly like Ssc=Srated/Zk%
Stot20_fault = Str_fault_hi*Str_fault/(Str_fault_hi+Str_fault);

n_res = sqrt(Stot20_fault/(Qcap + Q_cable));
f_res = n_res*f;

Ltot = Utr_lo^2/(2*pi*f*Stot20_fault);
Ccap = Qcap/(2*pi*f*Utr_lo^2);
Ctot = cable_C0*cable_l*cable_n + Ccap;

f1_res = 1/(2*pi*sqrt(Ltot*Ctot));

%% A) Impact of cable length:

% resonance frequency
cable_l_varA = 0:0.1:20;
Q_cable_varA = 2*pi*f*cable_C0*cable_l_varA*cable_n*Utr_lo^2;
n_res_varA = sqrt(Stot20_fault./(Qcap + Q_cable_varA));
f_res_varA = n_res_varA*f;

Ctot_varA = cable_C0*cable_l_varA*cable_n + Ccap;
f1_res_varA = 1./(2*pi*sqrt(Ltot*Ctot_varA));

figure(1)
plot(cable_l_varA, f_res_varA, cable_l_varA, f1_res_varA)

% impedance 
% Z = (jhwL/(1-h^2w^2LC))
w = 2*pi*f;
h = 5;
Z_varA = abs(1i*h*w*Ltot./(1-h^2*w^2*Ltot*Ctot_varA));

figure(2) 
hold on
plot(cable_l_varA, Z_varA)
% axis([0 20 0 100])

h = 6;
Z_varA = abs(1i*h*w*Ltot./(1-h^2*w^2*Ltot*Ctot_varA));
plot(cable_l_varA, Z_varA)
h = 7;
Z_varA = abs(1i*h*w*Ltot./(1-h^2*w^2*Ltot*Ctot_varA));
plot(cable_l_varA, Z_varA)
hold off

% different frequencies, same capacitance
h = 1:20;
Z_h = abs(1i*h*w*Ltot./(1-h.^2*w^2*Ltot*Ctot));

%% B) Impact of 


