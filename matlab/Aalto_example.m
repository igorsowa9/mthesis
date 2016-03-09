clear all; close all; clc

% system data:
LCL_L1 = 1.2; % H
LCL_R1 = 0.0; % ohm %%%%%%
LCL_L2 = 0.641; % H
LCL_R2 = 0.0; % ohm %%%%%%
LCL_C = 1.491e-7; % F
LCL_Rc = 0; % ohm %%%%%%
Tr1_L = 51.568e-3; % H

Cable33_L = 18.181e-3; % H
Cable33_R = 0.372; % ohm
Cable33_C1 = 5.709e-8; % F
Cable33_C2 = Cable33_C1; % F
Tr2_L = 38.676e-3; % H

Cable150_L = 1e-3; % H
Cable150_R = 0.056; % ohm
Cable150_C1 = 0.52e-6; % F
Cable150_C2 = Cable150_C1; % F
Tr3_L = 19.338e-3; % H
Tuned_C = 5.658e-6; % F
PhReact_L = 19.3e-3; % H

% -------
f = 50;
w = 2*pi*f;
H = 0:0.001:30;

% convertion to equivalent circuit at V=150kV
Vout = 150e3;
LCL_L1 = ind_equiv(LCL_L1,f,8e3,Vout);
LCL_R1 = res_equiv(LCL_R1,8e3,Vout);
LCL_L2 = ind_equiv(LCL_L2,f,8e3,Vout);
LCL_R2 = res_equiv(LCL_R2,8e3,Vout);
LCL_C  = cap_equiv(LCL_C,f,8e3,Vout);
LCL_Rc = res_equiv(LCL_Rc,8e3,Vout);
Tr1_L  = ind_equiv(Tr1_L,f,8e3,Vout);

Cable33_L  = ind_equiv(Cable33_L,f,33e3,Vout);
Cable33_R  = res_equiv(Cable33_R,33e3,Vout);
Cable33_C1 = cap_equiv(Cable33_C1,f,33e3,Vout);
Cable33_C2 = cap_equiv(Cable33_C2,f,33e3,Vout);
Tr2_L = ind_equiv(Tr2_L,f,33e3,Vout);

% -------
ymax = 5e3;
view = [0,H(length(H)),0,ymax];
% view = [9 11 0 100e4];

% Impedance - 1 leg (seen from WT)
ZabsM1 = zeros(length(H),1);
for hh=1:length(H)
    h=H(hh);
    s = 1i*h*w;
    Z1 = s*Tr3_L + imp_parallel(s*PhReact_L, 1/(s*Tuned_C));
    Z2 = imp_cable_pi(Cable150_R+s*Cable150_L, 1/(s*Cable150_C1));
    Z3 = s*Tr2_L + imp_cable_pi(Cable33_R+s*Cable33_L, 1/(s*Cable33_C1));
    Z4 = imp_parallel(LCL_R2 + s*LCL_L2, 1/(s*LCL_C)) + LCL_R1 + s*LCL_L1;

    Z = Z1+Z2+Z3+Z4;
    Zabs = abs(Z);
    ZabsM1(hh) = Zabs;
end
clear hh h Z1 Z2 Z3 Z4 Z Zabs

figure(1)
plot(H,ZabsM1, 'LineWidth', 1)
axis(view)

% Impedance - 2 legs (seen from WT)
ZabsM2 = zeros(length(H),1);
for hh=1:length(H)
    h=H(hh);
    s = 1i*h*w;
    
    Z1 = s*Tr3_L + imp_parallel(s*PhReact_L, 1/(s*Tuned_C));
    Z2 = imp_cable_pi(Cable150_R+s*Cable150_L, 1/(s*Cable150_C1));
    Z3 = s*Tr2_L + imp_cable_pi(Cable33_R+s*Cable33_L, 1/(s*Cable33_C1));
    Z4 = imp_parallel(LCL_R2 + s*LCL_L2, 1/(s*LCL_C)) + LCL_R1 + s*LCL_L1;

%     Z = imp_parallel(Z3+Z4,Z3+Z4)+Z2+Z1;
    Z = imp_parallel(Z1+Z2,Z3+Z4)+Z3+Z4;
    Zabs = abs(Z);
    ZabsM2(hh) = Zabs;
end
clear hh h Z1 Z2 Z3 Z4 Z Zabs

figure(2)
plot(H,ZabsM2, 'LineWidth', 1)
axis(view)

% Impedance - 4 legs (seen from WT)
ZabsM3 = zeros(length(H),1);
for hh=1:length(H)
    h=H(hh);
    s = 1i*h*w;
    
    Z1 = s*Tr3_L + imp_parallel(s*PhReact_L, 1/(s*Tuned_C));
    Z2 = imp_cable_pi(Cable150_R+s*Cable150_L, 1/(s*Cable150_C1));
    Z3 = s*Tr2_L + imp_cable_pi(Cable33_R+s*Cable33_L, 1/(s*Cable33_C1));
    Z4 = imp_parallel(LCL_R2 + s*LCL_L2, 1/(s*LCL_C)) + LCL_R1 + s*LCL_L1;

    Z5 = imp_parallel(Z3+Z4,Z3+Z4)+Z2;
    Z6 = imp_parallel(Z1,Z5)+Z2;
    Z7 = imp_parallel(Z6,Z3+Z4);
    Z = Z7 + Z3+Z4;

    Zabs = abs(Z);
    ZabsM3(hh) = Zabs;
end
clear hh h Z1 Z2 Z3 Z4 Z Zabs

figure(3)
plot(H,ZabsM3, 'LineWidth', 1)
axis(view)
