clear all; close all; clc

%% system data:
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

%% convertion to equivalent circuit at V=150kV
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

Cable150_L = ind_equiv(Cable150_L,f,150e3,Vout);
Cable150_R = res_equiv(Cable150_R,150e3,Vout);
Cable150_C1 = cap_equiv(Cable150_C1,f,150e3,Vout);
Cable150_C2 = cap_equiv(Cable150_C2,f,150e3,Vout);
Tr3_L = ind_equiv(Tr3_L,f,150e3,Vout);
Tuned_C = cap_equiv(Tuned_C,f,150e3,Vout);
PhReact_L = ind_equiv(PhReact_L,f,150e3,Vout);

% -------
ymax = 200;
view = [0,H(length(H)),0,ymax];
% view = [9 14 0 100]; % zoom
xlab = 'harmonics order';
ylab = 'impedance (ohms)';

%% Impedance - 1 leg (seen from WT)
ZabsM1 = zeros(length(H),1);
for hh=1:length(H)
    h=H(hh);
    s = 1i*h*w;
    
%     Z1 = s*Tr3_L + imp_parallel(s*PhReact_L, 1/(s*Tuned_C));
%     Z2 = imp_cable_pi(Cable150_R+s*Cable150_L, 1/(s*Cable150_C1));
%     Z3 = s*Tr2_L + imp_cable_pi(Cable33_R+s*Cable33_L, 1/(s*Cable33_C1));
%     Z4 = imp_parallel(LCL_R2 + s*LCL_L2, 1/(s*LCL_C));% + LCL_R1 +
%     s*LCL_L1;
%     Z = Z1+Z2+Z3+Z4;

    Z1 = imp_parallel(s*PhReact_L, 1/(s*Tuned_C)) + s*Tr3_L;
    Z2 = imp_parallel(Z1, 1/(s*Cable150_C2)) + Cable150_R+s*Cable150_L;
    Z3 = imp_parallel(Z2, 1/(s*Cable150_C1)) + s*Tr2_L;
    Z4 = imp_parallel(Z3, 1/(s*Cable33_C2)) + Cable33_R+s*Cable33_L;
    Z5 = imp_parallel(Z4, 1/(s*Cable33_C1)) + s*Tr1_L + LCL_R2+s*LCL_L2;
    Z6 = imp_parallel(Z5, 1/(s*LCL_C));
    Z7 = LCL_R1+s*LCL_L1;
    Z = imp_parallel(Z6, Z7);
    
    Zabs = abs(Z);
    ZabsM1(hh) = Zabs;
end
clear hh h Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z Zabs

figure(1)
plot(H,ZabsM1, 'LineWidth', 1)
% view1 = [0,H(length(H)),0,300];
axis(view)
title('1. case: 1 WT');
xlabel(xlab);
ylabel(ylab);

%% Impedance - 2 legs (seen from WT)
ZabsM2 = zeros(length(H),1);
for hh=1:length(H)
    h=H(hh);
    s = 1i*h*w;
    
%     Z1 = s*Tr3_L + imp_parallel(s*PhReact_L, 1/(s*Tuned_C));
%     Z2 = imp_cable_pi(Cable150_R+s*Cable150_L, 1/(s*Cable150_C1));
%     Z3 = s*Tr2_L + imp_cable_pi(Cable33_R+s*Cable33_L, 1/(s*Cable33_C1));
%     Z4 = imp_parallel(LCL_R2 + s*LCL_L2, 1/(s*LCL_C));% + LCL_R1 +
%     s*LCL_L1;
%     Z = imp_parallel(Z1+Z2,Z3+Z4)+Z3+Z4;

    Z1 = imp_parallel(s*PhReact_L, 1/(s*Tuned_C)) + s*Tr3_L;
    Z2 = imp_parallel(Z1, 1/(s*Cable150_C2)) + Cable150_R+s*Cable150_L;
    Z3 = imp_parallel(Z2, 1/(s*Cable150_C1));
    
    Z1B = imp_parallel(LCL_R1+s*LCL_L1,1/(s*LCL_C))+LCL_R2+s*LCL_L2+s*Tr1_L;
    Z2B = imp_parallel(Z1B, 1/(s*Cable33_C1)) + Cable33_R+s*Cable33_L;
    Z3B = imp_parallel(Z2B, 1/(s*Cable33_C2)) + s*Tr2_L;
    
    Z4 = imp_parallel(Z3B, Z3) + s*Tr2_L;
    Z5 = imp_parallel(Z4, 1/(s*Cable33_C2)) + Cable33_R+s*Cable33_L;
    Z6 = imp_parallel(Z5, 1/(s*Cable33_C1)) + s*Tr1_L + s*LCL_L2+LCL_R2;
    Z7 = imp_parallel(Z6, 1/(s*LCL_C));
    Z8 = LCL_R1+s*LCL_L1;
    Z = imp_parallel(Z7, Z8);
    
    Zabs = abs(Z);
    ZabsM2(hh) = Zabs;
end
clear hh h Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z1B Z2B Z3B Z Zabs

figure(2)
plot(H,ZabsM2, 'LineWidth', 1)
axis(view)
title('2. case: 2 WT');
xlabel(xlab);
ylabel(ylab);

%% Impedance - 4 legs (seen from WT)
ZabsM3 = zeros(length(H),1);
RM3 = zeros(length(H),1);
for hh=1:length(H)
    h=H(hh);
    s = 1i*h*w;
    
%     Z1 = s*Tr3_L + imp_parallel(s*PhReact_L, 1/(s*Tuned_C));
%     Z2 = imp_cable_pi(Cable150_R+s*Cable150_L, 1/(s*Cable150_C1));
%     Z3 = s*Tr2_L + imp_cable_pi(Cable33_R+s*Cable33_L, 1/(s*Cable33_C1));
%     Z4 = imp_parallel(LCL_R2 + s*LCL_L2, 1/(s*LCL_C));% + LCL_R1 + s*LCL_L1;
% 
%     Z5 = imp_parallel(Z3+Z4,Z3+Z4)+Z2;
%     Z6 = imp_parallel(Z1,Z5)+Z2;
%     Z7 = imp_parallel(Z6,Z3+Z4);
%     Z = Z7 + Z3+Z4;

    Z1 = imp_parallel(s*PhReact_L, 1/(s*Tuned_C)) + s*Tr3_L;
    
    Z1B = imp_parallel(LCL_R1+s*LCL_L1, 1/(s*LCL_C))+ LCL_R2+s*LCL_L2+s*Tr1_L;
    Z2B = imp_parallel(Z1B, 1/(s*Cable33_C1)) + Cable33_R+s*Cable33_L;
    Z3B = imp_parallel(Z2B, 1/(s*Cable33_C2)) + s*Tr2_L;
    
    Z1C = imp_parallel(Z3B, Z3B);
    Z2C = imp_parallel(Z1C, 1/(s*Cable150_C1)) + Cable150_R+s*Cable150_L;
    Z3C = imp_parallel(Z2C, 1/(s*Cable150_C2));
    
    Z1A = imp_parallel(Z1,Z3C);
    Z2 = imp_parallel(Z1A, 1/(s*Cable150_C2)) + Cable150_R+s*Cable150_L;
    Z3 = imp_parallel(Z2, 1/(s*Cable150_C1));
    
    Z4 = imp_parallel(Z3B, Z3) + s*Tr2_L;
    Z5 = imp_parallel(Z4, 1/(s*Cable33_C2)) + Cable33_R+s*Cable33_L;
    Z6 = imp_parallel(Z5, 1/(s*Cable33_C1)) + s*Tr1_L + s*LCL_L2+LCL_R2;
    Z7 = imp_parallel(Z6, 1/(s*LCL_C));
    Z8 = LCL_R1+s*LCL_L1;
    Z = imp_parallel(Z7, Z8);
    
    Zabs = abs(Z);
    ZabsM3(hh) = Zabs;
    RM3(hh) = real(Z);
end
clear hh h Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z1B Z2B Z3B Z1A Z Zabs

figure(3)
plot(H,ZabsM3, 'LineWidth', 1); %hold on
%plot(H,RM3, 'r'); hold off
axis(view)
title('3. case: 4 WT');
xlabel(xlab);
ylabel(ylab);

% Impedance - 4 legs (seen from WT) with R proportional to harmonics
ZabsM4 = zeros(length(H),1);
RM4 = zeros(length(H),1);
for hh=1:length(H)
    h=H(hh);
    s = 1i*h*w;
    
    Z1 = imp_parallel(s*PhReact_L, 1/(s*Tuned_C)) + s*Tr3_L;
    
    Z1B = imp_parallel(h*LCL_R1+s*LCL_L1, 1/(s*LCL_C))+ h*LCL_R2+s*LCL_L2+s*Tr1_L;
    Z2B = imp_parallel(Z1B, 1/(s*Cable33_C1)) + h*Cable33_R+s*Cable33_L;
    Z3B = imp_parallel(Z2B, 1/(s*Cable33_C2)) + s*Tr2_L;
    
    Z1C = imp_parallel(Z3B, Z3B);
    Z2C = imp_parallel(Z1C, 1/(s*Cable150_C1)) + h*Cable150_R+s*Cable150_L;
    Z3C = imp_parallel(Z2C, 1/(s*Cable150_C2));
    
    Z1A = imp_parallel(Z1,Z3C);
    Z2 = imp_parallel(Z1A, 1/(s*Cable150_C2)) + h*Cable150_R+s*Cable150_L;
    Z3 = imp_parallel(Z2, 1/(s*Cable150_C1));
    
    Z4 = imp_parallel(Z3B, Z3) + s*Tr2_L;
    Z5 = imp_parallel(Z4, 1/(s*Cable33_C2)) + h*Cable33_R+s*Cable33_L;
    Z6 = imp_parallel(Z5, 1/(s*Cable33_C1)) + s*Tr1_L + s*LCL_L2+h*LCL_R2;
    Z7 = imp_parallel(Z6, 1/(s*LCL_C));
    Z8 = h*LCL_R1+s*LCL_L1;
    Z = imp_parallel(Z7, Z8);
    
    Zabs = abs(Z);
    ZabsM4(hh) = Zabs;
    RM4(hh) = real(Z);
end
clear hh h Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z Zabs

figure(4)
plot(H,ZabsM4, 'LineWidth', 1); %hold on
%plot(H,RM4, 'r'); hold off
axis(view)
title('4. case: 4 WT with R proportional to h');
xlabel(xlab);
ylabel(ylab);