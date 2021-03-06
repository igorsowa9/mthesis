clear variables; close all; clc

% AALTO V2:
% DIFFERENT VALUES OF COMPONENTS (ORIGINAL DATA IS AT 150kV)
% UNKNOWN TYPE/TOPOLOGY OF AC TUNED FILTER ((((51.568e-3))))

%% system data:
LCL_L1o = 1.2; % H
LCL_R1o = 0.0; % ohm 
LCL_L2o = 0.641; % H
LCL_R2o = 0.0; % ohm 
LCL_Co = 1.491e-7; % F
LCL_Rco = 0; % ohm 
Tr1_Lo = 38.676e-3; % H %%%%%%%%%%%%%%%%%%%%%%%

Cable33_Lo = 18.181e-3; % H
Cable33_Ro = 0.372; % ohm
Cable33_C1o = 5.709e-8; % F
Cable33_C2o = Cable33_C1o; % F
Tr2_Lo = 38.676e-3; % H

Cable150_Lo = 1e-3; % H
Cable150_Ro = 0.056; % ohm
Cable150_C1o = 0.52e-6; % F
Cable150_C2o = Cable150_C1o; % F
Tr3_Lo = 19.338e-3; % H

%%% FILTER %%%%%%%%%%%%%%%%%%%%
Tuned_R1o = 21.63;
Tuned_C1o = 3.77e-6; % F
Tuned_L1o = 1.8e-3;
Tuned_R2o = 21.63;
Tuned_C2o = 1.89e-6; % F
Tuned_L2o = 0.883e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PhReact_Lo = 19.3e-3; % H

%% ------- settings -------
f = 50;
w = 2*pi*f;
H = 0:0.01:30;

ymax = 100e3;
view = [0,H(length(H)),0,ymax];
% view = [9 14 0 100]; % zoom
xlab = 'harmonics order';
ylab = 'impedance (ohms)';
plots = [false false true false];

load('four_cases_orginal_low.mat');

%% convertion to equivalent circuit at V=150kV
Vout = 150e3;

Vlow = 150e3;
Vmid = 150e3;
Vhigh = 150e3;

LCL_L1 = ind_equiv(LCL_L1o,f,Vlow,Vout);
LCL_R1 = res_equiv(LCL_R1o,Vlow,Vout);
LCL_L2 = ind_equiv(LCL_L2o,f,Vlow,Vout);
LCL_R2 = res_equiv(LCL_R2o,Vlow,Vout);
LCL_C  = cap_equiv(LCL_Co,f,Vlow,Vout);
LCL_Rc = res_equiv(LCL_Rco,Vlow,Vout);
Tr1_L  = ind_equiv(Tr1_Lo,f,Vlow,Vout);

Cable33_L  = ind_equiv(Cable33_Lo,f,Vmid,Vout);
Cable33_R  = res_equiv(Cable33_Ro,Vmid,Vout);
Cable33_C1 = cap_equiv(Cable33_C1o,f,Vmid,Vout);
Cable33_C2 = cap_equiv(Cable33_C2o,f,Vmid,Vout);
Tr2_L = ind_equiv(Tr2_Lo,f,Vmid,Vout);

Cable150_L = ind_equiv(Cable150_Lo,f,Vhigh,Vout);
Cable150_R = res_equiv(Cable150_Ro,Vhigh,Vout);
Cable150_C1 = cap_equiv(Cable150_C1o,f,Vhigh,Vout);
Cable150_C2 = cap_equiv(Cable150_C2o,f,Vhigh,Vout);
Tr3_L = ind_equiv(Tr3_Lo,f,Vhigh,Vout);

Tuned_R1 = res_equiv(Tuned_R1o,Vhigh,Vout);
Tuned_C1 = cap_equiv(Tuned_C1o,f,Vhigh,Vout);
Tuned_L1 = ind_equiv(Tuned_L1o,f,Vhigh,Vout);
Tuned_R2 = res_equiv(Tuned_R2o,Vhigh,Vout);
Tuned_C2 = cap_equiv(Tuned_C2o,f,Vhigh,Vout);
Tuned_L2 = ind_equiv(Tuned_L2o,f,Vhigh,Vout);

PhReact_L = ind_equiv(PhReact_Lo,f,Vhigh,Vout);

%% Impedance - 1 leg (seen from WT)
if plots(1)==true 
    ZabsM1 = zeros(length(H),1);
    for hh=1:length(H)
        h=H(hh);
        s = 1i*h*w;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ztuned1 = imp_parallel(Tuned_R1, s*Tuned_L1) + 1/(s*Tuned_C1);
        Ztuned2 = imp_parallel(Tuned_R2, s*Tuned_L2) + 1/(s*Tuned_C2);        
        Ztuned = Ztuned1+Ztuned2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Z1 = imp_parallel(s*PhReact_L, Ztuned) + s*Tr3_L;
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
    plot(H,ZabsM1, 'b', 'LineWidth', 2); hold on
%     plot(H,ZabsM1org_low, 'black--', 'LineWidth', 1) 
    hold off
    % view1 = [0,H(length(H)),0,300];
    axis(view)
    title('1. case: 1 WT');
    xlabel(xlab);
    ylabel(ylab);
end

%% Impedance - 2 legs (seen from WT)
if plots(2)==true 
    ZabsM2 = zeros(length(H),1);
    for hh=1:length(H)
        h=H(hh);
        s = 1i*h*w;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ztuned1 = imp_parallel(Tuned_R1, s*Tuned_L1) + 1/(s*Tuned_C1);
        Ztuned2 = imp_parallel(Tuned_R2, s*Tuned_L2) + 1/(s*Tuned_C2);        
        Ztuned = Ztuned1+Ztuned2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        Z1 = imp_parallel(s*PhReact_L, Ztuned) + s*Tr3_L;
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
end

%% Impedance - 4 legs (seen from WT)
if plots(3)==true 
    ZabsM3 = zeros(length(H),1);
    RM3 = zeros(length(H),1);
    for hh=1:length(H)
        h=H(hh);
        s = 1i*h*w;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Ztuned1 = Tuned_R1 + s*Tuned_L1 + 1/(s*Tuned_C1);
%         Ztuned2 = Ztuned1;
%         Ztuned3 = imp_parallel(Tuned_R2, s*Tuned_L2) + 1/(s*Tuned_C2);        
%         Ztuned = imp_parallel(imp_parallel(Ztuned1,Ztuned2),Ztuned3);
        Ztuned = imp_parallel(Tuned_R2, s*Tuned_L2) + 1/(s*Tuned_C2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         Ztuned = imp_parallel(s*Tuned_L1,Tuned_R1)+1/(s*Tuned_C1);

        Z1 = imp_parallel(s*PhReact_L, Ztuned) + s*Tr3_L;

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
    plot(H,ZabsM3, 'b', 'LineWidth', 1); hold on
%   plot(H,ZabsM3org_low, 'black', 'LineWidth', 1); hold off
    %plot(H,RM3, 'r'); hold off
    axis(view)
    title('3. case: 4 WT');
    xlabel(xlab);
    ylabel(ylab);
end

%% Impedance - 4 legs (seen from WT) with R proportional to harmonics
if plots(4)==true 
    ZabsM4 = zeros(length(H),1);
    RM4 = zeros(length(H),1);
    for hh=1:length(H)
        h=H(hh);
        s = 1i*h*w;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ztuned1 = imp_parallel(h*Tuned_R1, s*Tuned_L1) + 1/(s*Tuned_C1);
        Ztuned2 = imp_parallel(h*Tuned_R2, s*Tuned_L2) + 1/(s*Tuned_C2);        
        Ztuned = Ztuned1+Ztuned2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Z1 = imp_parallel(s*PhReact_L, Ztuned) + s*Tr3_L;

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
end
