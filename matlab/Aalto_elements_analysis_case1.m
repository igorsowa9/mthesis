clear all; close all; clc

%% system data:
% load('Aalto_data0to150org.mat') % only original values
LCL_L1o = 1.2; % H
LCL_R1o = 0.0; % ohm %%%%%%
LCL_L2o = 0.641; % H
LCL_R2o = 0.0; % ohm %%%%%%
LCL_Co = 1.491e-7; % F
LCL_Rco = 0; % ohm %%%%%%
Tr1_Lo = 51.568e-3; % H

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
Tuned_Co = 5.658e-6; % F
PhReact_Lo = 19.3e-3; % H

%% ------- settings -------
f = 50;
w = 2*pi*f;
H = 0:0.01:30;

%% convertion to equivalent circuit at V=150kV
Vout = 150e3;
Vlow = 8e3;
LCL_L1 = ind_equiv(LCL_L1o,f,Vlow,Vout);
LCL_R1 = res_equiv(LCL_R1o,Vlow,Vout);
LCL_L2 = ind_equiv(LCL_L2o,f,Vlow,Vout);
LCL_R2 = res_equiv(LCL_R2o,Vlow,Vout);
LCL_C  = cap_equiv(LCL_Co,f,Vlow,Vout);
LCL_Rc = res_equiv(LCL_Rco,Vlow,Vout);
Tr1_L  = ind_equiv(Tr1_Lo,f,Vlow,Vout);

Cable33_L  = ind_equiv(Cable33_Lo,f,33e3,Vout);
Cable33_R  = res_equiv(Cable33_Ro,33e3,Vout);
Cable33_C1 = cap_equiv(Cable33_C1o,f,33e3,Vout);
Cable33_C2 = cap_equiv(Cable33_C2o,f,33e3,Vout);
Tr2_L = ind_equiv(Tr2_Lo,f,33e3,Vout);

Cable150_L = ind_equiv(Cable150_Lo,f,150e3,Vout);
Cable150_R = res_equiv(Cable150_Ro,150e3,Vout);
Cable150_C1 = cap_equiv(Cable150_C1o,f,150e3,Vout);
Cable150_C2 = cap_equiv(Cable150_C2o,f,150e3,Vout);
Tr3_L = ind_equiv(Tr3_Lo,f,150e3,Vout);
Tuned_C = cap_equiv(Tuned_Co,f,150e3,Vout);
PhReact_L = ind_equiv(PhReact_Lo,f,150e3,Vout);

% if modified has to be read in Zdiff calculations !!!!!!!!!!!
save('Aalto_data0to150mod.mat','LCL_L1o','LCL_R1o','LCL_L2o','LCL_R2o',...
    'LCL_Co','LCL_Rco','Tr1_Lo','Cable33_Lo','Cable33_Ro','Cable33_C1o',...
    'Cable33_C2o','Tr2_Lo','Cable150_Lo','Cable150_Ro','Cable150_C1o',...
    'Cable150_C2o','Cable150_C1o','Tr3_Lo','Tuned_Co','PhReact_Lo',...
    'LCL_L1','LCL_R1','LCL_L2','LCL_R2',...
    'LCL_C','LCL_Rc','Tr1_L','Cable33_L','Cable33_R','Cable33_C1',...
    'Cable33_C2','Tr2_L','Cable150_L','Cable150_R','Cable150_C1',...
    'Cable150_C2','Cable150_C1','Tr3_L','Tuned_C','PhReact_L');

%% no harmonic for each element analysis
h = 13; s = 1i*h*w;
harmStr = ['Case 1 - Harmonic: ' num2str(h) '. - 1%-300% of '];
%% graphs
imp_max = 5000;
res = 0.01;
res_length = 300;
res_value0_idx = 100;

NAME = 1;
VOLT = 2;

valuesL = {{'LCL_L1'        8e3}    {'LCL_L2'       8e3}...
           {'Tr1_L'         8e3}    {'Cable33_L'    33e3}...
           {'Tr2_L'         33e3}   {'Cable150_L'   150e3}...
           {'Tr3_L'         150e3}  {'PhReact_L'    150e3}};
valuesC = {{'LCL_C'         8e3}    {'Cable33_C1'   33e3}...
           {'Cable33_C2'    33e3}   {'Cable150_C1'  150e3}...
           {'Cable150_C2'   150e3}  {'Tuned_C'      150e3}};
valuesR = {{'LCL_R1'        8e3}    {'LCL_R2'       8e3}...
           {'LCL_Rc'        8e3}    {'Cable33_R'    33e3}...
           {'Cable150_R'    150e3}  };

%% analysis for L elements varying (case 1)

ZabsVL = zeros(length(valuesL), res_length);

for vv = 1:length(valuesL)
    value0 = eval([valuesL{vv}{NAME} 'o']);
    E = res*value0 :res*value0: res_length*res*value0;

    for ee=1:length(E)
        sample = E(ee);
        ZabsVL(vv,ee) = case1_Zdiff(eval('valuesL{vv}{NAME}'),...
            ind_equiv(sample,f,valuesL{vv}{VOLT},Vout), s);
    end
    clear sample

    figure(vv)
    plot(E,ZabsVL(vv,:), 'LineWidth', 2);
    axis([0 res_length*res*value0 0 imp_max]);
    
    hx = graph2d.constantline(value0, 'LineStyle',':', 'Color','black');
    changedependvar(hx,'x');
    text(.6,.95, ['Current L value: ' num2str(value0)],...
        'FontSize',10,'Color','black','Units','normalized');
    Zrow = ZabsVL(vv,:);
    Znow = Zrow(res_value0_idx);
    hynow = graph2d.constantline(Znow, 'Color','green');
    changedependvar(hynow,'y');
    text(.6,.9,['Current Z value: ' num2str(Znow)],...
        'FontSize',10,'Color','green','Units','normalized');
    Zmax = max(Zrow);
    hymax = graph2d.constantline(Zmax, 'Color','red');
    changedependvar(hymax,'y');
    text(.6,.85,['Max. Z value: ' num2str(Zmax)],...
        'FontSize',10,'Color','red','Units','normalized');
    valmax = E(Zrow==Zmax);
    text(.65,.8,['for L=' num2str(valmax) ' (' num2str(valmax/value0*100) '% L0)'],...
        'FontSize',10,'Color','red','Units','normalized');
    
    title([harmStr valuesL{vv}{NAME}])
    xlabel('different values of element L [H]');
    ylabel('total impedance seen [ohms]'); 
end

clear E vv value value0 Znow Zmax hx hynow hymax

%% analysis for C elements varying (case 1)

ZabsVC = zeros(length(valuesC), res_length);    
       
for vv = 1:length(valuesC)
    value0 = eval([valuesC{vv}{NAME} 'o']);
    E = res*value0 :res*value0: res_length*res*value0;

    for ee=1:length(E)
        sample = E(ee);
        ZabsVC(vv,ee) = case1_Zdiff(eval('valuesC{vv}{NAME}'),...
            cap_equiv(sample,f,valuesC{vv}{VOLT},Vout), s);
    end
    clear sample

    figure(length(valuesL) + vv)
    plot(E,ZabsVC(vv,:), 'LineWidth', 2);
    axis([0 res_length*res*value0 0 imp_max]);
    
    hx = graph2d.constantline(value0, 'LineStyle',':', 'Color','black');
    changedependvar(hx,'x');
    text(.6,.95, ['Current C value: ' num2str(value0)],...
        'FontSize',10,'Color','black','Units','normalized');
    Zrow = ZabsVC(vv,:);
    Znow = Zrow(res_value0_idx);
    hynow = graph2d.constantline(Znow, 'Color','green');
    changedependvar(hynow,'y');
    text(.6,.9,['Current Z value: ' num2str(Znow)],...
        'FontSize',10,'Color','green','Units','normalized');
    Zmax = max(Zrow);
    hymax = graph2d.constantline(Zmax, 'Color','red');
    changedependvar(hymax,'y');
    text(.6,.85,['Max. Z value: ' num2str(Zmax)],...
        'FontSize',10,'Color','red','Units','normalized');
    valmax = E(Zrow==Zmax);
    text(.65,.8,['for C=' num2str(valmax) ' (' num2str(valmax/value0*100) '% C0)'],...
        'FontSize',10,'Color','red','Units','normalized');
    
    title([harmStr valuesC{vv}{NAME}]);
    xlabel('different values of element C [F]');
    ylabel('total impedance seen [ohms]');
end
clear E vv value value0 Znow Zmax hx hynow hymax       
break
%% analysis for R elements varying (case 1)

ZabsVR = zeros(length(valuesR), res_length);    
       
for vv = 1:length(valuesR)
    value0 = eval([valuesR{vv}{NAME} 'o']);
    if value0==0;
        fprintf('WARNING! no assigned value (=0) to element %s!!\n',valuesR{vv}{NAME});
        continue
    end
    E = res*value0 :res*value0: 100*res*value0;

    for ee=1:length(E)
        sample = E(ee);
        ZabsVR(vv,ee) = case1_Zdiff(eval('valuesC{vv}{NAME}'),...
            res_equiv(sample,valuesR{vv}{VOLT},Vout), s);
    end
    clear sample

    figure(length(valuesL)+length(valuesC)+  vv)
    plot(E,ZabsVR(vv,:), value0,0:imp_max);
    axis([0 5*value0 0 imp_max]);
    title([harmStr valuesR{vv}{NAME}]);
    xlabel('different values of element R [ohm]');
    ylabel('total impedance seen [ohms]');
    clear E value0
end
clear vv value value0      
    
