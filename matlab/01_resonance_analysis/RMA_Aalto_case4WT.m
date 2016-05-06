clear variables; close all; clc

%% data

%% system data:
LCL_L1o = 1.2; % H
LCL_R1o = 0.0; % ohm 
LCL_L2o = 0.641; % H
LCL_R2o = 0.0; % ohm 
LCL_Co = 1.491e-7; % F
LCL_Rco = 0; % ohm 
Tr1_Lo = 51.568e-3; % H %%%%%%%%%%%%%%%%%%%%%%%

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
res = 0.01;
H = res:res:40;
%               Zmodels f.sweep HMA     Stability-bode
calculate = [   true    true    true    true];

fprintf('--------------- CASE 3 (4 WT); resolution: %f --------------- \n', res);

ymax = 50e3;
view = [0,H(length(H)),0,ymax];
% view = [9 14 0 100]; % zoom
xlab = 'Harmonics Order';
ylab = 'Impedance [Ohm]';

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
Tuned_C = cap_equiv(Tuned_Co,f,Vhigh,Vout);
PhReact_L = ind_equiv(PhReact_Lo,f,Vhigh,Vout);

%% loading converters models
if calculate(1)==true
    S_WTconverter = 100e6;
    [Zwt_p, Zwt_n, ~,~, Zhvdc_p, Zhvdc_n, ~,~] = convertersImpedanceModel(H,S_WTconverter);
    fprintf('Converters models: WTs and HVDC loaded! \n');

    Hf = H*f;
    figure(1); hold on;
    
    subplot(2,2,1)
    pos = semilogx(Hf,mag2db(abs(Zwt_p)), 'LineWidth', 1); hold on
    neg = semilogx(Hf,mag2db(abs(Zwt_n)), 'r', 'LineWidth', 1); hold off
    legend([pos,neg], 'positive', 'negative','Location','SouthEast');
    title('Impedance of aggregated WT converter (100MW)');
    xlabel('Frequency [Hz]');
    ylabel('Impedance magnitude [dB]');
    V = axis;
    axis([Hf(1), Hf(end), V(3), V(4)]);
    clear pos neg V
    grid on
    
    subplot(2,2,2)
    pos = semilogx(Hf,mag2db(abs(Zhvdc_p)), 'LineWidth', 1); hold on
    neg = semilogx(Hf,mag2db(abs(Zhvdc_n)), 'r', 'LineWidth', 1); hold off
    legend([pos,neg], 'positive', 'negative','Location','SouthEast');
    title('Impedance of HVDC converter');
    xlabel('Frequency [Hz]');
    ylabel('Impedance magnitude [dB]');
    V = axis;
    axis([Hf(1), Hf(end), V(3), V(4)]);
    clear pos neg V
    grid on
    
    subplot(2,2,3)
    pos = semilogx(Hf,angle(rad2deg(Zwt_p)), 'LineWidth', 1); hold on
    neg = semilogx(Hf,angle(rad2deg(Zwt_n)), 'r', 'LineWidth', 1); hold off
    legend([pos,neg], 'positive', 'negative')
    title('Impedance of aggregated WT converter (100MW)');
    xlabel('Frequency [Hz]');
    ylabel('Impedance angle [deg.]');
    V = axis;
    axis([Hf(1), Hf(end), V(3), V(4)]);
    clear pos neg V
    grid on
    
    subplot(2,2,4)
    pos = semilogx(Hf,angle(rad2deg(Zhvdc_p)), 'LineWidth', 1); hold on
    neg = semilogx(Hf,angle(rad2deg(Zhvdc_n)), 'r', 'LineWidth', 1); hold off
    legend([pos,neg], 'positive', 'negative')
    title('Impedance of HVDC converter');
    xlabel('Frequency [Hz]');
    ylabel('Impedance angle [deg.]');
    V = axis;
    axis([Hf(1), Hf(end), V(3), V(4)]);
    clear pos neg V
    grid on
    
    hold off
end

%% Impedance - 4 legs (seen from middle LCL point)
if calculate(2)==true
    fprintf('\n\n-------- Frequency Sweep --------\n\n')
    ZabsM3 = zeros(length(H),1);
    for model = 1:2;
        for hh=1:length(H)
            h=H(hh);
            s = 1i*h*w;
            
            if model==1
                Z1 = imp_parallel(s*PhReact_L, 1/(s*Tuned_C)) + s*Tr3_L;
            elseif model==2
                Z1 = Zhvdc_p(hh) + s*Tr3_L;
            end
            
            if model==1
                Z1B = imp_parallel(LCL_R1+s*LCL_L1,1/(s*LCL_C))+LCL_R2+s*LCL_L2+s*Tr1_L;
            elseif model==2
                Z1B = imp_parallel(Zwt_p(hh),1/(s*LCL_C))+LCL_R2+s*LCL_L2+s*Tr1_L;
            end
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
            if model==1
                Z8 = LCL_R1+s*LCL_L1;
            elseif model==2
                Z8 = Zwt_p(hh);
            end
            Z = imp_parallel(Z7, Z8);

            Zabs = abs(Z);
            ZabsM3(hh,model) = Zabs;
        end
    end
    clear hh h Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z1B Z2B Z3B Z1A Z1C Z2C Z3C Z Zabs

    figure(2)
    model1 = plot(H,ZabsM3(:,1), 'LineWidth', 1,...
        'Color',[4/255, 164/255, 217/255]); hold on
    model2 = plot(H,ZabsM3(:,2), 'LineWidth', 1,...
        'Color',[226/255, 104/255, 24/255]);
    axis(view)
    title('Frequency sweep - seen from middle LCL filter point');
    xlabel(xlab);
    ylabel(ylab);
    legend([model1,model2], 'no conv. models', 'with conv. models');
    hold off
    clear model1 model2;
    
    figure(3)
    model1 = semilogy(H,ZabsM3(:,1), 'LineWidth', 1,...
        'Color',[4/255, 164/255, 217/255]); hold on
    model2 = semilogy(H,ZabsM3(:,2), 'LineWidth', 1,...
        'Color',[226/255, 104/255, 24/255]);
    title('Frequency sweep - seen from middle LCL filter point');
    xlabel(xlab);
    ylabel(strcat(ylab,' (log axis)'));
    legend([model1,model2], 'no conv. models', 'with conv. models');
    hold off
    clear model1 model2;
    
    [zz1,ww1] = findpeaks(ZabsM3(:,1));
    [zz2,ww2] = findpeaks(ZabsM3(:,2));
    fprintf('\nNo converter model:')
    fprintf('\n- for h.order: %f, peak impedance: %f',[ww1'*res; zz1'])
    fprintf('\nWith converter model:')
    fprintf('\n- for h.order: %f, peak impedance: %f',[ww2'*res; zz2'])
    
    clear zz1 ww1 zz2 ww2
end

%% harmonics modal analysis - model 1 (no conv. models)
if calculate(3)==true
    fprintf('\n\n-------- HMA - without converter models --------\n\n')
    % admittance matrix
    n_bus = 20;
    n_h = length(H);

    ZmodalMaxM = zeros(n_h,4); % harm order, mode of max imp, modal imp abs, angle
    ZmodalAllM = zeros(n_h,n_bus); % all modes (impedances)
    PFmodalM = zeros(n_h,n_bus+1); % PFs for each bus for each harmonic
    eM1 = zeros(n_bus,n_h);
    for hh = 1:length(H)
        h = H(hh);
        s = 1i*h*w;

        y1_1 = 1/(s*PhReact_L) + s*Tuned_C + 1/(s*Tr3_L);
        y2_2 = 1/(s*Tr3_L) + 2*s*Cable150_C2 + 2*1/(Cable150_R+s*Cable150_L);
        y3_3 = 1/(Cable150_R+s*Cable150_L) + s*Cable150_C1 + 2*1/(s*Tr2_L);

        y4_4 = 1/(s*Tr2_L) + s*Cable33_C2 + 1/(Cable33_R+s*Cable33_L);
        y5_5 = 1/(Cable33_R+s*Cable33_L) + s*Cable33_C2 + 1/(s*Tr1_L);
        y6_6 = 1/(s*Tr1_L) + 1/(LCL_R2+s*LCL_L2);
        y7_7 = 1/(LCL_R2+s*LCL_L2) + s*LCL_C+LCL_Rc + 1/(LCL_R1+s*LCL_L1);

        y8_8 = y4_4;
        y9_9 = y5_5;
        y10_10 = y6_6;
        y11_11 = y7_7;

        y12_12 = y3_3;

        y13_13 = y4_4;
        y14_14 = y5_5;
        y15_15 = y6_6;
        y16_16 = y7_7;

        y17_17 = y4_4;
        y18_18 = y5_5;
        y19_19 = y6_6;
        y20_20 = y7_7;

        y1_2 = 1/(s*Tr3_L);

        y2_3 = 1/(Cable150_R+s*Cable150_L);
        y2_12 = y2_3;

        y3_4 = 1/(s*Tr2_L);
        y3_8 = y3_4;
        y12_13 = y3_4;
        y12_17 = y3_4;

        y4_5 = 1/(Cable33_R+s*Cable33_L);
        y5_6 = 1/(s*Tr1_L);
        y6_7 = 1/(LCL_R2+s*LCL_L2);
        y8_9 = y4_5;
        y9_10 = y5_6;
        y10_11 = y6_7;
        y13_14 = y4_5;
        y14_15 = y5_6;
        y15_16 = y6_7;
        y17_18 = y4_5;
        y18_19 = y5_6;
        y19_20 = y6_7;

        Y = [y1_1 -y1_2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            -y1_2 y2_2 -y2_3 0 0 0 0 0 0 0 0 -y2_12 0 0 0 0 0 0 0 0;...
            0 -y2_3 y3_3 -y3_4 0 0 0 -y3_8 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 -y3_4 y4_4 -y4_5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 -y4_5 y5_5 -y5_6 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 -y5_6 y6_6 -y6_7 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 0 -y6_7 y7_7 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 -y3_8 0 0 0 0 y8_8 -y8_9 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0 -y8_9 y9_9 -y9_10 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0 0 -y9_10 y10_10 -y10_11 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 -y10_11 y11_11 0 0 0 0 0 0 0 0 0;...
            0 -y2_12 0 0 0 0 0 0 0 0 0 y12_12 -y12_13 0 0 0 -y12_17 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 -y12_13 y13_13 -y13_14 0 0 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 -y13_14 y14_14 -y14_15 0 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 -y14_15 y15_15 -y15_16 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 -y15_16 y16_16 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 -y12_17 0 0 0 0 y17_17 -y17_18 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -y17_18 y18_18 -y18_19 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -y18_19 y19_19 -y19_20;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -y19_20 y20_20];

        e = eig(Y); % eigenvalues
        eM1(:,hh) = e; 
        [T,A] = eig(Y); % T - rigth eigenvector matrix
        L = inv(T); % L - left eigenvector matrix

        ZmodalAll = abs(inv(A));
        for zz = 1:n_bus
            ZmodalAllM(hh,zz) = ZmodalAll(zz,zz);
        end    

        [lambdaMin,mode] = min(abs(e));
        ZmodalMax = 1/lambdaMin;
        em = e(mode);
        ang = rad2deg(angle(em));

        ZmodalMaxM(hh,1) = h;
        ZmodalMaxM(hh,2) = mode;
        ZmodalMaxM(hh,3) = ZmodalMax;
        ZmodalMaxM(hh,4) = ang;

        PFmodalM(hh,1) = h;
        for b=2:n_bus+1
            PFmodalM(hh,b) = abs(L(mode,b-1)*T(b-1,mode));
        end

    end
    figure(4);
    subplot(2,1,1)
    plot(H,ZmodalAllM);
    axis(view);
    title('HMA (all modes) - no conv. models');
    xlabel(xlab);
    ylabel(ylab);
    legend(num2str([1:n_bus]'),'Location','NorthEast');

    figure(5);
    subplot(2,1,1)
    plot(H,ZmodalMaxM(:,3));
    axis(view)
    title('HMA (only max. modes) - no conv. models');    
    xlabel(xlab);
    ylabel(ylab);

    [Z_peak,h_crit_idx] = findpeaks(ZmodalMaxM(:,3));
    h_crit = h_crit_idx * res;

    fprintf('harmonic order - critical mode - modal impedance(abs) - angle\n');
    ZmodalHcrit = ZmodalMaxM(h_crit_idx,:)

    fprintf('harmonic order - participation factors for all buses\n');
    PFmodalHcrit = PFmodalM(h_crit_idx,:);
    PFmodalHcrit = PFmodalHcrit(:,2:end)
    
    fprintf('--> Eigenvalues of critical frequencies\n');
    eM1(:,h_crit_idx)

    fprintf('greates participation factors:\n');
    for ff=1:length(PFmodalHcrit(:,1))
        M = PFmodalHcrit(ff,:);
        mx = max(M);
        tmx = find(M==mx);
        omx1 = 0.99999*mx;
        omx2 = 1.00001*mx;
        omx = setdiff(intersect(find(M>omx1),find(M<omx2)),tmx);

        fprintf('For harmonic: %f, bus: %s has greatest PF=%f.\n',...
            ZmodalHcrit(ff,1), num2str(tmx), mx);
        if ~isempty(omx)
            fprintf('\t(also same PF at buses: ');
            fprintf('%s', num2str(omx));
            fprintf(')\n');
        end
    end

    top_modes = unique(ZmodalMaxM(h_crit_idx,2));
    top_modes_impedances = ZmodalAllM(:,top_modes);
    figure(6)
    subplot(2,1,1)
    plot(H,top_modes_impedances);
    title('HMA (critical modes only) - with conv. models');
    legend(num2str(top_modes));
    xlabel(xlab);
    ylabel(ylab);

    %% harmonics modal analysis - model 2 (no conv. models)
    fprintf('\n\n-------- HMA - with converter models --------\n\n')
    ZmodalMaxM2 = zeros(n_h,4); % harm order, mode of max imp, modal imp abs, angle
    ZmodalAllM2 = zeros(n_h,n_bus); % all modes (impedances)
    PFmodalM2 = zeros(n_h,n_bus+1); % PFs for each bus for each harmonic
    eM2 = zeros(n_bus,n_h);
    for hh = 1:length(H)
        h = H(hh);
        s = 1i*h*w;

        y1_1 = 1/(Zhvdc_p(hh)) + 1/(s*Tr3_L);
        y2_2 = 1/(s*Tr3_L) + 2*s*Cable150_C2 + 2*1/(Cable150_R+s*Cable150_L);
        y3_3 = 1/(Cable150_R+s*Cable150_L) + s*Cable150_C1 + 2*1/(s*Tr2_L);

        y4_4 = 1/(s*Tr2_L) + s*Cable33_C2 + 1/(Cable33_R+s*Cable33_L);
        y5_5 = 1/(Cable33_R+s*Cable33_L) + s*Cable33_C2 + 1/(s*Tr1_L);
        y6_6 = 1/(s*Tr1_L) + 1/(LCL_R2+s*LCL_L2);
        y7_7 = 1/(LCL_R2+s*LCL_L2) + s*LCL_C+LCL_Rc + 1/(Zwt_p(hh));

        y8_8 = y4_4;
        y9_9 = y5_5;
        y10_10 = y6_6;
        y11_11 = y7_7;

        y12_12 = y3_3;

        y13_13 = y4_4;
        y14_14 = y5_5;
        y15_15 = y6_6;
        y16_16 = y7_7;

        y17_17 = y4_4;
        y18_18 = y5_5;
        y19_19 = y6_6;
        y20_20 = y7_7;

        y1_2 = 1/(s*Tr3_L);

        y2_3 = 1/(Cable150_R+s*Cable150_L);
        y2_12 = y2_3;

        y3_4 = 1/(s*Tr2_L);
        y3_8 = y3_4;
        y12_13 = y3_4;
        y12_17 = y3_4;

        y4_5 = 1/(Cable33_R+s*Cable33_L);
        y5_6 = 1/(s*Tr1_L);
        y6_7 = 1/(LCL_R2+s*LCL_L2);
        y8_9 = y4_5;
        y9_10 = y5_6;
        y10_11 = y6_7;
        y13_14 = y4_5;
        y14_15 = y5_6;
        y15_16 = y6_7;
        y17_18 = y4_5;
        y18_19 = y5_6;
        y19_20 = y6_7;

        Y = [y1_1 -y1_2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            -y1_2 y2_2 -y2_3 0 0 0 0 0 0 0 0 -y2_12 0 0 0 0 0 0 0 0;...
            0 -y2_3 y3_3 -y3_4 0 0 0 -y3_8 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 -y3_4 y4_4 -y4_5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 -y4_5 y5_5 -y5_6 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 -y5_6 y6_6 -y6_7 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 0 -y6_7 y7_7 0 0 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 -y3_8 0 0 0 0 y8_8 -y8_9 0 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0 -y8_9 y9_9 -y9_10 0 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0 0 -y9_10 y10_10 -y10_11 0 0 0 0 0 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 -y10_11 y11_11 0 0 0 0 0 0 0 0 0;...
            0 -y2_12 0 0 0 0 0 0 0 0 0 y12_12 -y12_13 0 0 0 -y12_17 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 -y12_13 y13_13 -y13_14 0 0 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 -y13_14 y14_14 -y14_15 0 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 -y14_15 y15_15 -y15_16 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 -y15_16 y16_16 0 0 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 -y12_17 0 0 0 0 y17_17 -y17_18 0 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -y17_18 y18_18 -y18_19 0;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -y18_19 y19_19 -y19_20;...
            0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -y19_20 y20_20];

        e = eig(Y); % eigenvalues
        eM2(:,hh) = e; 
        [T,A] = eig(Y); % T - rigth eigenvector matrix
        L = inv(T); % L - left eigenvector matrix

        ZmodalAll = abs(inv(A));
        for zz = 1:n_bus
            ZmodalAllM2(hh,zz) = ZmodalAll(zz,zz);
        end    

        [lambdaMin,mode] = min(abs(e));
        ZmodalMax = 1/lambdaMin;
        em = e(mode);
        ang = rad2deg(angle(em));

        ZmodalMaxM2(hh,1) = h;
        ZmodalMaxM2(hh,2) = mode;
        ZmodalMaxM2(hh,3) = ZmodalMax;
        ZmodalMaxM2(hh,4) = ang;

        PFmodalM2(hh,1) = h;
        for b=2:n_bus+1
            PFmodalM2(hh,b) = abs(L(mode,b-1)*T(b-1,mode));
        end

    end
    figure(4);
    subplot(2,1,2)
    plot(H,ZmodalAllM2);
    axis(view);
    title('HMA (all modes) - with conv. models');
    xlabel(xlab);
    ylabel(ylab);
    legend(num2str([1:n_bus]'),'Location','NorthEast');

    figure(5);
    subplot(2,1,2)
    plot(H,ZmodalMaxM2(:,3));
    axis(view)
    title('HMA (only max. modes) - with conv. models');    
    xlabel(xlab);
    ylabel(ylab);

    [Z_peak,h_crit_idx] = findpeaks(ZmodalMaxM2(:,3));
    h_crit = h_crit_idx * res;

    fprintf('harmonic order - critical mode - modal impedance(abs) - angle\n');
    ZmodalHcrit = ZmodalMaxM2(h_crit_idx,:)

    fprintf('harmonic order - participation factors for all buses\n');
    PFmodalHcrit = PFmodalM2(h_crit_idx,:);
    PFmodalHcrit = PFmodalHcrit(:,2:end)
    
    fprintf('--> Eigenvalues of critical frequencies\n');
    eM2(:,h_crit_idx)

    fprintf('greates participation factors:\n');
    for ff=1:length(PFmodalHcrit(:,1))
        M = PFmodalHcrit(ff,:);
        mx = max(M);
        tmx = find(M==mx);
        omx1 = 0.99999*mx;
        omx2 = 1.00001*mx;
        omx = setdiff(intersect(find(M>omx1),find(M<omx2)),tmx);

        fprintf('For harmonic: %f, bus: %s has greatest PF=%f.\n',...
            ZmodalHcrit(ff,1), num2str(tmx), mx);
        if ~isempty(omx)
            fprintf('\t(also same PF at buses: ');
            fprintf('%s', num2str(omx));
            fprintf(')\n');
        end
    end

    top_modes = unique(ZmodalMaxM2(h_crit_idx,2));
    top_modes_impedances = ZmodalAllM2(:,top_modes);
    figure(6)
    subplot(2,1,2)
    plot(H,top_modes_impedances);
    title('HMA (critical modes only) - with conv. models');
    legend(num2str(top_modes));
    xlabel(xlab);
    ylabel(ylab);
end

%% Stability
% PCC of the system to analyze stability is between TR3 and Cable 150kV.
% Thus TR3 belongs to Zr/Zl. LCL_L1 included in Zw. Zc is the rest of the grid.
% Zc (inner grid) has to be improved because now it includes only series
% impedances.
nyqview = 1.5;
Hstab = 0.2:0.01:40; % order
Hfstab = Hstab*f; % Hz

if calculate(4)==true
    fprintf('\n\n-------- Stability of the system --------\n')
    [Zwt_p, Zwt_n, ~, ~, Zhvdc_p, Zhvdc_n, ~, ~] = convertersImpedanceModel(Hstab,S_WTconverter);
    
    %% BODE
    ZspM = zeros(1,length(Hstab));
    ZsnM = zeros(1,length(Hstab));
    Ztr3M = zeros(1,length(Hstab));
    for hh = 1:length(Hstab)
        h = Hstab(hh);
        s = 1i*h*w;
        
        Z1p = imp_parallel(Zwt_p(hh),1/(s*LCL_C)) + ...
            LCL_R2 + s*(LCL_L2+Tr1_L);
        Z2p = imp_parallel(Z1p,1/(s*Cable33_C1))+Cable33_R+s*Cable33_L;
        Z3p = imp_parallel(Z2p,1/(s*Cable33_C2))+s*Tr2_L;
        Z3pB = imp_parallel(Z3p,Z3p); % incl. second branch
        Z4p = imp_parallel(Z3pB,1/(s*Cable150_C1))+Cable150_R+s*Cable150_L;
        Z5p = imp_parallel(Z4p,1/(s*Cable150_C2));
        Zsp = imp_parallel(Z5p,Z5p);
        ZspM(hh) = Zsp;
        
        Z1n = imp_parallel(Zwt_n(hh),1/(s*LCL_C)) + ...
            LCL_R2 + s*(LCL_L2+Tr1_L);
        Z2n = imp_parallel(Z1n,1/(s*Cable33_C1))+Cable33_R+s*Cable33_L;
        Z3n = imp_parallel(Z2n,1/(s*Cable33_C2))+s*Tr2_L;
        Z3nB = imp_parallel(Z3n,Z3n); % incl. second branch
        Z4n = imp_parallel(Z3nB,1/(s*Cable150_C1))+Cable150_R+s*Cable150_L;
        Z5n = imp_parallel(Z4n,1/(s*Cable150_C2));
        Zsn = imp_parallel(Z5n,Z5n);
        ZsnM(hh) = Zsn;
    
        Ztr3M(hh) = s*Tr3_L;
    end
    clear Z1p Z2p Z3p Z4p Zsp Z1n Z2n Z3n Z4n Zsn Z3pB Z3nB
    
    Zlp = Zhvdc_p + Ztr3M;
    Zln = Zhvdc_n + Ztr3M;
    
    figure(7)
    subplot(2,1,1)
    pos = semilogx(Hfstab,mag2db(abs(ZspM)), 'b','LineWidth', 1); hold on
    neg = semilogx(Hfstab,mag2db(abs(ZsnM)), 'b--', 'LineWidth', 1); 
    grdpos = semilogx(Hfstab,mag2db(abs(Zlp)), 'r', 'LineWidth', 1);
    grdneg = semilogx(Hfstab,mag2db(abs(Zln)), 'r--', 'LineWidth', 1); hold off
    grid on
        title('Bode diagram of frequency dependent impedance of "source" and "grid"');
    V = axis;
    axis([10, 2000, V(3), V(4)]);
    legend([pos,neg,grdpos,grdneg], 'pos. source', 'neg. source',...
        'pos. grid', 'neg. grid', 'Location', 'EastOutside');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [dB]');
    clear pos neg grdpos grdneg V
    
    subplot(2,1,2)
    pos = semilogx(Hfstab,rad2deg(angle(ZspM)), 'b','LineWidth', 1); hold on
    neg = semilogx(Hfstab,rad2deg(angle(ZsnM)), 'b--', 'LineWidth', 1); 
    grdpos = semilogx(Hfstab,rad2deg(angle(Zlp)), 'r', 'LineWidth', 1);
    grdneg = semilogx(Hfstab,rad2deg(angle(Zln)), 'r--', 'LineWidth', 1); hold off
    grid on
    V = axis;
    axis([10, 2000, V(3), V(4)]);
    legend([pos,neg,grdpos,grdneg], 'pos. source', 'neg. source',...
        'pos. grid', 'neg. grid', 'Location', 'EastOutside');
    xlabel('Frequency [Hz]');
    ylabel('Angle [deg]');
    clear pos neg grdpos grdneg V
    
end   