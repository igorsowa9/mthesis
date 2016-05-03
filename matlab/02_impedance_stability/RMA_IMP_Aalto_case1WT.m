clear variables; close all; clc

%% data
% system data:
LCL_L1o = 1.2; % H
LCL_R1o = 0.0; % ohm 
LCL_L2o = 0.641; % H
LCL_R2o = 0.0; % ohm 
LCL_Co = 1.491e-7; % F
LCL_Rco = 0; % ohm 
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
res = 0.01;
H = res:res:40;
%               Zmodels f.sweep HMA     Stability-bode   nyquist others
calculate = [   true    false    true    true        false      false];

fprintf('--------------- CASE 1 (1 WT); resolution: %f --------------- \n', res);

ymax = 50e3;
view = [0,H(length(H)),0,ymax];
% view = [9 14 0 100]; % zoom
xlab = 'harmonics order';
ylab = 'impedance (ohms)';

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
    S_WTconverter = 400e6;
    [Zwt_p, Zwt_n, ~,~, Zhvdc_p, Zhvdc_n, ~,~] = convertersImpedanceModel(H,S_WTconverter);
    fprintf('Converters models: WTs and HVDC loaded! \n');

    Hf = H*f;
    figure(1); hold on;
    
    subplot(2,2,1)
    pos = semilogx(Hf,mag2db(abs(Zwt_p)), 'LineWidth', 1); hold on
    neg = semilogx(Hf,mag2db(abs(Zwt_n)), 'r', 'LineWidth', 1); hold off
    legend([pos,neg], 'positive', 'negative')
    title('Impedance of WT converter rated: 400 MW');
    xlabel('Frequency [Hz] (log axis)');
    ylabel('Impedance of WT conv. (LCL_L1 incl.) [dB]');
    clear pos neg
    grid on
    
    subplot(2,2,2)
    pos = semilogx(Hf,mag2db(abs(Zhvdc_p)), 'LineWidth', 1); hold on
    neg = semilogx(Hf,mag2db(abs(Zhvdc_n)), 'r', 'LineWidth', 1); hold off
    legend([pos,neg], 'positive', 'negative')
    title('Impedance of HVDC converter');
    xlabel('Frequency [Hz] (log axis)');
    ylabel('Impedance of WT conv. (PhR.&tunedC incl.) [dB]');
    clear pos neg
    grid on
    
    subplot(2,2,3)
    pos = semilogx(Hf,angle(rad2deg(Zwt_p)), 'LineWidth', 1); hold on
    neg = semilogx(Hf,angle(rad2deg(Zwt_n)), 'r', 'LineWidth', 1); hold off
    legend([pos,neg], 'positive', 'negative')
    title('Impedance of WT converter rated: 400 MW');
    xlabel('Frequency [Hz] (log axis)');
    ylabel('Impedance phase of WT conv. (LCL_L1 incl.) [deg]');
    clear pos neg
    grid on
    
    subplot(2,2,4)
    pos = semilogx(Hf,angle(rad2deg(Zhvdc_p)), 'LineWidth', 1); hold on
    neg = semilogx(Hf,angle(rad2deg(Zhvdc_n)), 'r', 'LineWidth', 1); hold off
    legend([pos,neg], 'positive', 'negative')
    title('Impedance of WT converter rated: 400 MW');
    xlabel('Frequency [Hz] (log axis)');
    ylabel('Impedance phase of WT conv. (LCL_L1 incl.) [deg]');
    clear pos neg
    grid on
    
    hold off
end

%% frequency sweep - 1 leg (seen from middle LCL point)
if calculate(2)==true
    fprintf('\n\n-------- Frequency Sweep --------\n\n')
    ZabsM1 = zeros(length(H),2);
    for model = 1:2 % 1-no conv.impedance models; 2-with conv.impedance models
        for hh=1:length(H)
            h=H(hh);
            s = 1i*h*w;

            if model==1
                Z1 = imp_parallel(s*PhReact_L, 1/(s*Tuned_C)) + s*Tr3_L;
            elseif model==2
                Z1 = Zhvdc_p(hh) + s*Tr3_L;
            end
            Z2 = imp_parallel(Z1, 1/(s*Cable150_C2)) + Cable150_R+s*Cable150_L;
            Z3 = imp_parallel(Z2, 1/(s*Cable150_C1)) + s*Tr2_L;
            Z4 = imp_parallel(Z3, 1/(s*Cable33_C2)) + Cable33_R+s*Cable33_L;
            Z5 = imp_parallel(Z4, 1/(s*Cable33_C1)) + s*Tr1_L + LCL_R2+s*LCL_L2;
            Z6 = imp_parallel(Z5, 1/(s*LCL_C));

%             Z6 = s*Tr1_L + LCL_R2+s*LCL_L2 + Cable33_R+s*Cable33_L +...
%                 s*Tr2_L + Cable150_R+s*Cable150_L;
            if model==1
                Z7 = LCL_R1+s*LCL_L1;
            elseif model==2
                Z7 = Zwt_p(hh);
            end
            Z = imp_parallel(Z6, Z7);

            Zabs = abs(Z);
            ZabsM1(hh,model) = Zabs;
        end
        clear hh h Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z Zabs
    end

    figure(2)
    plot(H,ZabsM1, 'LineWidth', 1)
    axis(view)
    title('F.sweep: 1. case. From middle LCL. W/out & with converter models');
    xlabel(xlab);
    ylabel(ylab);
    legend(num2str([1:2]'));

    figure(3)
    semilogy(H,ZabsM1, 'LineWidth', 1);
    title('F.sweep: 1. case. From middle LCL. W & w-out converter models');
    xlabel(xlab);
    ylabel(strcat(ylab,' (Y-log)'));
    legend(num2str([1:2]'));
    
    [zz1,ww1] = findpeaks(ZabsM1(:,1));
    [zz2,ww2] = findpeaks(ZabsM1(:,2));
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
    n_bus = 7;
    n_h = length(H);

    ZmodalMaxM = zeros(n_h,4); % harm order, mode of max imp, modal imp abs, angle
    ZmodalAllM = zeros(n_h,n_bus); % all modes (impedances)
    PFmodalM = zeros(n_h,n_bus+1); % PFs for each bus for each harmonic
    eM1 = zeros(n_bus,n_h);
    for hh = 1:length(H)
        h = H(hh);
        s = 1i*h*w;

        y11 = 1/(s*PhReact_L) + s*Tuned_C + 1/(s*Tr3_L);
        y22 = 1/(s*Tr3_L) + s*Cable150_C2 + 1/(Cable150_R+s*Cable150_L);
        y33 = 1/(Cable150_R+s*Cable150_L) + s*Cable150_C1 + 1/(s*Tr2_L);
        y44 = 1/(s*Tr2_L) + s*Cable33_C2 + 1/(Cable33_R+s*Cable33_L);
        y55 = 1/(Cable33_R+s*Cable33_L) + s*Cable33_C2 + 1/(s*Tr1_L);
        y66 = 1/(s*Tr1_L) + 1/(LCL_R2+s*LCL_L2);
        y77 = 1/(LCL_R2+s*LCL_L2) + s*LCL_C+LCL_Rc + 1/(LCL_R1+s*LCL_L1);

        y12 = 1/(s*Tr3_L);
        y23 = 1/(Cable150_R+s*Cable150_L);
        y34 = 1/(s*Tr2_L);
        y45 = 1/(Cable33_R+s*Cable33_L);
        y56 = 1/(s*Tr1_L);
        y67 = 1/(LCL_R2+s*LCL_L2);

        Y = [y11    -y12    0       0       0       0       0   ;...
            -y12    y22     -y23    0       0       0       0   ;...
            0       -y23    y33     -y34    0       0       0   ;...
            0       0       -y34    y44     -y45    0       0   ;...
            0       0       0       -y45    y55     -y56    0   ;...
            0       0       0       0       -y56    y66     -y67;...
            0       0       0       0       0       -y67    y77 ];

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
    title('Harmonics Modal Analysis (ALL modes): 1. case - no converter models');
    xlabel(xlab);
    ylabel(ylab);

    figure(5);
    subplot(2,1,1)
    plot(H,ZmodalMaxM(:,3));
    axis(view)
    title('Harmonics Modal Analysis (for only MAX impedances): 1. case - no converter models');
    xlabel(xlab);
    ylabel(ylab);

    [Z_peak,h_crit_idx] = findpeaks(ZmodalMaxM(:,3));
    h_crit = h_crit_idx * res;

    fprintf('--> harmonic order - critical mode - modal impedance(abs) - angle\n');
    ZmodalHcrit = ZmodalMaxM(h_crit_idx,:)

    fprintf('--> harmonic order - participation factors for all buses\n');
    PFmodalHcrit = PFmodalM(h_crit_idx,:);
    PFmodalHcrit = PFmodalHcrit(:,2:end)
    
    fprintf('--> Eigenvalues of critical frequencies\n');
    eM1(:,h_crit_idx)

    fprintf('--> greatest participation factors:\n');
    for ff=1:length(PFmodalHcrit(:,1))
        m = PFmodalHcrit(ff,:);
        fprintf('For harmonic: %f, bus: %f has greatest PF=%f\n',...
            ZmodalHcrit(ff,1), find(m==max(m)), max(m));
    end

    top_modes = unique(ZmodalMaxM(h_crit_idx,2));
    top_modes_impedances = ZmodalAllM(:,top_modes);
    figure(6)
    subplot(2,1,1)
    plot(H,top_modes_impedances);
    title('Harmonics Modal Analysis (Critical modes only): 1. case (1 WT)');
    legend(num2str(top_modes));
    xlabel(xlab);
    ylabel(ylab);

    %% harmonics modal analysis - model 2 (with conv. models)
    fprintf('\n\n-------- HMA - with converter models --------\n\n')
    ZmodalMaxM2 = zeros(n_h,4); % harm order, mode of max imp, modal imp abs, angle
    ZmodalAllM2 = zeros(n_h,n_bus); % all modes (impedances)
    PFmodalM2 = zeros(n_h,n_bus+1); % PFs for each bus for each harmonic
    eM2 = zeros(n_bus,n_h);
    for hh = 1:length(H)
        h = H(hh);
        s = 1i*h*w;

        y11 = 1/(Zhvdc_p(hh)) + 1/(s*Tr3_L);
        y22 = 1/(s*Tr3_L) + s*Cable150_C2 + 1/(Cable150_R+s*Cable150_L);
        y33 = 1/(Cable150_R+s*Cable150_L) + s*Cable150_C1 + 1/(s*Tr2_L);
        y44 = 1/(s*Tr2_L) + s*Cable33_C2 + 1/(Cable33_R+s*Cable33_L);
        y55 = 1/(Cable33_R+s*Cable33_L) + s*Cable33_C2 + 1/(s*Tr1_L);
        y66 = 1/(s*Tr1_L) + 1/(LCL_R2+s*LCL_L2);
        y77 = 1/(LCL_R2+s*LCL_L2) + s*LCL_C+LCL_Rc + 1/(Zwt_p(hh));

        y12 = 1/(s*Tr3_L);
        y23 = 1/(Cable150_R+s*Cable150_L);
        y34 = 1/(s*Tr2_L);
        y45 = 1/(Cable33_R+s*Cable33_L);
        y56 = 1/(s*Tr1_L);
        y67 = 1/(LCL_R2+s*LCL_L2);

        Y = [y11    -y12    0       0       0       0       0   ;...
            -y12    y22     -y23    0       0       0       0   ;...
            0       -y23    y33     -y34    0       0       0   ;...
            0       0       -y34    y44     -y45    0       0   ;...
            0       0       0       -y45    y55     -y56    0   ;...
            0       0       0       0       -y56    y66     -y67;...
            0       0       0       0       0       -y67    y77 ];

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
    title('HMA (ALL modes): 1. case - with converter models');
    xlabel(xlab);
    ylabel(ylab);

    figure(5);
    subplot(2,1,2)
    plot(H,ZmodalMaxM2(:,3));
    axis(view)
    title('HMA (for only MAX impedances): 1. case - with converter models');
    xlabel(xlab);
    ylabel(ylab);

    [Z_peak,h_crit_idx] = findpeaks(ZmodalMaxM2(:,3));
    h_crit = h_crit_idx * res;

    fprintf('--> harmonic order - critical mode - modal impedance(abs) - angle\n');
    ZmodalHcrit = ZmodalMaxM2(h_crit_idx,:)

    fprintf('--> harmonic order - participation factors for all buses\n');
    PFmodalHcrit = PFmodalM2(h_crit_idx,:);
    PFmodalHcrit = PFmodalHcrit(:,2:end)
    
    fprintf('--> Eigenvalues of critical frequencies\n');
    eM2(:,h_crit_idx)
    
    fprintf('--> greatest participation factors:\n');
    for ff=1:length(PFmodalHcrit(:,1))
        m = PFmodalHcrit(ff,:);
        fprintf('For harmonic: %f, bus: %f has greatest PF=%f\n',...
            ZmodalHcrit(ff,1), find(m==max(m)), max(m));
    end

    top_modes = unique(ZmodalMaxM2(h_crit_idx,2));
    top_modes_impedances = ZmodalAllM2(:,top_modes);
    figure(6)
    subplot(2,1,2)
    plot(H,top_modes_impedances);
    title('HMA (Critical modes only): 1. case - with converter models');
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
        Z4p = imp_parallel(Z3p,1/(s*Cable150_C1))+Cable150_R+s*Cable150_L;
        Zsp = imp_parallel(Z4p,1/(s*Cable150_C2));
        ZspM(hh) = Zsp;
        
        Z1n = imp_parallel(Zwt_n(hh),1/(s*LCL_C)) + ...
            LCL_R2 + s*(LCL_L2+Tr1_L);
        Z2n = imp_parallel(Z1n,1/(s*Cable33_C1))+Cable33_R+s*Cable33_L;
        Z3n = imp_parallel(Z2n,1/(s*Cable33_C2))+s*Tr2_L;
        Z4n = imp_parallel(Z3n,1/(s*Cable150_C1))+Cable150_R+s*Cable150_L;
        Zsn = imp_parallel(Z4n,1/(s*Cable150_C2));
        ZsnM(hh) = Zsn;
    
        Ztr3M(hh) = s*Tr3_L;
    end
    clear Z1p Z2p Z3p Z4p Zsp Z1n Z2n Z3n Z4n Zsn
    
    Zlp = Zhvdc_p + Ztr3M;
    Zln = Zhvdc_n + Ztr3M;
    
    figure(7)
    title('Bode diagram of frequency dependent impedance of source nad load');
    subplot(2,1,1)
    pos = semilogx(Hfstab,mag2db(abs(ZspM)), 'b','LineWidth', 1); hold on
    neg = semilogx(Hfstab,mag2db(abs(ZsnM)), 'b--', 'LineWidth', 1); 
    grdpos = semilogx(Hfstab,mag2db(abs(Zlp)), 'r', 'LineWidth', 1);
    grdneg = semilogx(Hfstab,mag2db(abs(Zln)), 'r--', 'LineWidth', 1); hold off
    grid on
    V = axis;
    axis([10, 2000, V(3), V(4)]);
    legend([pos,neg,grdpos,grdneg], 'pos. source', 'neg. source',...
        'pos. grid', 'neg. grid');
    xlabel('Frequency [Hz] (log axis)');
    ylabel('Impedances magnitude [dB]');
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
        'pos. grid', 'neg. grid');
    xlabel('Frequency [Hz] (log axis)');
    ylabel('Angle [deg]');
    clear pos neg grdpos grdneg V
    
end   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if calculate(5)==true
    
    fprintf('\n\n-------- Stability of the system --------\n')
    [Zwt_p, Zwt_n, TFwt_p, TFwt_n, Zhvdc_p, Zhvdc_n, TFhvdc_p, TFhvdc_n] = convertersImpedanceModel(Hstab,S_WTconverter);
        
    %% TF - Nyquist - 3 elements
    
    s = tf('s');
    TFtr3 = s*Tr3_L;
    TFr_p = minreal(TFhvdc_p+TFtr3);
    TFr_n = minreal(TFhvdc_n+TFtr3);
    
    TFw_p = TFwt_p;
    TFw_n = TFwt_n;
    
    TFc = (LCL_R2+Cable33_R+Cable150_R)+...
        (LCL_L2+Tr1_L+Cable33_L+Tr2_L+Cable150_L)*s; % series impedances only !!
    
    TF_p = minreal(TFr_p/(TFw_p+TFc));
    TF_n = minreal(TFr_n/(TFw_n+TFc));
    clear s
    
    figure(7)
    subplot(1,2,1)
    title('With only series impedances as impedance of the grid')
    nyquist(TF_p); hold on
    nyquist(TF_n); 
    circle(0,0,1); 
    axis([-nyqview nyqview -nyqview nyqview]);
    hold off
    
    % Only positive sequence
    [REp,IMp,W] = nyquist(TF_p);
    
    REp = REp(:);
    IMp = IMp(:);
    
    distancesToZero = sqrt(REp.^2+IMp.^2);
    distancesToCircle = abs(sqrt(REp.^2+IMp.^2)-1);
    
    closest_idx = intersect(find(distancesToZero>0.99),find(distancesToZero<1.01));
    if isempty(closest_idx)
        closest = min(distancesToCircle);
        closest_idx = find(abs(distancesToCircle) == closest);
    end
    WcircleHz = W(closest_idx)/(2*pi)
    WcircleOrder = WcircleHz/f
    
    %% Nyquist simplified - 2 elements
    s = tf('s');
    TF1 = minreal(minreal(imp_parallel(TFwt_p,1/(s*LCL_C))) + LCL_R2 + (LCL_L2+Tr1_L)*s);
    TF2 = minreal(minreal(imp_parallel(TF1,1/(s*Cable33_C1))) + Cable33_R + Cable33_L*s);
    TF3 = minreal(minreal(imp_parallel(TF2,1/(s*Cable33_C2))) + Tr2_L*s);
    TF4 = minreal(minreal(imp_parallel(TF3,1/(s*Cable150_C1))) + Cable150_R + Cable150_L*s);
    TFs = TF4;
    %     TFs = minreal(imp_parallel(TF4,1/(s*Cable150_C2)));
    
    % !!!!!!!!!!!!!!!!!!!!! too many zeros and poles !!!!!!!!!!!!!!!!!!!!!!
    
    TFtr3 = s*Tr3_L;
    TFl_p = minreal(TFhvdc_p+TFtr3);
    
    TF_p = minreal(TFs/TFl_p);
    
    clear s
    
    figure(7)
    subplot(1,2,2)
    nyquist(TF_p); hold on
    circle(0,0,1); 
%     axis([-nyqview nyqview -nyqview nyqview]);
    hold off
    
    %% Nyquist with Zin Zout
    if calculate(6)==true
        % ABCD two-port matrix
        s = tf('s');
        E = [1 0;s*LCL_C 1]*[1 LCL_R2+s*(LCL_L2+Tr1_L);0 1]*...
            [1 0;s*Cable33_C1 1]*[1 Cable33_R+s*Cable33_L;0 1]*...
            [1 0;s*Cable33_C2 1]*[1 s*Tr2_L;0 1]*[1 0;s*Cable150_C1 1]*...
            [1 Cable150_R+s*Cable150_L;0 1]*[1 0;s*Cable150_C2 1];
        A = E(1,1);
        B = E(1,2);
        C = E(2,1);
        D = E(2,2);

        Z11 = A/C;
        Z12 = (A*D-B*C)/C;
        Z21 = 1/C;
        Z22 = D/C;

        TFc_in = Z11 - (Z12*Z21/(Z22+TFhvdc_p)); % System input impedance with
    %     impedance of HVDC as load
        TFc_out = Z22 - (Z12*Z21/(Z11+TFwt_p)); % System output impedance with
    %     impedance of WT as a load

        TF_p_in = minreal(TFr_p/(TFw_p+TFc_in));
        TF_n_in = minreal(TFr_n/(TFw_n+TFc_in));
        TF_p_out = minreal(TFr_p/(TFw_p+TFc_out));
        TF_n_out = minreal(TFr_n/(TFw_n+TFc_out));
        clear s

        figure(8)
        title('with Z of the grid as Zin - HVDC impedance as "load"')
        nyquist(TF_p_in); hold on
        nyquist(TF_n_in); 
        circle(0,0,1); 
        axis([-1.5 1.5 -1.5 1.5]);
        hold off

        figure(9)
        title('with Z of the grid as Zout - WT impedance as "load"')
        nyquist(TF_p_out); hold on
        nyquist(TF_n_out); 
        circle(0,0,1); 
        axis([-1.5 1.5 -1.5 1.5]);
        hold off
    end
    
%     figure(8)
%     nyquist(TF_p,{0.001,W(closest_idx)}); hold on
%     circle(0,0,1); 
%     axis([-1.5 1.5 -1.5 1.5]);
%     hold off
    
%     figure(9)
%     nyquist(TF_p,{0.001,3412}); hold on
%     circle(0,0,1); 
%     axis([-1.5 1.5 -1.5 1.5]);
%     hold off
    
%     TR3_Z = Tr3_L * 1i*Hstab*w;
%     Zr_p = Zhvdc_p + TR3_Z; 
%     Zw_p = Zwt_p;
%     Zr_n = Zhvdc_n + TR3_Z; 
%     Zw_n = Zwt_n; 
%       
%     Zc = (LCL_R2+Cable33_R+Cable150_R)*ones(1,length(Hstab)) +...
%         (LCL_L2+Tr1_L+Cable33_L+Tr2_L+Cable150_L) * 1i*Hstab*w; % series impedances only !!
% 
%     nyq_p = Zr_p./(Zw_p+Zc);
%     nyq_n = Zr_n./(Zw_n+Zc);
%     figure(8)
%     plot(nyq_p, 'b'); hold on
%     plot(nyq_n, 'r');
%     grid on
%     circle(0,0,1); 
% %     axis([-1.5 1.5 -1.5 1.5]);
%     hold off
    
end

