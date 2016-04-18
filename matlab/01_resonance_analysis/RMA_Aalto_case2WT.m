clear variables; close all;

%% data

% system data:
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
H = res:res:30;
fprintf('--------------- CASE 2 (2 WT); resolution: %f --------------- \n', res);

ymax = 5e4;
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

%% Impedance - 2 legs (seen from middle LCL point)
ZabsM2 = zeros(length(H),1);
for hh=1:length(H)
    h=H(hh);
    s = 1i*h*w;

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
clear hh h Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z1B Z2B Z3B Z1A Z Zabs

figure(1)
plot(H,ZabsM2, 'b', 'LineWidth', 2);
axis(view)
title('Frequency Sweep: 2. case (2 WTs). Seen from middle LCL.');
xlabel(xlab);
ylabel(ylab);

%% harmonics modal analysis

% admittance matrix
n_bus = 11;
n_h = length(H);

ZmodalMaxM = zeros(n_h,4); % harm order, mode of max imp, modal imp abs, angle
ZmodalAllM = zeros(n_h,n_bus); % all modes (impedances)
PFmodalM = zeros(n_h,n_bus+1); % PFs for each bus for each harmonic
for hh = 1:length(H)
    h = H(hh);
    s = 1i*h*w;
    
    y1_1 = 1/(s*PhReact_L) + s*Tuned_C + 1/(s*Tr3_L);
    y2_2 = 1/(s*Tr3_L) + s*Cable150_C2 + 1/(Cable150_R+s*Cable150_L);
    y3_3 = 1/(Cable150_R+s*Cable150_L) + s*Cable150_C1 + 2*1/(s*Tr2_L);
    
    y4_4 = 1/(s*Tr2_L) + s*Cable33_C2 + 1/(Cable33_R+s*Cable33_L);
    y5_5 = 1/(Cable33_R+s*Cable33_L) + s*Cable33_C2 + 1/(s*Tr1_L);
    y6_6 = 1/(s*Tr1_L) + 1/(LCL_R2+s*LCL_L2);
    y7_7 = 1/(LCL_R2+s*LCL_L2) + s*LCL_C+LCL_Rc + 1/(LCL_R1+s*LCL_L1);
    
    y8_8 = y4_4;
    y9_9 = y5_5;
    y10_10 = y6_6;
    y11_11 = y7_7;

    y1_2 = 1/(s*Tr3_L);    
    y2_3 = 1/(Cable150_R+s*Cable150_L);
    
    y3_4 = 1/(s*Tr2_L);
    y3_8 = y3_4;
    
    y4_5 = 1/(Cable33_R+s*Cable33_L);
    y5_6 = 1/(s*Tr1_L);
    y6_7 = 1/(LCL_R2+s*LCL_L2);
    y8_9 = y4_5;
    y9_10 = y5_6;
    y10_11 = y6_7;
    
    Y = [y1_1 -y1_2 0 0 0 0 0 0 0 0 0;...
        -y1_2 y2_2 -y2_3 0 0 0 0 0 0 0 0;...
        0 -y2_3 y3_3 -y3_4 0 0 0 -y3_8 0 0 0;...
        0 0 -y3_4 y4_4 -y4_5 0 0 0 0 0 0;...
        0 0 0 -y4_5 y5_5 -y5_6 0 0 0 0 0;...
        0 0 0 0 -y5_6 y6_6 -y6_7 0 0 0 0;...
        0 0 0 0 0 -y6_7 y7_7 0 0 0 0;...
        0 0 -y3_8 0 0 0 0 y8_8 -y8_9 0 0;...
        0 0 0 0 0 0 0 -y8_9 y9_9 -y9_10 0;...
        0 0 0 0 0 0 0 0 -y9_10 y10_10 -y10_11;...
        0 0 0 0 0 0 0 0 0 -y10_11 y11_11];
    
    e = eig(Y); % eigenvalues
    [T,A] = eig(Y); % T - rigth eigenvector matrix
    L = inv(T); % L - left eigenvector matrix
    
    ZmodalAll = abs(inv(A));
    for zz = 1:n_bus
        ZmodalAllM(hh,zz) = ZmodalAll(zz,zz);
    end
    
    [lambda,mode] = min(abs(e));
    ZmodalMax = 1/lambda;
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
figure(2);
plot(H,ZmodalAllM);
axis(view);
title('Harmonics Modal Analysis (ALL modes): 2. case (2 WT)');
xlabel(xlab);
ylabel(ylab);

figure(3);
plot(H,ZmodalMaxM(:,3));
axis(view)
title('Harmonics Modal Analysis (for only MAX impedances): 2. case (2 WT)');
xlabel(xlab);
ylabel(ylab);

[Z_peak,h_crit_idx] = findpeaks(ZmodalMaxM(:,3));
h_crit = h_crit_idx * res;

fprintf('harmonic order - critical mode - modal impedance(abs) - angle\n');
ZmodalHcrit = ZmodalMaxM(h_crit_idx,:)

fprintf('harmonic order - participation factors for all buses\n');
PFmodalHcrit = PFmodalM(h_crit_idx,:);
PFmodalHcrit = PFmodalHcrit(:,2:end)

fprintf('greates participation factors:\n');
for f=1:length(PFmodalHcrit(:,1))
    M = PFmodalHcrit(f,:);
    mx = max(M);
    tmx = find(M==mx);
    omx1 = 0.99999*mx;
    omx2 = 1.00001*mx;
    omx = setdiff(intersect(find(M>omx1),find(M<omx2)),tmx);
    
    fprintf('For harmonic: %f, bus: %s has greatest PF=%f.\n',...
        ZmodalHcrit(f,1), num2str(tmx), mx);
    if ~isempty(omx)
        fprintf('\t(also same PF at buses: ');
        fprintf('%s', num2str(omx));
        fprintf(')\n');
    end
end

top_modes = unique(ZmodalMaxM(h_crit_idx,2));
top_modes_impedances = ZmodalAllM(:,top_modes);
figure(4)
plot(H,top_modes_impedances);
% axis(view)
title('Harmonics Modal Analysis (Critical modes only): 2. case (2 WT)');
legend(num2str(top_modes));
xlabel(xlab);
ylabel(ylab);
