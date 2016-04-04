clear all; clc; close all;

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
H = res:res:30;

ymax = 50e3;
view = [0,H(length(H)),0,ymax];
% view = [9 14 0 100]; % zoom
xlab = 'harmonics order';
ylab = 'impedance (ohms)';
plots = [true false false false];

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
Tuned_C = cap_equiv(Tuned_Co,f,Vhigh,Vout);
PhReact_L = ind_equiv(PhReact_Lo,f,Vhigh,Vout);

%% Impedance - 1 leg (seen from WT)
if plots(1)==true 
    ZabsM1 = zeros(length(H),1);
    for hh=1:length(H)
        h=H(hh);
        s = 1i*h*w;

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
    plot(H,ZabsM1, 'b', 'LineWidth', 2); hold on
    hold off
    axis(view)
    title('1. case: 1 WT');
    xlabel(xlab);
    ylabel(ylab);
end


%% harmonics modal analysis

% admittance matrix
n_bus = 7;
n_h = length(H);

ZmodalM = zeros(n_h,4); % harm order, mode, modal imp abs, angle
PFmodalM = zeros(n_h,n_bus+1); % PFs for each bus for each harmonic
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
    [T,A] = eig(Y); % T - rigth eigenvector matrix
    L = inv(T); % L - left eigenvector matrix
    
    [lambda,mode] = min(abs(e));
    Zmodal = 1/lambda;
    em = e(mode);
    ang = rad2deg(angle(em));
    
    ZmodalM(hh,1) = h;
    ZmodalM(hh,2) = mode;
    ZmodalM(hh,3) = Zmodal;
    ZmodalM(hh,4) = ang;
    
    PFmodalM(hh,1) = h;
    for b=2:n_bus+1
        PFmodalM(hh,b) = abs(L(mode,b-1)*T(b-1,mode));
    end
    
end
figure(2);
plot(H,ZmodalM(:,3));
axis(view);

[Z_peak,h_crit_idx] = findpeaks(ZmodalM(:,3));
h_crit = h_crit_idx * res;

fprintf('harmonic order - critical mode - modal impedance(abs) - angle\n');
ZmodalHcrit = ZmodalM(h_crit_idx,:)

fprintf('harmonic order - participation factors for all buses\n');
PFmodalHcrit = PFmodalM(h_crit_idx,:);
PFmodalHcrit = PFmodalHcrit(:,2:end)

fprintf('greates participation factors:\n');
for f=1:length(PFmodalHcrit(:,1))
    m = PFmodalHcrit(f,:);
    fprintf('For harmonic: %f, bus: %f has greatest PF=%f\n',...
        ZmodalHcrit(f,1), find(m==max(m)), max(m));
end


break

%     % bus number choice <<<<<<<<<<<<<<<<<<<<<<
%     bus = 1;
%     I = zeros(n_bus,1);
%     I(bus) = 1;
% 
%     J = T*I;
%     U = A\J; % kolejne wiesze to mody! nie bus'y
% 
%     Tcrit = T(mode,:); % rows taken
% 
% 
%     V = L/A*T*I; % shows the harmonic voltages at certain 
%     % bus as a result of harmonic current injection
% 
% 
%     % modal current J1 linear projection of the physical 
%     % currents in the direction of the first eigenvector
%     J1 = sum(J);
% 
%     
%     Ainv_aprox = zeros(n_bus,n_bus);
%     Ainv_aprox(mode,mode) = invA(mode,mode); % only critical impedance left
% 
%     Vaprox = L*Ainv_aprox*T*I;
%     
%     invA = inv(A);
%     ZmodalM(hh) = abs(invA(mode,mode));
%     [abs(A) rad2deg(angle(A))]
% 
% figure(2)
% plot(H,ZM4)
% axis([0 60 0 500]);


