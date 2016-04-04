clear all; clc; close all;

%% data

Rsys = 0.04;
Xsys = 0.3;
RL = 0.835;
XL = 4; 

G1 = 0;
G2 = G1;
B1 = 0.0013;
B2 = B1;

%% frequency sweep

H = 0.01:0.01:60;

% from bus1
ZM1 = zeros(1,length(H));
for hh = 1:length(H)
    h = H(hh);
    
    Z1 = imp_parallel((1/(G2+h*1i*B2) + 2*(RL+h*1i*XL)), 1/(G1+h*1i*B1));
    Z = imp_parallel(Z1, (Rsys+h*1i*Xsys));
    ZM1(hh)=abs(Z);
end
clear Z1 Z hh h;

% from bus2
ZM2 = zeros(1,length(H));
for hh = 1:length(H)
    h = H(hh);
    
    Z1 = 1/(G2+h*1i*B2) + RL+h*1i*XL;
    Z2 = imp_parallel(Rsys+h*1i*Xsys, 1/(G1+h*1i*B1)) + RL+h*1i*XL;
    Z = imp_parallel(Z1, Z2);
    ZM2(hh)=abs(Z);
end
clear Z1 Z2 Z hh h;

% from bus3
ZM3 = zeros(1,length(H));
for hh = 1:length(H)
    h = H(hh);
    
    Z1 = 1/(G2+h*1i*B2);
    Z2 = imp_parallel(Rsys+h*1i*Xsys, 1/(G1+h*1i*B1)) + 2*(RL+h*1i*XL);
    Z = imp_parallel(Z1, Z2);
    ZM3(hh)=abs(Z);
end
clear Z1 Z2 Z hh h;

figure(1); hold on
plot(H, ZM1);
plot(H, ZM2, 'r');
plot(H, ZM3, 'g');
axis([0 60 0 500]);
hold off

%% harmonics modal analysis

% admittance matrix
n_bus = 3;
res = 0.001;

H = res:res:60;
n_h = length(H);

% mode 1
ZmodalM = zeros(n_h,4); % harm order, mode, modal imp abs, angle
PFmodalM = zeros(n_h,n_bus+1); % PFs for buses
for hh = 1:length(H)
    h = H(hh);
    
    y11 = 1/(Rsys+h*1i*Xsys) + G1+h*1i*B1 + 1/(RL+h*1i*XL);
    y22 = 1/(RL+h*1i*XL) + 1/(RL+h*1i*XL);
    y33 = 1/(RL+h*1i*XL) + G2+h*1i*B2;

    y12 = 1/(RL+h*1i*XL);
    y23 = 1/(RL+h*1i*XL);

    Y = [y11 -y12 0;...
        -y12 y22 -y23;...
        0 -y23 y33];
    
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
    PFmodalM(hh,2) = abs(L(1,mode)*T(mode,1));
    PFmodalM(hh,3) = abs(L(2,mode)*T(mode,2));
    PFmodalM(hh,4) = abs(L(3,mode)*T(mode,3));
    
end
figure(2);
plot(H,ZmodalM(:,3));
axis([0 60 0 500]);

[Z_peak,h_crit_idx] = findpeaks(ZmodalM(:,3));
h_crit = h_crit_idx * res;

ZmodalM(h_crit_idx,:)
PFmodalM(h_crit_idx,:)

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


