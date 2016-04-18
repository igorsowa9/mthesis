clear variables; close all

%% WT converter
Sr = 400e6;

% V0 = 1500; % dc bus (original) !! WHY?
V0 = 150e3*sqrt(2)/0.9; 
% V1 = 563; % ph. volt. ampl. (original) !! WHY?
Vll = 150e3;
V1 = Vll/sqrt(3);
% I1 = 236e3; % ph. curr. ampl. (original) (S/sqrt(3)/Vll)
I1 = Sr/(sqrt(3)*Vll);

L = 1.2;%0.526e-6; % ph. imp.
% Lf = 0.5e-6; % filter imp.
% Cf = 0.207; % filter cap.
% Rf = 2.75e-3; % filter res.

% Hi(s) current control compensator
Kp_i = 0.44e-6;
Ki_i = 0.55e-3;

Kd_i = 22.5e-6; % current decoupling coefficient
fi = 250; % current control bandwidth

% Hp(s) PLL compensator
Kp_p = 0.239;
Ki_p = 45;
fpll = 30; % PLL bandwidth

f1 = 50;
w1 = 2*pi*f1;
res = 0.01;
H = 0.2:res:200; % 10Hz - 10kHz
Hf = H*f1;

Hi = @(s) Kp_i+Ki_i/s;
Hp = @(s) (Kp_p+Ki_p/s)/s;

ZpM = zeros(1,length(H));
ZnM = zeros(1,length(H));
for hh=1:length(H)
    h=H(hh);
    s = h*1i*w1;
    
    Tpll_subt1 = V1*Hp(s-1i*w1)/(2*(1+V1*Hp(s-1i*w1)));
    ZpM(hh) = (Hi(s-1i*w1)*V0 + (s-1i*w1)*L)/...
        (1-Tpll_subt1*(1+Hi(s-1i*w1)*I1*V0/V1));
    
    Tpll_add1 = V1*Hp(s+1i*w1)/(2*(1+V1*Hp(s+1i*w1)));
    ZnM(hh) = (Hi(s+1i*w1)*V0 + (s+1i*w1)*L)/...
        (1-Tpll_add1*(1+Hi(s+1i*w1)*I1*V0/V1));
end

% S = tf('s') specifies the transfer function H(s) = s (Laplace variable).
%      Z = tf('z',TS) specifies H(z) = z with sample time TS.
%      You can then specify transfer functions directly as expressions in S
%      or Z, for example,
%         s = tf('s');  H = exp(-s)*(s+1)/(s^2+3*s+1)
s = tf('s');
Tf_p = ((Kp_i+Ki_i/(s-1i*w1))*V0 + (s-1i*w1)*L)/...
    (1-V1*((Kp_p+Ki_p/(s-1i*w1))/(s-1i*w1))/(2*(1+V1*((Kp_p+Ki_p/(s-1i*w1))/(s-1i*w1))))*(1+(Kp_i+Ki_i/(s-1i*w1))*I1*V0/V1));
Tf_n = ((Kp_i+Ki_i/(s+1i*w1))*V0 + (s+1i*w1)*L)/...
    (1-V1*((Kp_p+Ki_p/(s+1i*w1))/(s+1i*w1))/(2*(1+V1*((Kp_p+Ki_p/(s+1i*w1))/(s+1i*w1))))*(1+(Kp_i+Ki_i/(s+1i*w1))*I1*V0/V1));
clear s
figure(1)
nyquist(Tf_p); hold on
nyquist(Tf_n)
circle(0,0,1);
% axis([-1.5 1.5 -1.5 1.5])
hold off
return
figure(1)
semilogx(Hf,mag2db(abs(ZpM))); hold on
semilogx(Hf,mag2db(abs(ZnM)), 'r'); 
hold off

figure(2)
semilogx(Hf,rad2deg(angle(ZpM))); hold on
semilogx(Hf,rad2deg(angle(ZnM)), 'r'); 
hold off



%%%%