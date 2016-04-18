clear variables;clc;close all

%% HVDC converter

Vdc = 300e3; % HVDC DC link voltage
Vm = 122e3; % ph. voltage ampl.
L = 17.9e-3; % ph. inductance

Lf = 0.43e-6; % filter imp.
Cf = 17.7e-6; % filter cap.
Rf = 10; % filter res.

% Hi(s) current control compensator
Kp_i = 0.075e-3;
Ki_i = 0.094;

Kid = 23e-6; % current decoupling coefficient
fi = 250; % current control bandwidth

% Hp(s) PLL compensator
Kp_v = 11.1e-3;
Ki_v = 8.388;

Kvd = 5.56e-3; % voltage decoupling coefficient
fv = 100; % PLL bandwidth

f1 = 50;
w1 = 2*pi*f1;
res = 0.01;
H = 0.2:res:200; % 10Hz - 10kHz
Hf = H*f1;

Hi = @(s) Kp_i+Ki_i/s;
Hv = @(s) Kp_v+Ki_v/s;

ZpM = zeros(1,length(H));
ZnM = zeros(1,length(H));
for hh=1:length(H)
    h=H(hh);
    s = h*1i*w1;
    
%     Yf = imp_parallel(s*Cf,((1/Rf)+(1/(s*Lf))));
%     Yf = 1/(imp_parallel(Rf,s*Lf)+1/(s*Cf)); (the same)
    Yf = s*Cf;
    Tp = (Hi(s-1i*w1)+1i*Kid)*Hv(s-1i*w1)*Vdc;
    ZpM(hh) = (Hi(s-1i*w1)*Vdc + s*L)/...
        (1+Yf*(Hi(s-1i*w1)*Vdc+s*L)+Tp);
    
    Tn = (Hi(s+1i*w1)-1i*Kid)*Hv(s+1i*w1)*Vdc;
    ZnM(hh) = (Hi(s+1i*w1)*Vdc + s*L)/...
        (1+Yf*(Hi(s+1i*w1)*Vdc+s*L)+Tn);
end

s = tf('s');
Tf_p = ((Kp_i+Ki_i/(s-1i*w1))*Vdc + s*L)/...
    (1+(s*Cf)*((Kp_i+Ki_i/(s-1i*w1))*Vdc+s*L)+...
    ((Kp_i+Ki_i/(s-1i*w1))+1i*Kid)*(Kp_v+Ki_v/(s-1i*w1))*Vdc);
Tf_n = ((Kp_i+Ki_i/(s+1i*w1))*Vdc + s*L)/...
    (1+(s*Cf)*((Kp_i+Ki_i/(s+1i*w1))*Vdc+s*L)+...
    ((Kp_i+Ki_i/(s+1i*w1))+1i*Kid)*(Kp_v+Ki_v/(s+1i*w1))*Vdc);
clear s

figure(1)
nyquist(Tf_p, {0.1,306}); hold on
% nyquist(Tf_n, {0.1,306})
circle(0,0,1);
axis([-1.5 1.5 -1.5 1.5])
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